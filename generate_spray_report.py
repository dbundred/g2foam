#!/opt/homebrew/opt/python@3.13/bin/python3.13
"""
Generate a comprehensive 6-page PDF report for the OpenFOAM ECN Spray G2 case,
including direct comparisons against official Sandia Engine Combustion Network (ECN)
experimental datasets (liquid & vapour penetration, spray widths at Z = 15, 25, 35 mm,
and 2D Diffused Back-Illumination / Schlieren spray envelopes).

Pages:
  1. Parcels coloured by droplet radius (um) on Gas Temperature (K) background (0.1 - 1.5 ms)
     + ECN experimental liquid penetration markers
  2. Parcels coloured by droplet temperature (K) on Gas Temperature (K) background (0.1 - 1.5 ms)
     + ECN experimental liquid penetration markers
  3. Parcels coloured by droplet radius (um) on Fuel Vapour Mass Fraction (IC8H18) background (0.1 - 1.5 ms)
     + ECN experimental Schlieren vapour penetration markers
  4. Side-by-Side Comparison at Key Spray Stages (0.2, 0.5, 0.8, 1.2 ms)
     + ECN 2D experimental DBI liquid boundary (PLV = 0.2e-3) and Schlieren 50% vapour boundary contours
  5. ECN Experimental Validation Benchmark:
     - Axial Liquid Penetration (OpenFOAM PLV & 95%/99% mass vs. Sandia ECN G2-cold & G2)
     - Axial Vapour Penetration (OpenFOAM Y_IC8H18 >= 0.1% vs. Sandia ECN Schlieren)
     - Liquid Spray Width at Z = 15, 25, 35 mm (OpenFOAM vs. Sandia ECN DBI)
     - Vapour Spray Width at Z = 15, 25, 35 mm (OpenFOAM vs. Sandia ECN Schlieren)
  6. Global Spray Physics Summary (Fuel mass balance, penetraton, Sauter Mean Radius D32/2, Gas/Droplet T)

Usage:
    /opt/homebrew/opt/python@3.13/bin/python3.13 generate_spray_report.py
"""

import glob
import os
import re
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.lines import Line2D


# =============================================================================
# 1. OpenFOAM Binary & Mesh Readers
# =============================================================================

def parse_blockmesh_info(blockmesh_path="system/blockMeshDict"):
    """Read domain bounds and grid dimensions from system/blockMeshDict."""
    with open(blockmesh_path, "r") as f:
        text = f.read()

    scale_m = re.search(r"scale\s+([0-9.eE+-]+)\s*;", text)
    scale = float(scale_m.group(1)) if scale_m else 1.0

    verts_m = re.search(r"vertices\s*\((.*?)\)\s*;", text, re.S)
    verts_raw = re.findall(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", verts_m.group(1))
    verts = np.array([[float(a), float(b), float(c)] for a, b, c in verts_raw]) * scale

    hex_m = re.search(r"hex\s*\([^)]+\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", text)
    n1, n2, n3 = int(hex_m.group(1)), int(hex_m.group(2)), int(hex_m.group(3))
    nz, nx, ny = n1, n2, n3

    xmin, xmax = verts[:, 0].min(), verts[:, 0].max()
    ymin, ymax = verts[:, 1].min(), verts[:, 1].max()
    zmin, zmax = verts[:, 2].min(), verts[:, 2].max()

    return {
        "nx": nx, "ny": ny, "nz": nz,
        "ncells": nx * ny * nz,
        "xmin": xmin, "xmax": xmax,
        "ymin": ymin, "ymax": ymax,
        "zmin": zmin, "zmax": zmax,
    }


def read_vol_scalar(path, n_cells):
    """Read an OpenFOAM binary or ASCII volScalarField internalField."""
    with open(path, "rb") as f:
        raw = f.read()

    m_uni = re.search(rb"internalField\s+uniform\s+([0-9.eE+-]+)\s*;", raw)
    if m_uni:
        return np.full(n_cells, float(m_uni.group(1)), dtype=np.float64)

    m_bin = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(", raw)
    if m_bin:
        n = int(m_bin.group(1))
        return np.frombuffer(raw, dtype="<f8", count=n, offset=m_bin.end()).copy()

    raise ValueError(f"Could not parse internalField in {path}")


def read_lagrangian_scalar(path):
    """Read an OpenFOAM binary or ASCII lagrangian scalarField."""
    if not os.path.exists(path):
        return np.array([], dtype=np.float64)
    with open(path, "rb") as f:
        raw = f.read()

    m_uni = re.search(rb"// \* \* \*.*?\n\s*(\d+)\s*\{([0-9.eE+-]+)\}", raw, re.S)
    if m_uni:
        return np.full(int(m_uni.group(1)), float(m_uni.group(2)), dtype=np.float64)

    m_bin = re.search(rb"// \* \* \*.*?\n\s*(\d+)\s*\(", raw, re.S)
    if m_bin:
        n = int(m_bin.group(1))
        if n == 0:
            return np.array([], dtype=np.float64)
        return np.frombuffer(raw, dtype="<f8", count=n, offset=m_bin.end()).copy()

    return np.array([], dtype=np.float64)


def read_lagrangian_positions(path):
    """Read OpenFOAM v2606 binary lagrangian positions file: (x y z cellI)."""
    if not os.path.exists(path):
        return np.zeros((0, 3), dtype=np.float64)
    with open(path, "rb") as f:
        raw = f.read()

    m = re.search(rb"// \* \* \*.*?\n\s*(\d+)\s*\(\n", raw, re.S)
    if not m:
        return np.zeros((0, 3), dtype=np.float64)
    n = int(m.group(1))
    if n == 0:
        return np.zeros((0, 3), dtype=np.float64)

    dt = np.dtype([
        ("lp", "S1"),
        ("x", "<f8"),
        ("y", "<f8"),
        ("z", "<f8"),
        ("cellI", "<i4"),
        ("rp", "S2"),
    ])
    recs = np.frombuffer(raw, dtype=dt, count=n, offset=m.end())
    return np.column_stack([recs["x"], recs["y"], recs["z"]])


def load_case_snapshot(time_dir, mesh_info, slice_half_width_m=0.0025, y_nozzle_mm=99.5):
    """
    Load Eulerian z=0 slice, Lagrangian parcels, and ECN-standardized penetration & width
    metrics for a given OpenFOAM time directory.
    """
    nx, ny, nz = mesh_info["nx"], mesh_info["ny"], mesh_info["nz"]
    ncells = mesh_info["ncells"]

    T_gas_3d = read_vol_scalar(os.path.join(time_dir, "T"), ncells).reshape((ny, nx, nz))
    Y_fuel_3d = read_vol_scalar(os.path.join(time_dir, "IC8H18"), ncells).reshape((ny, nx, nz))

    # Central z=0 slice: average the two cells straddling z = 0
    iz0 = nz // 2 - 1
    iz1 = nz // 2
    T_slice = 0.5 * (T_gas_3d[:, :, iz0] + T_gas_3d[:, :, iz1])
    Y_slice = 0.5 * (Y_fuel_3d[:, :, iz0] + Y_fuel_3d[:, :, iz1])

    # Load Lagrangian parcels
    cloud_dir = os.path.join(time_dir, "lagrangian", "sprayCloud")
    pos = read_lagrangian_positions(os.path.join(cloud_dir, "positions"))
    d = read_lagrangian_scalar(os.path.join(cloud_dir, "d"))
    Tp = read_lagrangian_scalar(os.path.join(cloud_dir, "T"))
    npart = read_lagrangian_scalar(os.path.join(cloud_dir, "nParticle"))

    dx_mm = (mesh_info["xmax"] - mesh_info["xmin"]) * 1000.0 / nx
    dy_mm = (mesh_info["ymax"] - mesh_info["ymin"]) * 1000.0 / ny
    x_centers = mesh_info["xmin"] * 1000.0 + (np.arange(nx) + 0.5) * dx_mm
    y_centers = mesh_info["ymin"] * 1000.0 + (np.arange(ny) + 0.5) * dy_mm
    z_axial_cells = y_nozzle_mm - y_centers  # Axial distance downstream of injector tip [mm]

    if len(pos) > 0:
        pos_mm_all = pos * 1000.0
        d_mm_all = d * 1000.0
        vol_p_mm3 = npart * (np.pi / 6.0) * (d_mm_all ** 3)
        z_axial_p = y_nozzle_mm - pos_mm_all[:, 1]

        # 1. Axial 95% and 99% liquid mass penetration
        order_z = np.argsort(z_axial_p)
        cum_z = np.cumsum(vol_p_mm3[order_z]) / np.sum(vol_p_mm3)
        z_liq_95 = float(z_axial_p[order_z[np.searchsorted(cum_z, 0.95)]])
        z_liq_99 = float(z_axial_p[order_z[np.searchsorted(cum_z, 0.99)]])

        # 2. ECN Projected Liquid Volume (PLV in mm^3 liquid / mm^2 area) projected along z line-of-sight
        xbins = np.linspace(mesh_info["xmin"] * 1000.0, mesh_info["xmax"] * 1000.0, nx + 1)
        zbins = np.linspace(0.0, (mesh_info["ymax"] - mesh_info["ymin"]) * 1000.0, ny + 1)
        H_vol, _, _ = np.histogram2d(z_axial_p, pos_mm_all[:, 0], bins=[zbins, xbins], weights=vol_p_mm3)
        PLV = H_vol / (dx_mm * dy_mm)
        z_bin_centers = 0.5 * (zbins[:-1] + zbins[1:])

        has_low = (PLV >= 0.2e-3).any(axis=1)
        has_high = (PLV >= 2.0e-3).any(axis=1)
        plv_pen_low = float(z_bin_centers[has_low].max()) if has_low.any() else 0.0
        plv_pen_high = float(z_bin_centers[has_high].max()) if has_high.any() else 0.0

        # Liquid spray width at Z = 15, 25, 35 mm from PLV >= 0.2e-3
        w_liq = {}
        for zt in [15.0, 25.0, 35.0]:
            iz_b = np.argmin(np.abs(z_bin_centers - zt))
            xs = x_centers[PLV[iz_b, :] >= 0.2e-3]
            w_liq[zt] = float(xs.max() - xs.min() + dx_mm) if len(xs) > 0 else 0.0

        # Filter slice parcels
        mask = np.abs(pos[:, 2]) <= slice_half_width_m
        pos_slice = pos_mm_all[mask]
        r_um = 0.5 * d[mask] * 1e6
        Tp_slice = Tp[mask]
    else:
        pos_slice = np.zeros((0, 3))
        r_um = np.array([])
        Tp_slice = np.array([])
        z_liq_95, z_liq_99 = 0.0, 0.0
        plv_pen_low, plv_pen_high = 0.0, 0.0
        w_liq = {15.0: 0.0, 25.0: 0.0, 35.0: 0.0}

    # 3. ECN Vapour Penetration & Width (where projected fuel mass fraction Y_IC8H18 >= 0.001 = 0.1%)
    Y_proj = Y_fuel_3d.max(axis=2)  # Line-of-sight maximum across z
    has_vap = (Y_proj >= 0.001).any(axis=1)
    vap_pen = float(z_axial_cells[has_vap].max()) if has_vap.any() else 0.0

    w_vap = {}
    for zt in [15.0, 25.0, 35.0]:
        iy_c = np.argmin(np.abs(z_axial_cells - zt))
        xs = x_centers[Y_proj[iy_c, :] >= 0.001]
        w_vap[zt] = float(xs.max() - xs.min() + dx_mm) if len(xs) > 0 else 0.0

    return {
        "time_s": float(time_dir),
        "time_ms": float(time_dir) * 1000.0,
        "T_slice": T_slice,
        "Y_slice": Y_slice,
        "pos_mm": pos_slice,
        "r_um": r_um,
        "Tp": Tp_slice,
        "total_parcels": len(pos),
        "slice_parcels": len(pos_slice),
        "z_liq_95": z_liq_95,
        "z_liq_99": z_liq_99,
        "plv_pen_low": plv_pen_low,
        "plv_pen_high": plv_pen_high,
        "vap_pen": vap_pen,
        "w_liq": w_liq,
        "w_vap": w_vap,
    }


# =============================================================================
# 2. ECN Sandia Experimental Data Loader
# =============================================================================

def load_ecn_experimental_data(exp_dir="ECN_exp_data"):
    """
    Load Sandia ECN Spray G2-cold (T_amb = 333 K, p_amb = 50 kPa) and G2 (T_amb = 573 K)
    experimental datasets downloaded from https://ecn.sandia.gov/data/sandia-spray-g-data/.
    """
    data = {}
    for key, prefix in [("G2C", "ic8_G2C"), ("G2", "ic8_G2")]:
        liq_low_path = os.path.join(exp_dir, f"{prefix}_liquid_length", "liquid_length", "liquid_data_2e-4_thr.mat")
        liq_high_path = os.path.join(exp_dir, f"{prefix}_liquid_length", "liquid_length", "liquid_data_2e-3_thr.mat")
        vap_path = os.path.join(exp_dir, f"{prefix}_vapor_length", "vapor_length", "vapor_data.mat")
        b_low_path = os.path.join(exp_dir, f"{prefix}_liquid_length", "liquid_length", "liquid_boundary_2e-4_womax.mat")
        sch_path = os.path.join(exp_dir, f"{prefix}_vapor_length", "vapor_length", "Schlieren_probability.mat")

        if not os.path.exists(liq_low_path):
            continue

        l_low = sio.loadmat(liq_low_path)["liquid_data_low_thr"]
        l_high = sio.loadmat(liq_high_path)["liquid_data_high_thr"]
        vap = sio.loadmat(vap_path)["vapor_data"]

        # Mask valid window (0 <= t <= 1.6 ms)
        m_l = (l_low[:, 0] >= 0.0) & (l_low[:, 0] <= 1.6)
        m_v = (vap[:, 0] >= 0.0) & (vap[:, 0] <= 1.6)

        entry = {
            "t_liq": l_low[m_l, 0],
            # Column 16 (index 15): ensemble-averaged liquid penetration length [mm], Column 17 (index 16): std dev [mm]
            "pen_liq_low": np.where(l_low[m_l, 15] > 0.5, l_low[m_l, 15], np.nan),
            "pen_liq_low_std": l_low[m_l, 16],
            "pen_liq_high": np.where(l_high[m_l, 15] > 0.5, l_high[m_l, 15], np.nan),
            "pen_liq_high_std": l_high[m_l, 16],
            # Widths at Z = 15, 25, 35 mm (cols 2, 6, 10 -> indices 1, 5, 9; std at indices 2, 6, 10)
            "w15_liq": np.where(l_low[m_l, 1] > 0.5, l_low[m_l, 1], np.nan),
            "w15_liq_std": l_low[m_l, 2],
            "w25_liq": np.where(l_low[m_l, 5] > 0.5, l_low[m_l, 5], np.nan),
            "w25_liq_std": l_low[m_l, 6],
            "w35_liq": np.where(l_low[m_l, 9] > 0.5, l_low[m_l, 9], np.nan),
            "w35_liq_std": l_low[m_l, 10],
            # Vapour data: col 1 (index 0) time [ms], col 12 (index 11) penetration [mm], col 13 (index 12) std [mm]
            "t_vap": vap[m_v, 0],
            "pen_vap": np.where(vap[m_v, 11] > 0.5, vap[m_v, 11], np.nan),
            "pen_vap_std": vap[m_v, 12],
            # Vapour widths at Z = 15, 25, 35 mm (cols 2, 5, 8 -> indices 1, 4, 7; std at indices 2, 5, 8)
            "w15_vap": np.where(vap[m_v, 1] > 0.5, vap[m_v, 1], np.nan),
            "w15_vap_std": vap[m_v, 2],
            "w25_vap": np.where(vap[m_v, 4] > 0.5, vap[m_v, 4], np.nan),
            "w25_vap_std": vap[m_v, 5],
            "w35_vap": np.where(vap[m_v, 7] > 0.5, vap[m_v, 7], np.nan),
            "w35_vap_std": vap[m_v, 8],
            "raw_t_liq": l_low[:, 0],
            "raw_t_vap": vap[:, 0],
        }

        if key == "G2C":
            if os.path.exists(b_low_path):
                entry["liq_boundary_imgs"] = sio.loadmat(b_low_path, squeeze_me=True)["outline1_wo_outline"]
            if os.path.exists(sch_path):
                entry["schlieren_imgs"] = sio.loadmat(sch_path, squeeze_me=True)["Schlieren_process"]

        data[key] = entry
    return data


def parse_sprayfoam_log(log_path="log.sprayFoam"):
    """Extract high-frequency time-history metrics from log.sprayFoam."""
    if not os.path.exists(log_path):
        return None
    with open(log_path, "r") as f:
        text = f.read()

    times_ms, liq_mass_mg, evap_mass_mg, pen_mm, d32_um, t_gas_min, t_gas_max = [], [], [], [], [], [], []
    blocks = re.split(r"\nTime = ([0-9.eE+-]+)\n", text)
    for i in range(1, len(blocks) - 1, 2):
        t_val = float(blocks[i]) * 1000.0
        b = blocks[i + 1]
        m_mass = re.search(r"Current mass in system\s*=\s*([0-9.eE+-]+)", b)
        m_evap = re.search(r"Mass transfer phase change\s*=\s*([0-9.eE+-]+)", b)
        m_d32 = re.search(r"D10, D32, Dmax \(mu\)\s*=\s*[0-9.eE+-]+,\s*([0-9.eE+-]+),", b)
        m_pen = re.search(r"Liquid penetration 95% mass \(m\)\s*=\s*([0-9.eE+-]+)", b)
        m_tgas = re.search(r"T gas min/max\s+([0-9.eE+-]+),\s*([0-9.eE+-]+)", b)

        if m_mass and m_evap and m_d32 and m_pen and m_tgas:
            times_ms.append(t_val)
            liq_mass_mg.append(float(m_mass.group(1)) * 1e6)
            evap_mass_mg.append(float(m_evap.group(1)) * 1e6)
            d32_um.append(float(m_d32.group(1)))
            pen_mm.append(float(m_pen.group(1)) * 1000.0)
            t_gas_min.append(float(m_tgas.group(1)))
            t_gas_max.append(float(m_tgas.group(2)))

    return {
        "t_ms": np.array(times_ms),
        "liq_mass_mg": np.array(liq_mass_mg),
        "evap_mass_mg": np.array(evap_mass_mg),
        "d32_um": np.array(d32_um),
        "pen_mm": np.array(pen_mm),
        "t_gas_min": np.array(t_gas_min),
        "t_gas_max": np.array(t_gas_max),
    }


# =============================================================================
# 3. Plotting Pages
# =============================================================================

def plot_time_series_grid(
    pdf,
    snapshots,
    mesh_info,
    ecn_data,
    bg_key,
    bg_cmap,
    bg_norm,
    bg_label,
    parcel_key,
    parcel_cmap,
    parcel_norm,
    parcel_label,
    page_title,
    page_subtitle,
    ref_type="liquid",
    xlim=(-55, 55),
    ylim=(15, 102),
    y_nozzle_mm=99.5,
):
    """Render a 3x5 grid of all 15 time steps with ECN experimental penetration reference lines."""
    fig, axes = plt.subplots(3, 5, figsize=(16.5, 10.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.055, right=0.88, bottom=0.07, top=0.88, wspace=0.08, hspace=0.16)

    fig.suptitle(page_title, fontsize=14.5, fontweight="bold", y=0.96, x=0.47)
    fig.text(0.47, 0.925, page_subtitle, ha="center", fontsize=10, color="#333333")

    extent = [
        mesh_info["xmin"] * 1000.0,
        mesh_info["xmax"] * 1000.0,
        mesh_info["ymin"] * 1000.0,
        mesh_info["ymax"] * 1000.0,
    ]

    g2c = ecn_data.get("G2C")

    for idx, snap in enumerate(snapshots):
        row, col = divmod(idx, 5)
        ax = axes[row, col]

        ax.imshow(
            snap[bg_key],
            origin="lower",
            extent=extent,
            cmap=bg_cmap,
            norm=bg_norm,
            interpolation="bilinear",
            aspect="equal",
        )

        ax.axvline(0.0, color="white", linestyle="--", linewidth=0.45, alpha=0.35)
        ax.plot(0.0, y_nozzle_mm, marker="v", color="white", markeredgecolor="black", markersize=5, zorder=7)

        # Overlay Lagrangian parcels in slice
        pos = snap["pos_mm"]
        if len(pos) > 0:
            p_val = snap[parcel_key]
            sizes = np.clip(1.8 + 0.32 * snap["r_um"], 2.0, 22.0)
            ax.scatter(
                pos[:, 0],
                pos[:, 1],
                c=p_val,
                cmap=parcel_cmap,
                norm=parcel_norm,
                s=sizes,
                alpha=0.88,
                edgecolors="black",
                linewidths=0.18,
                zorder=5,
            )

        # Overlay ECN experimental penetration reference line
        if g2c is not None:
            t_ms = snap["time_ms"]
            if ref_type == "liquid":
                pen_exp = np.interp(t_ms, g2c["t_liq"], g2c["pen_liq_low"])
                if not np.isnan(pen_exp) and pen_exp > 1.0:
                    y_exp = y_nozzle_mm - pen_exp
                    ax.axhline(y_exp, color="#00ffff", linestyle="--", linewidth=1.2, alpha=0.9, zorder=6)
                    ax.text(xlim[0] + 3, y_exp + 1.5, f"ECN Liq: {pen_exp:.0f}mm", color="#00ffff",
                            fontsize=6.8, fontweight="bold", bbox=dict(boxstyle="square,pad=0.15", facecolor="black", alpha=0.5, edgecolor="none"))
            else:
                pen_exp = np.interp(t_ms, g2c["t_vap"], g2c["pen_vap"])
                if not np.isnan(pen_exp) and pen_exp > 1.0:
                    y_exp = y_nozzle_mm - pen_exp
                    ax.axhline(y_exp, color="#ff9900", linestyle="--", linewidth=1.2, alpha=0.95, zorder=6)
                    ax.text(xlim[0] + 3, y_exp + 1.5, f"ECN Vap: {pen_exp:.0f}mm", color="#ff9900",
                            fontsize=6.8, fontweight="bold", bbox=dict(boxstyle="square,pad=0.15", facecolor="black", alpha=0.5, edgecolor="none"))

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        inj_status = "INJ ON" if snap["time_ms"] <= 0.78 else "POST-INJ"
        ax.set_title(
            f"t = {snap['time_ms']:.1f} ms ({inj_status}) | n={snap['slice_parcels']}",
            fontsize=8.6,
            fontweight="semibold",
            pad=4,
        )

        if row == 2:
            ax.set_xlabel("x [mm]", fontsize=9)
        if col == 0:
            ax.set_ylabel("y [mm] (Injector at 99.5)", fontsize=9)
        ax.tick_params(labelsize=8)

    cbar_ax1 = fig.add_axes([0.895, 0.50, 0.016, 0.36])
    cbar1 = fig.colorbar(cm.ScalarMappable(norm=bg_norm, cmap=bg_cmap), cax=cbar_ax1)
    cbar1.set_label(bg_label, fontsize=10, fontweight="semibold")
    cbar1.ax.tick_params(labelsize=8.5)

    cbar_ax2 = fig.add_axes([0.895, 0.08, 0.016, 0.36])
    cbar2 = fig.colorbar(cm.ScalarMappable(norm=parcel_norm, cmap=parcel_cmap), cax=cbar_ax2)
    cbar2.set_label(parcel_label, fontsize=10, fontweight="semibold")
    cbar2.ax.tick_params(labelsize=8.5)

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def plot_comparative_stages_page(pdf, snapshots, mesh_info, ecn_data, bg_norm_T, bg_norm_Y, p_norm_r, p_norm_T, y_nozzle_mm=99.5):
    """
    Render a 3x4 side-by-side comparison page at 4 key physical stages (0.2, 0.5, 0.8, 1.2 ms)
    with direct 2D Sandia ECN experimental DBI liquid boundary & Schlieren vapour boundary overlays!
    """
    target_times = [0.2, 0.5, 0.8, 1.2]
    selected = [s for s in snapshots if any(abs(s["time_ms"] - tt) < 0.01 for tt in target_times)]

    fig, axes = plt.subplots(3, 4, figsize=(16.5, 10.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.88, bottom=0.08, top=0.86, wspace=0.08, hspace=0.15)

    fig.suptitle(
        "4. Side-by-Side Stages + Sandia ECN 2D Experimental Envelope Overlays (z = 0 Slice)",
        fontsize=14.5,
        fontweight="bold",
        y=0.96,
        x=0.47,
    )
    fig.text(
        0.47,
        0.922,
        "Cyan Contour: Sandia ECN DBI Liquid Envelope (PLV = 0.2e-3 mm³/mm²)  |  Orange Contour: Sandia ECN Schlieren 50% Vapour Envelope",
        ha="center",
        fontsize=10.2,
        color="#222222",
    )

    extent = [
        mesh_info["xmin"] * 1000.0,
        mesh_info["xmax"] * 1000.0,
        mesh_info["ymin"] * 1000.0,
        mesh_info["ymax"] * 1000.0,
    ]

    g2c = ecn_data.get("G2C")

    row_configs = [
        ("T_slice", "inferno", bg_norm_T, "r_um", "turbo", p_norm_r, "Row 1: Radius on Gas T"),
        ("T_slice", "inferno", bg_norm_T, "Tp", "coolwarm", p_norm_T, "Row 2: Droplet T on Gas T"),
        ("Y_slice", "viridis", bg_norm_Y, "r_um", "plasma", p_norm_r, "Row 3: Radius on Y_IC8H18"),
    ]

    for r_idx, (bg_k, bg_c, bg_n, p_k, p_c, p_n, row_title) in enumerate(row_configs):
        for c_idx, snap in enumerate(selected):
            ax = axes[r_idx, c_idx]
            ax.imshow(
                snap[bg_k],
                origin="lower",
                extent=extent,
                cmap=bg_c,
                norm=bg_n,
                interpolation="bilinear",
                aspect="equal",
            )
            ax.axvline(0.0, color="white", linestyle="--", linewidth=0.5, alpha=0.35)
            ax.plot(0.0, y_nozzle_mm, marker="v", color="white", markeredgecolor="black", markersize=5, zorder=7)

            pos = snap["pos_mm"]
            if len(pos) > 0:
                sizes = np.clip(2.2 + 0.38 * snap["r_um"], 2.5, 26.0)
                ax.scatter(
                    pos[:, 0],
                    pos[:, 1],
                    c=snap[p_k],
                    cmap=p_c,
                    norm=p_n,
                    s=sizes,
                    alpha=0.9,
                    edgecolors="black",
                    linewidths=0.2,
                    zorder=5,
                )

            # Overlay 2D ECN experimental boundary contour from the 512x512 camera frames
            if g2c is not None:
                t_ms = snap["time_ms"]
                if r_idx in (0, 1) and "liq_boundary_imgs" in g2c:
                    idx_f = int(np.argmin(np.abs(g2c["raw_t_liq"] - t_ms)))
                    b_img = g2c["liq_boundary_imgs"][idx_f]
                    # Pixel calibration: row 260 is x=0, col 48 is Z_axial=0, scale = 0.1908 mm/px
                    rows_px = np.arange(b_img.shape[0])
                    cols_px = np.arange(b_img.shape[1])
                    X_cam = (rows_px - 260.0) * 0.1908
                    Y_cam = y_nozzle_mm - (cols_px - 48.0) * 0.1908
                    Yg, Xg = np.meshgrid(Y_cam, X_cam)
                    ax.contour(Xg, Yg, b_img, levels=[0.5], colors=["#00ffff"], linewidths=1.4, linestyles="--", zorder=6)
                elif r_idx == 2 and "schlieren_imgs" in g2c:
                    idx_v = int(np.argmin(np.abs(g2c["raw_t_vap"] - t_ms)))
                    s_img = g2c["schlieren_imgs"][idx_v]
                    # Schlieren pixel calibration: row 251 is x=0, col 40 is Z_axial=0, scale = 0.218 mm/px
                    rows_px = np.arange(s_img.shape[0])
                    cols_px = np.arange(s_img.shape[1])
                    X_cam = (rows_px - 251.5) * 0.218
                    Y_cam = y_nozzle_mm - (cols_px - 40.0) * 0.218
                    Yg, Xg = np.meshgrid(Y_cam, X_cam)
                    ax.contour(Xg, Yg, s_img, levels=[50.0], colors=["#ff9900"], linewidths=1.5, linestyles="--", zorder=6)
                    if "liq_boundary_imgs" in g2c:
                        idx_f = int(np.argmin(np.abs(g2c["raw_t_liq"] - t_ms)))
                        b_img = g2c["liq_boundary_imgs"][idx_f]
                        X_l = (np.arange(b_img.shape[0]) - 260.0) * 0.1908
                        Y_l = y_nozzle_mm - (np.arange(b_img.shape[1]) - 48.0) * 0.1908
                        Ylg, Xlg = np.meshgrid(Y_l, X_l)
                        ax.contour(Xlg, Ylg, b_img, levels=[0.5], colors=["#00ffff"], linewidths=1.1, linestyles=":", zorder=6)

            ax.set_xlim(-55, 55)
            ax.set_ylim(15, 102)
            if r_idx == 0:
                ax.set_title(f"t = {snap['time_ms']:.1f} ms (n={snap['slice_parcels']} in slice)", fontsize=10.5, fontweight="bold")
            if r_idx == 2:
                ax.set_xlabel("x [mm]", fontsize=9.5)
            if c_idx == 0:
                ax.set_ylabel(f"{row_title}\ny [mm]", fontsize=9, fontweight="semibold")

    cb_specs = [
        ([0.895, 0.68, 0.014, 0.18], bg_norm_T, "inferno", "Gas T [K]"),
        ([0.895, 0.46, 0.014, 0.18], p_norm_T, "coolwarm", "Droplet T [K]"),
        ([0.895, 0.26, 0.014, 0.16], bg_norm_Y, "viridis", "Y_IC8H18 [-]"),
        ([0.895, 0.07, 0.014, 0.16], p_norm_r, "turbo", "Radius [um]"),
    ]
    for pos_rect, norm_obj, cmap_name, label_str in cb_specs:
        cax = fig.add_axes(pos_rect)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm_obj, cmap=cmap_name), cax=cax)
        cb.set_label(label_str, fontsize=9, fontweight="semibold")
        cb.ax.tick_params(labelsize=8)

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def plot_ecn_validation_page(pdf, snapshots, ecn_data, log_data):
    """
    Page 5: Comprehensive 2x2 quantitative validation benchmark against Sandia ECN datasets:
      - Top Left: Axial Liquid Penetration vs. ECN G2-cold (333 K) & G2 (573 K)
      - Top Right: Axial Vapour Penetration vs. ECN G2-cold & G2 Schlieren
      - Bottom Left: Liquid Spray Width at Z = 15, 25, 35 mm vs. ECN DBI
      - Bottom Right: Vapour Spray Width at Z = 15, 25, 35 mm vs. ECN Schlieren
    """
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 10.0))
    fig.subplots_adjust(left=0.065, right=0.96, bottom=0.075, top=0.84, wspace=0.20, hspace=0.28)

    fig.suptitle(
        "5. Quantitative Validation against Sandia ECN Spray G2 Experimental Database",
        fontsize=15,
        fontweight="bold",
        y=0.955,
    )
    fig.text(
        0.5,
        0.902,
        "Reference Data: Sandia National Laboratories Engine Combustion Network (ECN) AV67-028 Iso-Octane Dataset (Hwang & Pickett)\n"
        "Comparing OpenFOAM (T_amb = 333 K, p_amb = 50 kPa) against ECN G2-cold (T_amb = 333 K, p_amb = 50 kPa) and ECN G2 Hot (T_amb = 573 K)",
        ha="center",
        fontsize=9.8,
        color="#1e293b",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )

    t_sim = np.array([s["time_ms"] for s in snapshots])
    z_liq_95 = np.array([s["z_liq_95"] for s in snapshots])
    z_liq_99 = np.array([s["z_liq_99"] for s in snapshots])
    plv_low = np.array([s["plv_pen_low"] for s in snapshots])
    plv_high = np.array([s["plv_pen_high"] for s in snapshots])
    vap_sim = np.array([s["vap_pen"] for s in snapshots])

    g2c = ecn_data.get("G2C")
    g2 = ecn_data.get("G2")

    # --- Panel 1: Axial Liquid Penetration ---
    ax = axes[0, 0]
    if g2c is not None:
        ax.plot(g2c["t_liq"], g2c["pen_liq_low"], "-", color="#0284c7", linewidth=2.2,
                label="ECN Exp G2-cold (333 K) PLV >= 0.2e-3")
        ax.fill_between(
            g2c["t_liq"],
            g2c["pen_liq_low"] - g2c["pen_liq_low_std"],
            g2c["pen_liq_low"] + g2c["pen_liq_low_std"],
            color="#0284c7",
            alpha=0.18,
            label="ECN G2-cold ±1σ Band (300 injections)",
        )
        ax.plot(g2c["t_liq"], g2c["pen_liq_high"], "--", color="#0369a1", linewidth=1.6,
                label="ECN Exp G2-cold (333 K) PLV >= 2.0e-3")
    if g2 is not None:
        ax.plot(g2["t_liq"], g2["pen_liq_low"], ":", color="#64748b", linewidth=1.8,
                label="ECN Exp G2 Hot (573 K) PLV >= 0.2e-3")

    ax.plot(t_sim, z_liq_95, "o-", color="#dc2626", linewidth=2.2, markersize=5,
            label="OpenFOAM Axial 95% Liquid Mass Penetration")
    ax.plot(t_sim, z_liq_99, "s--", color="#ea580c", linewidth=1.8, markersize=4.5,
            label="OpenFOAM Axial 99% Liquid Mass Penetration")
    ax.plot(t_sim, plv_low, "^:", color="#b91c1c", linewidth=1.5, markersize=4.5,
            label="OpenFOAM Projected PLV >= 0.2e-3 (2 mm grid)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Axial Liquid Penetration Length vs. ECN Experiments", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Axial Liquid Penetration Z [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 95)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.0, loc="upper left")

    # --- Panel 2: Axial Vapour Penetration ---
    ax = axes[0, 1]
    if g2c is not None:
        ax.plot(g2c["t_vap"], g2c["pen_vap"], "-", color="#16a34a", linewidth=2.3,
                label="ECN Exp G2-cold (333 K) Schlieren Vapour")
        ax.fill_between(
            g2c["t_vap"],
            g2c["pen_vap"] - g2c["pen_vap_std"],
            g2c["pen_vap"] + g2c["pen_vap_std"],
            color="#16a34a",
            alpha=0.18,
            label="ECN G2-cold Schlieren ±1σ Band",
        )
    if g2 is not None:
        ax.plot(g2["t_vap"], g2["pen_vap"], ":", color="#64748b", linewidth=1.8,
                label="ECN Exp G2 Hot (573 K) Schlieren Vapour")

    ax.plot(t_sim, vap_sim, "o-", color="#7c3aed", linewidth=2.2, markersize=5.5,
            label="OpenFOAM Vapour Penetration (Y_IC8H18 >= 0.1%)")
    if log_data is not None:
        ax.plot(log_data["t_ms"], log_data["pen_mm"], "--", color="#94a3b8", linewidth=1.4,
                label="OpenFOAM 3D Radial 95% Parcel Distance (|x - x0|)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Axial Vapour Penetration Length vs. ECN Schlieren", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Axial Penetration Z [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 102)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="upper left")

    # --- Panel 3: Liquid Spray Width at Z = 15, 25, 35 mm ---
    ax = axes[1, 0]
    colors_z = {15.0: "#2563eb", 25.0: "#d97706", 35.0: "#16a34a"}
    if g2c is not None:
        for zt, k_w, k_s in [(15.0, "w15_liq", "w15_liq_std"), (25.0, "w25_liq", "w25_liq_std"), (35.0, "w35_liq", "w35_liq_std")]:
            c = colors_z[zt]
            ax.plot(g2c["t_liq"], g2c[k_w], "-", color=c, linewidth=2.0, label=f"ECN Exp Liquid Width @ Z={int(zt)} mm")
            ax.fill_between(g2c["t_liq"], g2c[k_w] - g2c[k_s], g2c[k_w] + g2c[k_s], color=c, alpha=0.14)

    for zt in [15.0, 25.0, 35.0]:
        c = colors_z[zt]
        w_sim = np.array([s["w_liq"][zt] for s in snapshots], dtype=float)
        w_sim[w_sim == 0] = np.nan
        ax.plot(t_sim, w_sim, "o--", color=c, linewidth=1.8, markersize=5, markeredgecolor="black",
                label=f"OpenFOAM Liquid Width @ Z={int(zt)} mm")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Liquid Spray Width at Axial Planes Z = 15, 25, 35 mm", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Spray Plume Outer Width [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8, ncol=2, loc="upper left")

    # --- Panel 4: Vapour Spray Width at Z = 15, 25, 35 mm ---
    ax = axes[1, 1]
    if g2c is not None:
        for zt, k_w, k_s in [(15.0, "w15_vap", "w15_vap_std"), (25.0, "w25_vap", "w25_vap_std"), (35.0, "w35_vap", "w35_vap_std")]:
            c = colors_z[zt]
            ax.plot(g2c["t_vap"], g2c[k_w], "-", color=c, linewidth=2.0, label=f"ECN Exp Vapour Width @ Z={int(zt)} mm")
            ax.fill_between(g2c["t_vap"], g2c[k_w] - g2c[k_s], g2c[k_w] + g2c[k_s], color=c, alpha=0.14)

    for zt in [15.0, 25.0, 35.0]:
        c = colors_z[zt]
        w_sim = np.array([s["w_vap"][zt] for s in snapshots], dtype=float)
        w_sim[w_sim == 0] = np.nan
        ax.plot(t_sim, w_sim, "s--", color=c, linewidth=1.8, markersize=5, markeredgecolor="black",
                label=f"OpenFOAM Vapour Width @ Z={int(zt)} mm")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Vapour Spray Width at Axial Planes Z = 15, 25, 35 mm", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Vapour Plume Outer Width [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 110)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8, ncol=2, loc="upper left")

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def plot_summary_metrics_page(pdf, log_data, snapshots, mesh_info, ecn_data):
    """Page 6: Global spray physics summary with ECN penetration reference curves overlaid."""
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 9.5))
    fig.subplots_adjust(left=0.07, right=0.96, bottom=0.08, top=0.84, wspace=0.22, hspace=0.28)

    fig.suptitle(
        "6. ECN Spray G2 — Global Thermodynamic, Breakup & Penetration Summary",
        fontsize=15,
        fontweight="bold",
        y=0.95,
    )
    summary_text = (
        f"Domain: {(mesh_info['xmax']-mesh_info['xmin'])*1000:.0f} x {(mesh_info['ymax']-mesh_info['ymin'])*1000:.0f} x "
        f"{(mesh_info['zmax']-mesh_info['zmin'])*1000:.0f} mm  |  "
        f"Grid: {mesh_info['nx']} x {mesh_info['ny']} x {mesh_info['nz']} ({mesh_info['ncells']:,} cells, 2 mm uniform)  |  "
        f"Ambient: p = 50 kPa, T = 333 K (N2)  |  Fuel: IC8H18 @ T0 = 363 K, 10 mg total (8 holes, 0.78 ms SOI->EOI)"
    )
    fig.text(0.5, 0.90, summary_text, ha="center", fontsize=10, color="#222222",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f4f8", edgecolor="#cbd5e1"))

    t_ms = log_data["t_ms"]
    snap_t = np.array([s["time_ms"] for s in snapshots])

    # 1. Liquid vs Evaporated Fuel Mass
    ax = axes[0, 0]
    ax.plot(t_ms, log_data["liq_mass_mg"], label="Active Liquid Mass in Chamber", color="#0284c7", linewidth=2.2)
    ax.plot(t_ms, log_data["evap_mass_mg"], label="Evaporated Fuel Vapour Mass", color="#ea580c", linewidth=2.2)
    ax.plot(t_ms, log_data["liq_mass_mg"] + log_data["evap_mass_mg"], "--", label="Total Injected Fuel Mass", color="#334155", linewidth=1.6)
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Fuel Mass Balance (Liquid vs. Vapour)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Mass [mg]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.8)

    # 2. Penetration Length with ECN Reference Curves
    ax = axes[0, 1]
    ax.plot(t_ms, log_data["pen_mm"], "--", color="#64748b", linewidth=1.6, label="OpenFOAM 3D Radial 95% Mass Dist (|x - x0|)")
    ax.plot(snap_t, [s["z_liq_95"] for s in snapshots], "o-", color="#dc2626", linewidth=2.0, markersize=4.5,
            label="OpenFOAM Axial 95% Liquid Mass Pen (Z)")
    ax.plot(snap_t, [s["vap_pen"] for s in snapshots], "s-", color="#7c3aed", linewidth=2.0, markersize=4.5,
            label="OpenFOAM Axial Vapour Pen (Y >= 0.1%)")
    g2c = ecn_data.get("G2C")
    if g2c is not None:
        ax.plot(g2c["t_liq"], g2c["pen_liq_low"], "-", color="#0284c7", linewidth=2.2,
                label="Sandia ECN G2-cold Liquid Pen (PLV 0.2e-3)")
        ax.fill_between(g2c["t_liq"], g2c["pen_liq_low"] - g2c["pen_liq_low_std"],
                        g2c["pen_liq_low"] + g2c["pen_liq_low_std"], color="#0284c7", alpha=0.15)
        ax.plot(g2c["t_vap"], g2c["pen_vap"], "-", color="#16a34a", linewidth=2.0,
                label="Sandia ECN G2-cold Schlieren Vapour Pen")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Spray Penetration Length (OpenFOAM vs. Sandia ECN)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Penetration Distance from Nozzle [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 102)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.0, loc="upper left")

    # 3. Sauter Mean Diameter (SMD / D32) & Mean Slice Radius
    ax = axes[1, 0]
    ax.plot(t_ms, log_data["d32_um"] * 0.5, color="#7c3aed", linewidth=2.2, label="Cloud Sauter Mean Radius (D32 / 2)")
    snap_r_mean = [s["r_um"].mean() if len(s["r_um"]) else 0 for s in snapshots]
    ax.plot(snap_t, snap_r_mean, "o-", color="#db2777", markersize=4.5, linewidth=1.5, label="Slice Mean Parcel Radius (r_p)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Droplet Size Evolution (Flash-Boiling & KH-RT Breakup)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Droplet Radius [um]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    # 4. Gas & Droplet Temperature Extremes
    ax = axes[1, 1]
    ax.plot(t_ms, log_data["t_gas_min"], color="#2563eb", linewidth=2.0, label="Min Gas Temperature (Evaporative Cooling)")
    ax.plot(t_ms, log_data["t_gas_max"], color="#dc2626", linewidth=2.0, label="Max Gas Temperature")
    snap_tp_min = [s["Tp"].min() if len(s["Tp"]) else np.nan for s in snapshots]
    snap_tp_mean = [s["Tp"].mean() if len(s["Tp"]) else np.nan for s in snapshots]
    ax.plot(snap_t, snap_tp_mean, "s--", color="#d97706", markersize=4, linewidth=1.5, label="Mean Droplet Temperature (Slice)")
    ax.plot(snap_t, snap_tp_min, "^:", color="#0891b2", markersize=4, linewidth=1.5, label="Min Droplet Temperature (Flash-Cooled)")
    ax.axhline(300.6, color="#0f766e", linestyle="-.", linewidth=1.1, label="IC8H18 Boiling Point @ 50 kPa (300.6 K)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Gas & Droplet Temperature Evolution", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Temperature [K]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.0, loc="best")

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


# =============================================================================
# Main Execution
# =============================================================================

def main():
    mesh_info = parse_blockmesh_info("system/blockMeshDict")
    time_dirs = sorted(
        [d for d in glob.glob("0.0*") if os.path.isdir(d) and d != "0.orig"],
        key=lambda x: float(x),
    )
    if not time_dirs:
        raise RuntimeError("No output time directories (0.0001 .. 0.0015) found!")

    slice_half_width_m = 0.0025
    print(f"Loading {len(time_dirs)} OpenFOAM time steps on {mesh_info['nx']}x{mesh_info['ny']}x{mesh_info['nz']} mesh...")
    snapshots = [load_case_snapshot(td, mesh_info, slice_half_width_m=slice_half_width_m) for td in time_dirs]
    log_data = parse_sprayfoam_log("log.sprayFoam")
    ecn_data = load_ecn_experimental_data("ECN_exp_data")

    bg_norm_T = Normalize(vmin=326.0, vmax=340.0)
    bg_norm_Y = Normalize(vmin=0.0, vmax=0.45)
    p_norm_r = Normalize(vmin=2.0, vmax=45.0)
    p_norm_T = Normalize(vmin=295.0, vmax=348.0)

    out_pdf = "ECN_SprayG2_Report.pdf"

    with PdfPages(out_pdf) as pdf:
        print("Rendering Page 1: Parcels (Radius) on Gas Temperature + ECN Liquid Pen...")
        plot_time_series_grid(
            pdf=pdf,
            snapshots=snapshots,
            mesh_info=mesh_info,
            ecn_data=ecn_data,
            bg_key="T_slice",
            bg_cmap="inferno",
            bg_norm=bg_norm_T,
            bg_label="Gas Temperature T [K]",
            parcel_key="r_um",
            parcel_cmap="turbo",
            parcel_norm=p_norm_r,
            parcel_label="Parcel Droplet Radius r [um]",
            page_title="1. Parcels Coloured by Radius on Gas Temperature Background (z = 0 Slice)",
            page_subtitle="Cyan dashed line: Sandia ECN G2-cold Experimental Liquid Penetration (PLV >= 0.2e-3 mm³/mm²) | Injector at (0, 99.5) mm",
            ref_type="liquid",
        )

        print("Rendering Page 2: Parcels (Temperature) on Gas Temperature + ECN Liquid Pen...")
        plot_time_series_grid(
            pdf=pdf,
            snapshots=snapshots,
            mesh_info=mesh_info,
            ecn_data=ecn_data,
            bg_key="T_slice",
            bg_cmap="inferno",
            bg_norm=bg_norm_T,
            bg_label="Gas Temperature T [K]",
            parcel_key="Tp",
            parcel_cmap="coolwarm",
            parcel_norm=p_norm_T,
            parcel_label="Parcel Droplet Temperature T_p [K]",
            page_title="2. Parcels Coloured by Droplet Temperature on Gas Temperature Background (z = 0 Slice)",
            page_subtitle="Cyan dashed line: Sandia ECN G2-cold Experimental Liquid Penetration | Flash-boiling cools droplets toward T_sat ≈ 300.6 K",
            ref_type="liquid",
        )

        print("Rendering Page 3: Parcels on Fuel Vapour Mass Fraction + ECN Vapour Pen...")
        plot_time_series_grid(
            pdf=pdf,
            snapshots=snapshots,
            mesh_info=mesh_info,
            ecn_data=ecn_data,
            bg_key="Y_slice",
            bg_cmap="viridis",
            bg_norm=bg_norm_Y,
            bg_label="Iso-Octane Vapour Mass Fraction Y_IC8H18 [-]",
            parcel_key="r_um",
            parcel_cmap="plasma",
            parcel_norm=p_norm_r,
            parcel_label="Parcel Droplet Radius r [um]",
            page_title="3. Parcels Coloured by Radius on Fuel Vapour Mass Fraction (IC8H18) Background (z = 0 Slice)",
            page_subtitle="Orange dashed line: Sandia ECN G2-cold Experimental Schlieren Vapour Penetration | Injector at (0, 99.5) mm",
            ref_type="vapor",
        )

        print("Rendering Page 4: Side-by-side stages + 2D ECN experimental boundary contours...")
        plot_comparative_stages_page(
            pdf=pdf,
            snapshots=snapshots,
            mesh_info=mesh_info,
            ecn_data=ecn_data,
            bg_norm_T=bg_norm_T,
            bg_norm_Y=bg_norm_Y,
            p_norm_r=p_norm_r,
            p_norm_T=p_norm_T,
        )

        print("Rendering Page 5: Quantitative ECN Experimental Validation Benchmark...")
        plot_ecn_validation_page(pdf, snapshots, ecn_data, log_data)

        if log_data is not None and len(log_data["t_ms"]) > 0:
            print("Rendering Page 6: Global spray physics & penetration summary...")
            plot_summary_metrics_page(pdf, log_data, snapshots, mesh_info, ecn_data)

    print(f"Successfully generated 6-page PDF report: {os.path.abspath(out_pdf)}")


if __name__ == "__main__":
    main()
