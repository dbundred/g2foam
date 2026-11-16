#!/opt/homebrew/opt/python@3.13/bin/python3.13
"""
Side-by-side multi-case comparison PDF report generator for OpenFOAM ECN Spray G2 cases.

Compares N OpenFOAM cases (e.g., grid convergence, breakup/evaporation model variants,
or parametric sweeps) side-by-side against each other and against Sandia ECN experimental data.

Usage:
    # Compare two or more cases with optional custom labels (path:label):
    /opt/homebrew/opt/python@3.13/bin/python3.13 compare_cases_report.py \
        results_2mm_baseline:"2 mm Baseline (320k)" \
        .:"1 mm Fine Grid (2.56M)"

    # Specify custom output path or key slice timestamps (ms):
    /opt/homebrew/opt/python@3.13/bin/python3.13 compare_cases_report.py \
        caseA:"Case A" caseB:"Case B" caseC:"Case C" \
        --times 0.2 0.5 0.8 1.2 --out My_Comparison.pdf
"""

import argparse
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


# =============================================================================
# 1. Mesh & OpenFOAM Binary Readers
# =============================================================================

def parse_mesh_info(case_dir):
    """
    Determine domain bounds and (nx, ny, nz) from log.blockMesh (preferred, as it
    matches the actual mesh in case_dir) or system/blockMeshDict.
    """
    log_bm = os.path.join(case_dir, "log.blockMesh")
    if os.path.exists(log_bm):
        with open(log_bm, "r") as f:
            text = f.read()
        bb_m = re.search(
            r"boundingBox:\s*\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)\s*"
            r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)",
            text,
        )
        nc_m = re.search(r"nCells:\s*(\d+)", text)
        dx_m = re.search(r"i\s*:\s*([0-9.eE+-]+)", text)
        dz_m = re.search(r"j\s*:\s*([0-9.eE+-]+)", text)
        dy_m = re.search(r"k\s*:\s*([0-9.eE+-]+)", text)
        if bb_m and nc_m and dx_m and dy_m and dz_m:
            xmin, ymin, zmin, xmax, ymax, zmax = [float(bb_m.group(i)) for i in range(1, 7)]
            dx, dz, dy = float(dx_m.group(1)), float(dz_m.group(1)), float(dy_m.group(1))
            nx = int(round((xmax - xmin) / dx))
            ny = int(round((ymax - ymin) / dy))
            nz = int(round((zmax - zmin) / dz))
            return {
                "nx": nx, "ny": ny, "nz": nz,
                "ncells": int(nc_m.group(1)),
                "xmin": xmin, "xmax": xmax,
                "ymin": ymin, "ymax": ymax,
                "zmin": zmin, "zmax": zmax,
                "dx_mm": dx * 1000.0,
            }

    bmd_path = os.path.join(case_dir, "system", "blockMeshDict")
    if os.path.exists(bmd_path):
        with open(bmd_path, "r") as f:
            text = f.read()
        scale_m = re.search(r"scale\s+([0-9.eE+-]+)\s*;", text)
        scale = float(scale_m.group(1)) if scale_m else 1.0
        verts_m = re.search(r"vertices\s*\((.*?)\)\s*;", text, re.S)
        verts_raw = re.findall(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", verts_m.group(1))
        verts = np.array([[float(a), float(b), float(c)] for a, b, c in verts_raw]) * scale
        hex_m = re.search(r"hex\s*\([^)]+\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", text)
        nz, nx, ny = int(hex_m.group(1)), int(hex_m.group(2)), int(hex_m.group(3))
        xmin, xmax = verts[:, 0].min(), verts[:, 0].max()
        ymin, ymax = verts[:, 1].min(), verts[:, 1].max()
        zmin, zmax = verts[:, 2].min(), verts[:, 2].max()
        return {
            "nx": nx, "ny": ny, "nz": nz,
            "ncells": nx * ny * nz,
            "xmin": xmin, "xmax": xmax,
            "ymin": ymin, "ymax": ymax,
            "zmin": zmin, "zmax": zmax,
            "dx_mm": (xmax - xmin) * 1000.0 / nx,
        }

    raise FileNotFoundError(f"Could not find log.blockMesh or system/blockMeshDict in '{case_dir}'")


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
    """Read an OpenFOAM binary lagrangian scalarField."""
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
    dt = np.dtype([("lp", "S1"), ("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("cellI", "<i4"), ("rp", "S2")])
    recs = np.frombuffer(raw, dtype=dt, count=n, offset=m.end())
    return np.column_stack([recs["x"], recs["y"], recs["z"]])


def parse_sprayfoam_log(log_path):
    """Extract time-history metrics from log.sprayFoam."""
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


def load_case(case_dir, label, slice_half_width_m=0.0025, y_nozzle_mm=99.5):
    """Load all time steps and log metrics for a single OpenFOAM case."""
    mesh_info = parse_mesh_info(case_dir)
    nx, ny, nz = mesh_info["nx"], mesh_info["ny"], mesh_info["nz"]
    ncells = mesh_info["ncells"]

    # Find numeric time directories (> 0)
    time_dirs = []
    for name in os.listdir(case_dir):
        full_p = os.path.join(case_dir, name)
        if not os.path.isdir(full_p) or name in ("0", "0.orig"):
            continue
        try:
            t_val = float(name)
            if t_val > 0 and os.path.exists(os.path.join(full_p, "T")):
                time_dirs.append((t_val, name, full_p))
        except ValueError:
            continue
    time_dirs.sort(key=lambda x: x[0])

    if not time_dirs:
        raise RuntimeError(f"No valid time directories found in '{case_dir}'")

    dx_mm = (mesh_info["xmax"] - mesh_info["xmin"]) * 1000.0 / nx
    dy_mm = (mesh_info["ymax"] - mesh_info["ymin"]) * 1000.0 / ny
    x_centers = mesh_info["xmin"] * 1000.0 + (np.arange(nx) + 0.5) * dx_mm
    y_centers = mesh_info["ymin"] * 1000.0 + (np.arange(ny) + 0.5) * dy_mm
    z_axial_cells = y_nozzle_mm - y_centers

    snapshots = []
    for t_val, name, td_path in time_dirs:
        T_gas_3d = read_vol_scalar(os.path.join(td_path, "T"), ncells).reshape((ny, nx, nz))
        Y_fuel_3d = read_vol_scalar(os.path.join(td_path, "IC8H18"), ncells).reshape((ny, nx, nz))

        iz0, iz1 = nz // 2 - 1, nz // 2
        T_slice = 0.5 * (T_gas_3d[:, :, iz0] + T_gas_3d[:, :, iz1])
        Y_slice = 0.5 * (Y_fuel_3d[:, :, iz0] + Y_fuel_3d[:, :, iz1])

        cloud_dir = os.path.join(td_path, "lagrangian", "sprayCloud")
        pos = read_lagrangian_positions(os.path.join(cloud_dir, "positions"))
        d = read_lagrangian_scalar(os.path.join(cloud_dir, "d"))
        Tp = read_lagrangian_scalar(os.path.join(cloud_dir, "T"))
        npart = read_lagrangian_scalar(os.path.join(cloud_dir, "nParticle"))

        if len(pos) > 0:
            pos_mm_all = pos * 1000.0
            d_mm_all = d * 1000.0
            vol_p_mm3 = npart * (np.pi / 6.0) * (d_mm_all ** 3)
            z_axial_p = y_nozzle_mm - pos_mm_all[:, 1]

            order_z = np.argsort(z_axial_p)
            cum_z = np.cumsum(vol_p_mm3[order_z]) / np.sum(vol_p_mm3)
            z_liq_95 = float(z_axial_p[order_z[np.searchsorted(cum_z, 0.95)]])
            z_liq_99 = float(z_axial_p[order_z[np.searchsorted(cum_z, 0.99)]])

            xbins = np.linspace(mesh_info["xmin"] * 1000.0, mesh_info["xmax"] * 1000.0, nx + 1)
            zbins = np.linspace(0.0, (mesh_info["ymax"] - mesh_info["ymin"]) * 1000.0, ny + 1)
            H_vol, _, _ = np.histogram2d(z_axial_p, pos_mm_all[:, 0], bins=[zbins, xbins], weights=vol_p_mm3)
            PLV = H_vol / (dx_mm * dy_mm)
            z_bin_centers = 0.5 * (zbins[:-1] + zbins[1:])

            has_low = (PLV >= 0.2e-3).any(axis=1)
            plv_pen_low = float(z_bin_centers[has_low].max()) if has_low.any() else 0.0

            w_liq = {}
            for zt in [15.0, 25.0, 35.0]:
                iz_b = np.argmin(np.abs(z_bin_centers - zt))
                xs = x_centers[PLV[iz_b, :] >= 0.2e-3]
                w_liq[zt] = float(xs.max() - xs.min() + dx_mm) if len(xs) > 0 else 0.0

            mask = np.abs(pos[:, 2]) <= slice_half_width_m
            pos_slice = pos_mm_all[mask]
            r_um = 0.5 * d[mask] * 1e6
            Tp_slice = Tp[mask]
        else:
            pos_slice = np.zeros((0, 3))
            r_um = np.array([])
            Tp_slice = np.array([])
            z_liq_95, z_liq_99, plv_pen_low = 0.0, 0.0, 0.0
            w_liq = {15.0: 0.0, 25.0: 0.0, 35.0: 0.0}

        Y_proj = Y_fuel_3d.max(axis=2)
        has_vap = (Y_proj >= 0.001).any(axis=1)
        vap_pen = float(z_axial_cells[has_vap].max()) if has_vap.any() else 0.0

        w_vap = {}
        for zt in [15.0, 25.0, 35.0]:
            iy_c = np.argmin(np.abs(z_axial_cells - zt))
            xs = x_centers[Y_proj[iy_c, :] >= 0.001]
            w_vap[zt] = float(xs.max() - xs.min() + dx_mm) if len(xs) > 0 else 0.0

        snapshots.append({
            "time_s": t_val,
            "time_ms": t_val * 1000.0,
            "T_slice": T_slice,
            "Y_slice": Y_slice,
            "pos_mm": pos_slice,
            "r_um": r_um,
            "Tp": Tp_slice,
            "slice_parcels": len(pos_slice),
            "z_liq_95": z_liq_95,
            "z_liq_99": z_liq_99,
            "plv_pen_low": plv_pen_low,
            "vap_pen": vap_pen,
            "w_liq": w_liq,
            "w_vap": w_vap,
            "tp_mean": float(Tp_slice.mean()) if len(Tp_slice) else np.nan,
            "tp_min": float(Tp_slice.min()) if len(Tp_slice) else np.nan,
        })

    log_data = parse_sprayfoam_log(os.path.join(case_dir, "log.sprayFoam"))

    return {
        "dir": case_dir,
        "label": label,
        "mesh_info": mesh_info,
        "snapshots": snapshots,
        "log_data": log_data,
    }


def load_ecn_experimental_data(exp_dir="ECN_exp_data"):
    """Load Sandia ECN Spray G2-cold and G2 experimental datasets."""
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
        m_l = (l_low[:, 0] >= 0.0) & (l_low[:, 0] <= 1.6)
        m_v = (vap[:, 0] >= 0.0) & (vap[:, 0] <= 1.6)

        entry = {
            "t_liq": l_low[m_l, 0],
            "pen_liq_low": np.where(l_low[m_l, 15] > 0.5, l_low[m_l, 15], np.nan),
            "pen_liq_low_std": l_low[m_l, 16],
            "pen_liq_high": np.where(l_high[m_l, 15] > 0.5, l_high[m_l, 15], np.nan),
            "w15_liq": np.where(l_low[m_l, 1] > 0.5, l_low[m_l, 1], np.nan),
            "w15_liq_std": l_low[m_l, 2],
            "w25_liq": np.where(l_low[m_l, 5] > 0.5, l_low[m_l, 5], np.nan),
            "w25_liq_std": l_low[m_l, 6],
            "w35_liq": np.where(l_low[m_l, 9] > 0.5, l_low[m_l, 9], np.nan),
            "w35_liq_std": l_low[m_l, 10],
            "t_vap": vap[m_v, 0],
            "pen_vap": np.where(vap[m_v, 11] > 0.5, vap[m_v, 11], np.nan),
            "pen_vap_std": vap[m_v, 12],
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


# =============================================================================
# 2. Multi-Case Comparison Pages
# =============================================================================

CASE_COLORS = ["#dc2626", "#2563eb", "#16a34a", "#9333ea", "#d97706", "#0891b2", "#db2777"]
CASE_MARKERS = ["o", "s", "^", "D", "v", "P", "*"]


def plot_penetration_and_width_comparison(pdf, cases, ecn_data):
    """Page 1: N-case comparison of penetration & spray widths vs. Sandia ECN."""
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 10.0))
    fig.subplots_adjust(left=0.065, right=0.96, bottom=0.075, top=0.84, wspace=0.20, hspace=0.28)

    case_summary = "  vs.  ".join(
        f"{c['label']} ({c['mesh_info']['nx']}x{c['mesh_info']['ny']}x{c['mesh_info']['nz']})"
        for c in cases
    )
    fig.suptitle(
        f"1. Multi-Case Penetration & Width Comparison ({len(cases)} Cases vs. Sandia ECN)",
        fontsize=15,
        fontweight="bold",
        y=0.955,
    )
    fig.text(
        0.5,
        0.905,
        case_summary,
        ha="center",
        fontsize=10.2,
        color="#1e293b",
        bbox=dict(boxstyle="round,pad=0.38", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )

    g2c = ecn_data.get("G2C")

    # 1. Axial Liquid Penetration
    ax = axes[0, 0]
    if g2c is not None:
        ax.plot(g2c["t_liq"], g2c["pen_liq_low"], "-", color="#0f172a", linewidth=2.4,
                label="Sandia ECN G2-cold (PLV >= 0.2e-3)", zorder=3)
        ax.fill_between(g2c["t_liq"], g2c["pen_liq_low"] - g2c["pen_liq_low_std"],
                        g2c["pen_liq_low"] + g2c["pen_liq_low_std"], color="#94a3b8", alpha=0.25,
                        label="ECN ±1σ Band")

    for i, c in enumerate(cases):
        col = CASE_COLORS[i % len(CASE_COLORS)]
        mk = CASE_MARKERS[i % len(CASE_MARKERS)]
        ts = [s["time_ms"] for s in c["snapshots"]]
        ax.plot(ts, [s["z_liq_95"] for s in c["snapshots"]], f"{mk}-", color=col, linewidth=2.0,
                markersize=4.5, label=f"{c['label']} (Axial 95% mass)")
        ax.plot(ts, [s["z_liq_99"] for s in c["snapshots"]], "--", color=col, linewidth=1.3,
                alpha=0.75, label=f"{c['label']} (Axial 99% mass)")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Axial Liquid Penetration Length", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Axial Liquid Penetration Z [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 95)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.0, loc="upper left")

    # 2. Axial Vapour Penetration
    ax = axes[0, 1]
    if g2c is not None:
        ax.plot(g2c["t_vap"], g2c["pen_vap"], "-", color="#0f172a", linewidth=2.4,
                label="Sandia ECN G2-cold (Schlieren)", zorder=3)
        ax.fill_between(g2c["t_vap"], g2c["pen_vap"] - g2c["pen_vap_std"],
                        g2c["pen_vap"] + g2c["pen_vap_std"], color="#94a3b8", alpha=0.25,
                        label="ECN Schlieren ±1σ Band")

    for i, c in enumerate(cases):
        col = CASE_COLORS[i % len(CASE_COLORS)]
        mk = CASE_MARKERS[i % len(CASE_MARKERS)]
        ts = [s["time_ms"] for s in c["snapshots"]]
        ax.plot(ts, [s["vap_pen"] for s in c["snapshots"]], f"{mk}-", color=col, linewidth=2.0,
                markersize=4.8, label=f"{c['label']} (Y_IC8H18 >= 0.1%)")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Axial Vapour Penetration Length", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Axial Vapour Penetration Z [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 102)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="upper left")

    # 3. Liquid Spray Width at Z = 15 & 35 mm
    ax = axes[1, 0]
    if g2c is not None:
        ax.plot(g2c["t_liq"], g2c["w15_liq"], "-", color="#334155", linewidth=2.0, label="ECN Exp Liquid Width @ Z=15 mm")
        ax.plot(g2c["t_liq"], g2c["w35_liq"], "--", color="#334155", linewidth=2.0, label="ECN Exp Liquid Width @ Z=35 mm")

    for i, c in enumerate(cases):
        col = CASE_COLORS[i % len(CASE_COLORS)]
        mk = CASE_MARKERS[i % len(CASE_MARKERS)]
        ts = [s["time_ms"] for s in c["snapshots"]]
        w15 = np.array([s["w_liq"][15.0] for s in c["snapshots"]], dtype=float)
        w35 = np.array([s["w_liq"][35.0] for s in c["snapshots"]], dtype=float)
        w15[w15 == 0] = np.nan
        w35[w35 == 0] = np.nan
        ax.plot(ts, w15, f"{mk}-", color=col, linewidth=1.8, markersize=4.5, label=f"{c['label']} @ Z=15 mm")
        ax.plot(ts, w35, f"{mk}:", color=col, linewidth=1.8, markersize=4.5, label=f"{c['label']} @ Z=35 mm")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Liquid Spray Width at Z = 15 mm & Z = 35 mm", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Spray Plume Outer Width [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8, ncol=2, loc="upper left")

    # 4. Vapour Spray Width at Z = 15 & 35 mm
    ax = axes[1, 1]
    if g2c is not None:
        ax.plot(g2c["t_vap"], g2c["w15_vap"], "-", color="#334155", linewidth=2.0, label="ECN Exp Vapour Width @ Z=15 mm")
        ax.plot(g2c["t_vap"], g2c["w35_vap"], "--", color="#334155", linewidth=2.0, label="ECN Exp Vapour Width @ Z=35 mm")

    for i, c in enumerate(cases):
        col = CASE_COLORS[i % len(CASE_COLORS)]
        mk = CASE_MARKERS[i % len(CASE_MARKERS)]
        ts = [s["time_ms"] for s in c["snapshots"]]
        w15 = np.array([s["w_vap"][15.0] for s in c["snapshots"]], dtype=float)
        w35 = np.array([s["w_vap"][35.0] for s in c["snapshots"]], dtype=float)
        w15[w15 == 0] = np.nan
        w35[w35 == 0] = np.nan
        ax.plot(ts, w15, f"{mk}-", color=col, linewidth=1.8, markersize=4.5, label=f"{c['label']} @ Z=15 mm")
        ax.plot(ts, w35, f"{mk}:", color=col, linewidth=1.8, markersize=4.5, label=f"{c['label']} @ Z=35 mm")

    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Vapour Spray Width at Z = 15 mm & Z = 35 mm", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Vapour Plume Outer Width [mm]")
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 110)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.8, ncol=2, loc="upper left")

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def plot_global_physics_comparison(pdf, cases):
    """Page 2: N-case comparison of evaporation, SMD, and gas/droplet temperatures."""
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 10.0))
    fig.subplots_adjust(left=0.065, right=0.96, bottom=0.075, top=0.86, wspace=0.20, hspace=0.28)

    fig.suptitle(
        f"2. Multi-Case Thermodynamic, Evaporation & Breakup Comparison ({len(cases)} Cases)",
        fontsize=15,
        fontweight="bold",
        y=0.955,
    )

    # 1. Active Liquid & Evaporated Vapour Mass
    ax = axes[0, 0]
    for i, c in enumerate(cases):
        ld = c["log_data"]
        if ld is None or len(ld["t_ms"]) == 0:
            continue
        col = CASE_COLORS[i % len(CASE_COLORS)]
        ax.plot(ld["t_ms"], ld["liq_mass_mg"], "-", color=col, linewidth=2.0, label=f"{c['label']} — Liquid")
        ax.plot(ld["t_ms"], ld["evap_mass_mg"], "--", color=col, linewidth=2.0, label=f"{c['label']} — Vapour")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Fuel Mass Balance (Active Liquid vs. Evaporated Vapour)", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Mass [mg]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="best")

    # 2. Cloud Sauter Mean Radius (D32 / 2)
    ax = axes[0, 1]
    for i, c in enumerate(cases):
        ld = c["log_data"]
        if ld is None or len(ld["t_ms"]) == 0:
            continue
        col = CASE_COLORS[i % len(CASE_COLORS)]
        ax.plot(ld["t_ms"], ld["d32_um"] * 0.5, "-", color=col, linewidth=2.0, label=f"{c['label']} (D32 / 2)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2, label="EOI (0.78 ms)")
    ax.set_title("Cloud Sauter Mean Radius (D32 / 2)", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Sauter Mean Radius [um]")
    ax.set_ylim(15, 75)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5, loc="best")

    # 3. Gas Temperature Min / Max
    ax = axes[1, 0]
    for i, c in enumerate(cases):
        ld = c["log_data"]
        if ld is None or len(ld["t_ms"]) == 0:
            continue
        col = CASE_COLORS[i % len(CASE_COLORS)]
        ax.plot(ld["t_ms"], ld["t_gas_min"], "-", color=col, linewidth=2.0, label=f"{c['label']} — Min Gas T")
        ax.plot(ld["t_ms"], ld["t_gas_max"], "--", color=col, linewidth=1.6, alpha=0.8, label=f"{c['label']} — Max Gas T")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Eulerian Gas Temperature Extremes (Evaporative Cooling)", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Gas Temperature [K]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="best")

    # 4. Droplet Temperature Evolution (Slice Mean & Min)
    ax = axes[1, 1]
    for i, c in enumerate(cases):
        col = CASE_COLORS[i % len(CASE_COLORS)]
        mk = CASE_MARKERS[i % len(CASE_MARKERS)]
        ts = [s["time_ms"] for s in c["snapshots"]]
        ax.plot(ts, [s["tp_mean"] for s in c["snapshots"]], f"{mk}-", color=col, linewidth=1.8, markersize=4.5,
                label=f"{c['label']} — Mean T_p")
        ax.plot(ts, [s["tp_min"] for s in c["snapshots"]], f"{mk}:", color=col, linewidth=1.5, markersize=4.0,
                alpha=0.8, label=f"{c['label']} — Min T_p")
    ax.axhline(300.6, color="#0f766e", linestyle="-.", linewidth=1.1, label="IC8H18 Boiling Pt @ 50 kPa (300.6 K)")
    ax.axvline(0.78, color="gray", linestyle=":", linewidth=1.2)
    ax.set_title("Lagrangian Droplet Temperature (z = 0 Slice)", fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Time aSOI [ms]")
    ax.set_ylabel("Droplet Temperature [K]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.2, loc="best")

    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def plot_side_by_side_slice_pages(
    pdf,
    cases,
    ecn_data,
    target_times,
    bg_key,
    bg_cmap,
    bg_norm,
    bg_label,
    parcel_key,
    parcel_cmap,
    parcel_norm,
    parcel_label,
    title_prefix,
    subtitle,
    overlay_mode="liquid",
    y_nozzle_mm=99.5,
    max_cases_per_page=4,
):
    """
    Render side-by-side slice comparison pages:
      Rows = Cases (up to 4 per page, automatically paginated for larger N)
      Cols = Key timestamps (target_times)
    """
    g2c = ecn_data.get("G2C")
    n_cols = len(target_times)

    for start_idx in range(0, len(cases), max_cases_per_page):
        chunk = cases[start_idx : start_idx + max_cases_per_page]
        n_rows = len(chunk)

        fig_h = max(6.5, 2.75 * n_rows + 1.8)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16.5, fig_h), sharex=True, sharey=True, squeeze=False)
        fig.subplots_adjust(left=0.075, right=0.88, bottom=0.075, top=0.86, wspace=0.07, hspace=0.12)

        page_tag = f" (Cases {start_idx+1}-{start_idx+n_rows} of {len(cases)})" if len(cases) > max_cases_per_page else ""
        fig.suptitle(f"{title_prefix}{page_tag}", fontsize=14.5, fontweight="bold", y=0.96, x=0.47)
        fig.text(0.47, 0.915, subtitle, ha="center", fontsize=10.0, color="#222222")

        for r_idx, case in enumerate(chunk):
            mi = case["mesh_info"]
            extent = [mi["xmin"] * 1000.0, mi["xmax"] * 1000.0, mi["ymin"] * 1000.0, mi["ymax"] * 1000.0]

            for c_idx, tt in enumerate(target_times):
                ax = axes[r_idx, c_idx]
                snap = min(case["snapshots"], key=lambda s: abs(s["time_ms"] - tt))

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

                pos = snap["pos_mm"]
                if len(pos) > 0:
                    sizes = np.clip(2.0 + 0.34 * snap["r_um"], 2.2, 24.0)
                    ax.scatter(
                        pos[:, 0],
                        pos[:, 1],
                        c=snap[parcel_key],
                        cmap=parcel_cmap,
                        norm=parcel_norm,
                        s=sizes,
                        alpha=0.88,
                        edgecolors="black",
                        linewidths=0.18,
                        zorder=5,
                    )

                if g2c is not None:
                    t_ms = snap["time_ms"]
                    if overlay_mode == "liquid" and "liq_boundary_imgs" in g2c:
                        idx_f = int(np.argmin(np.abs(g2c["raw_t_liq"] - t_ms)))
                        b_img = g2c["liq_boundary_imgs"][idx_f]
                        Xg, Yg = np.meshgrid(
                            y_nozzle_mm - (np.arange(b_img.shape[1]) - 48.0) * 0.1908,
                            (np.arange(b_img.shape[0]) - 260.0) * 0.1908,
                        )
                        # Swap back to (X_cam, Y_cam)
                        X_cam = (np.arange(b_img.shape[0]) - 260.0) * 0.1908
                        Y_cam = y_nozzle_mm - (np.arange(b_img.shape[1]) - 48.0) * 0.1908
                        Yg2, Xg2 = np.meshgrid(Y_cam, X_cam)
                        ax.contour(Xg2, Yg2, b_img, levels=[0.5], colors=["#00ffff"], linewidths=1.3, linestyles="--", zorder=6)
                    elif overlay_mode == "vapor" and "schlieren_imgs" in g2c:
                        idx_v = int(np.argmin(np.abs(g2c["raw_t_vap"] - t_ms)))
                        s_img = g2c["schlieren_imgs"][idx_v]
                        X_cam = (np.arange(s_img.shape[0]) - 251.5) * 0.218
                        Y_cam = y_nozzle_mm - (np.arange(s_img.shape[1]) - 40.0) * 0.218
                        Yg2, Xg2 = np.meshgrid(Y_cam, X_cam)
                        ax.contour(Xg2, Yg2, s_img, levels=[50.0], colors=["#ff9900"], linewidths=1.4, linestyles="--", zorder=6)

                ax.set_xlim(-55, 55)
                ax.set_ylim(15, 102)
                if r_idx == 0:
                    inj_tag = "INJ ON" if snap["time_ms"] <= 0.78 else "POST-INJ"
                    ax.set_title(f"t = {snap['time_ms']:.1f} ms ({inj_tag})", fontsize=10.5, fontweight="bold", pad=5)
                if r_idx == n_rows - 1:
                    ax.set_xlabel("x [mm]", fontsize=9.5)
                if c_idx == 0:
                    ax.set_ylabel(f"{case['label']}\ny [mm]", fontsize=9.5, fontweight="bold")
                ax.tick_params(labelsize=8.2)

        cbar_ax1 = fig.add_axes([0.895, 0.50, 0.015, 0.35])
        cb1 = fig.colorbar(cm.ScalarMappable(norm=bg_norm, cmap=bg_cmap), cax=cbar_ax1)
        cb1.set_label(bg_label, fontsize=9.8, fontweight="semibold")
        cb1.ax.tick_params(labelsize=8.5)

        cbar_ax2 = fig.add_axes([0.895, 0.09, 0.015, 0.35])
        cb2 = fig.colorbar(cm.ScalarMappable(norm=parcel_norm, cmap=parcel_cmap), cax=cbar_ax2)
        cb2.set_label(parcel_label, fontsize=9.8, fontweight="semibold")
        cb2.ax.tick_params(labelsize=8.5)

        pdf.savefig(fig, dpi=200)
        plt.close(fig)


# =============================================================================
# 3. CLI & Entry Point
# =============================================================================

def parse_case_spec(spec):
    """Parse 'path:label' or 'path' into (path, label)."""
    if ":" in spec:
        path, label = spec.split(":", 1)
    else:
        path = spec
        norm_p = os.path.normpath(path)
        label = "Current Case" if norm_p == "." else os.path.basename(norm_p)
    return path, label


def main():
    parser = argparse.ArgumentParser(description="Generate a side-by-side multi-case OpenFOAM ECN Spray G2 PDF report.")
    parser.add_argument(
        "cases",
        nargs="*",
        help="Case directories as 'path' or 'path:Label'. Defaults to comparing results_2mm_baseline and current directory.",
    )
    parser.add_argument("--ecn", default="ECN_exp_data", help="Path to ECN experimental data directory.")
    parser.add_argument("--times", nargs="+", type=float, default=[0.2, 0.5, 0.8, 1.2], help="Key timestamps [ms] for slice comparisons.")
    parser.add_argument("--out", default="ECN_SprayG2_SideBySide_Report.pdf", help="Output PDF filename.")
    args = parser.parse_args()

    if not args.cases:
        default_specs = []
        if os.path.isdir("results_2mm_baseline"):
            default_specs.append('results_2mm_baseline:2 mm Baseline (320k)')
        default_specs.append('.:1 mm Fine Grid (2.56M)')
        case_specs = [parse_case_spec(s) for s in default_specs]
    else:
        case_specs = [parse_case_spec(s) for s in args.cases]

    cases = []
    for path, label in case_specs:
        print(f"Loading case '{label}' from {path}...")
        cases.append(load_case(path, label))

    ecn_data = load_ecn_experimental_data(args.ecn)

    bg_norm_T = Normalize(vmin=320.0, vmax=340.0)
    bg_norm_Y = Normalize(vmin=0.0, vmax=0.45)
    p_norm_r = Normalize(vmin=2.0, vmax=45.0)
    p_norm_T = Normalize(vmin=295.0, vmax=348.0)

    with PdfPages(args.out) as pdf:
        print("Rendering Page 1: Penetration & Width Comparison...")
        plot_penetration_and_width_comparison(pdf, cases, ecn_data)

        print("Rendering Page 2: Global Physics & Thermodynamics Comparison...")
        plot_global_physics_comparison(pdf, cases)

        print("Rendering Page 3: Side-by-Side Slice — Parcels (Radius) on Gas T...")
        plot_side_by_side_slice_pages(
            pdf=pdf,
            cases=cases,
            ecn_data=ecn_data,
            target_times=args.times,
            bg_key="T_slice",
            bg_cmap="inferno",
            bg_norm=bg_norm_T,
            bg_label="Gas Temperature T [K]",
            parcel_key="r_um",
            parcel_cmap="turbo",
            parcel_norm=p_norm_r,
            parcel_label="Parcel Droplet Radius r [um]",
            title_prefix="3. Side-by-Side Slice: Parcels (Radius) on Gas Temperature (z = 0)",
            subtitle="Rows = Cases | Cols = Time aSOI | Cyan dashed contour: Sandia ECN DBI Liquid Envelope (PLV = 0.2e-3)",
            overlay_mode="liquid",
        )

        print("Rendering Page 4: Side-by-Side Slice — Parcels (Droplet T) on Gas T...")
        plot_side_by_side_slice_pages(
            pdf=pdf,
            cases=cases,
            ecn_data=ecn_data,
            target_times=args.times,
            bg_key="T_slice",
            bg_cmap="inferno",
            bg_norm=bg_norm_T,
            bg_label="Gas Temperature T [K]",
            parcel_key="Tp",
            parcel_cmap="coolwarm",
            parcel_norm=p_norm_T,
            parcel_label="Parcel Droplet Temperature T_p [K]",
            title_prefix="4. Side-by-Side Slice: Parcels (Droplet T) on Gas Temperature (z = 0)",
            subtitle="Rows = Cases | Cols = Time aSOI | Cyan dashed contour: Sandia ECN DBI Liquid Envelope (PLV = 0.2e-3)",
            overlay_mode="liquid",
        )

        print("Rendering Page 5: Side-by-Side Slice — Parcels (Radius) on Fuel Mass Fraction...")
        plot_side_by_side_slice_pages(
            pdf=pdf,
            cases=cases,
            ecn_data=ecn_data,
            target_times=args.times,
            bg_key="Y_slice",
            bg_cmap="viridis",
            bg_norm=bg_norm_Y,
            bg_label="Iso-Octane Vapour Mass Fraction Y_IC8H18 [-]",
            parcel_key="r_um",
            parcel_cmap="plasma",
            parcel_norm=p_norm_r,
            parcel_label="Parcel Droplet Radius r [um]",
            title_prefix="5. Side-by-Side Slice: Parcels (Radius) on Fuel Vapour Mass Fraction (z = 0)",
            subtitle="Rows = Cases | Cols = Time aSOI | Orange dashed contour: Sandia ECN Schlieren 50% Vapour Envelope",
            overlay_mode="vapor",
        )

    print(f"Successfully generated side-by-side comparison report: {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()

