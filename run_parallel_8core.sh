#!/usr/bin/env bash
set -e

cd /Users/dell/Documents/OpenFOAM/ECN_SprayG2

echo "================================================================================" > progress.log
echo " ECN Spray G2 — 1 mm Fine Grid Convergence Run (160x160x100 = 2.56M cells, 8 cores)" >> progress.log
echo " Target End Time: 1.5000 ms (0.0015 s)" >> progress.log
echo "================================================================================" >> progress.log
echo "[1/3] Running blockMesh (generating 2.56M cell mesh)..." >> progress.log

blockMesh > log.blockMesh 2>&1

echo "[2/3] Running decomposePar (splitting mesh across 8 cores)..." >> progress.log
decomposePar -force > log.decomposePar 2>&1

echo "[3/3] Starting parallel sprayFoam on 8 cores..." >> progress.log
echo "--------------------------------------------------------------------------------" >> progress.log
printf "%-10s | %-8s | %-9s | %-9s | %-8s | %-11s | %-10s\n" \
    "Time [ms]" "Progress" "Parcels" "Pen95[mm]" "SMD[um]" "T_gas [K]" "WallClock" >> progress.log
echo "--------------------------------------------------------------------------------" >> progress.log

mpirun -np 8 sprayFoam -parallel 2>&1 | tee log.sprayFoam | awk '
/^Time = / {
    t_s = $3 + 0.0
    t_ms = t_s * 1000.0
    pct = (t_s / 0.0015) * 100.0
}
/Current number of parcels/ {
    parcels = $NF
}
/D10, D32, Dmax/ {
    gsub(/,/, "", $(NF-1))
    d32 = $(NF-1) + 0.0
}
/Liquid penetration 95% mass/ {
    pen_mm = ($NF + 0.0) * 1000.0
}
/T gas min\/max/ {
    gsub(/,/, "", $(NF-1))
    tmin = $(NF-1) + 0.0
    tmax = $NF + 0.0
}
/^ExecutionTime = / {
    clock_s = $(NF-1)
    printf "%-10.4f | %6.1f%%  | %-9s | %9.1f | %8.1f | %5.1f-%5.1f | %7s s\n", \
        t_ms, pct, parcels, pen_mm, d32, tmin, tmax, clock_s >> "progress.log"
    fflush("progress.log")
}
/^End/ {
    echo_end = 1
}
'

echo "--------------------------------------------------------------------------------" >> progress.log
echo "[Post] Reconstructing parallel time directories (reconstructPar)..." >> progress.log
reconstructPar > log.reconstructPar 2>&1
echo "[Done] 1 mm fine-grid simulation & reconstruction complete!" >> progress.log

