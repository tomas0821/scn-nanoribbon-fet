#!/bin/bash
#SBATCH --job-name=scn-a5w
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%x_%j.out
# A5-2 job 2 + A3-3 (PBE twin part): the Wannier90 variants of the hold-outs, side by side. The steps are
# dft/A4b/run_a4b_w90.sh's (serial wannier90.x per variant; the variant .win from the windows step; the .amn/.mmn/.eig
# of the matching projection set linked in), with three differences:
#   - each argument is <dir>:<ribbon>, <dir> relative to scn_fet (the A5 twin and the A3 slab live in different
#     directories), so one node serves all of them;
#   - the variants are W1 W2 by default; W3 (set S) only where amendment 1's rule applies (VARIANTS="W1 W2 W3");
#   - 8 h limit: serial Wannier90 with 130–192 WFs is ≈ 5–9× A4b-1's (estimate 1.1–3.5 h for the longest run).
# Usage: sbatch run_a5_w90.sh A5/twin:r110_8 A5/twin:r100_6 A3:slab3
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
BASE=${SCN_FET:-/path/to/scn_fet}
VARIANTS=${VARIANTS:-W1 W2}
[ $# -gt 0 ] || { echo "usage: run_a5_w90.sh <dir>:<ribbon>..."; exit 1; }
for arg in "$@"; do
  D=${arg%%:*}; R=${arg##*:}
  for v in $VARIANTS; do
    d="${BASE}/${D}/w90/${R}/${v}"
    [ -s "${d}/${R}_${v}.win" ] || { echo "[FAIL] ${R} ${v}: no .win (run the windows step)"; continue; }
    if [ -s "${d}/${R}_${v}.wout" ] && grep -q "All done" "${d}/${R}_${v}.wout"; then echo "[skip] ${R} ${v}"; continue; fi
    case "$v" in W3) set_=S ;; *) set_=B ;; esac
    for ext in amn mmn eig; do ln -sf "../${R}_${set_}.${ext}" "${d}/${R}_${v}.${ext}"; done
    echo "$(date '+%F %T') start ${D} ${R} ${v} (set ${set_})"
    ( cd "$d" && srun --mpi=pmix --exclusive -n 1 wannier90.x "${R}_${v}" < /dev/null > w90.log 2>&1 ) &
  done
done
wait
for arg in "$@"; do
  D=${arg%%:*}; R=${arg##*:}
  for v in $VARIANTS; do
    f="${BASE}/${D}/w90/${R}/${v}/${R}_${v}.wout"
    if [ -s "$f" ] && grep -q "All done" "$f"; then echo "[done] ${R} ${v}"; else echo "[FAIL] ${R} ${v}"; fi
  done
done
echo "$(date '+%F %T') finished"
