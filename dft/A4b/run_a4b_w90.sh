#!/bin/bash
#SBATCH --job-name=scn-a4b1w
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_%x_%j.out
# A4b-1 job 2: the Wannier90 variants (W1, W2 on set B; W3 on set S) for the given ribbons, all at once (serial
# wannier90.x runs side by side). The variant .win files come from `parse_a4b.py windows` (w90/<r>/<v>/<r>_<v>.win);
# this script links the .amn/.mmn/.eig of the matching set into each variant directory.
# Usage: sbatch run_a4b_w90.sh r100_3 r100_3h r110_4 r110_4h
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd "${A4B_DIR:-/path/to/scn_fet/A4b}"
[ $# -gt 0 ] || { echo "usage: run_a4b_w90.sh <ribbon>..."; exit 1; }
for R in "$@"; do
  for v in W1 W2 W3; do
    d="w90/${R}/${v}"
    [ -s "${d}/${R}_${v}.win" ] || { echo "[FAIL] ${R} ${v}: no .win (run parse_a4b.py windows)"; continue; }
    if [ -s "${d}/${R}_${v}.wout" ] && grep -q "All done" "${d}/${R}_${v}.wout"; then echo "[skip] ${R} ${v}"; continue; fi
    case "$v" in W3) set_=S ;; *) set_=B ;; esac
    for ext in amn mmn eig; do ln -sf "../${R}_${set_}.${ext}" "${d}/${R}_${v}.${ext}"; done
    echo "$(date '+%F %T') start ${R} ${v} (set ${set_})"
    ( cd "$d" && srun --mpi=pmix --exclusive -n 1 wannier90.x "${R}_${v}" < /dev/null > w90.log 2>&1 ) &
  done
done
wait
for R in "$@"; do
  for v in W1 W2 W3; do
    f="w90/${R}/${v}/${R}_${v}.wout"
    if [ -s "$f" ] && grep -q "All done" "$f"; then echo "[done] ${R} ${v}"; else echo "[FAIL] ${R} ${v}"; fi
  done
done
echo "$(date '+%F %T') finished"
