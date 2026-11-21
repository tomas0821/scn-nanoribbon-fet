#!/bin/bash
#SBATCH --job-name=scn-a4b2w
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_%x_%j.out
# A4b-2 job 2: the HSE06 Wannier90 runs, W2 (the recipe) and W1 (reported), for the given ribbons, side by side;
# plus W3 (set S) where w90/<r>/H/<r>_HS.win exists (amendment 2026-11-14, ribbons whose H3 fired): its .amn is made
# first (wannier90 -pp + pw2wannier90 on the HSE open_grid wavefunctions); .mmn/.eig are linked from <r>_H.
# The variant .win files come from `parse_a4b.py hwindows` (w90/<r>/H/<v>/<r>_H<v>.win).
# Usage: sbatch run_a4b2_w90.sh r100_3 r100_3h r110_4 r110_4h
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd "${A4B_DIR:-/path/to/scn_fet/A4b}"
[ $# -gt 0 ] || { echo "usage: run_a4b2_w90.sh <ribbon>..."; exit 1; }
for R in "$@"; do
  H="w90/${R}/H"
  if [ -s "${H}/${R}_HS.win" ] && ! { [ -s "${H}/${R}_HS.amn" ] && grep -q "JOB DONE" "${H}/${R}_HS.pw2wan.out" 2>/dev/null; }; then
    echo "$(date '+%F %T') start ${R} set S amn"
    ( cd "$H" && srun --mpi=pmix -n 1 wannier90.x -pp "${R}_HS" < /dev/null > w90pp_S.log 2>&1 \
        && srun --mpi=pmix -n 16 pw2wannier90.x -in "${R}_HS.pw2wan" < /dev/null > "${R}_HS.pw2wan.out" 2>&1 )
    if [ -s "${H}/${R}_HS.amn" ] && grep -q "JOB DONE" "${H}/${R}_HS.pw2wan.out"; then echo "[done] ${R} set S amn"; else echo "[FAIL] ${R} set S amn"; fi
  fi
done
for R in "$@"; do
  for v in W1 W2 W3; do
    d="w90/${R}/H/${v}"
    s="${R}_H${v}"
    if [ ! -s "${d}/${s}.win" ]; then
      [ "$v" = W3 ] || echo "[FAIL] ${R} ${v}: no .win (run parse_a4b.py hwindows)"
      continue
    fi
    if [ -s "${d}/${s}.wout" ] && grep -q "All done" "${d}/${s}.wout"; then echo "[skip] ${R} ${v}"; continue; fi
    if [ "$v" = W3 ]; then
      [ -s "w90/${R}/H/${R}_HS.amn" ] || { echo "[FAIL] ${R} W3: no set-S .amn"; continue; }
      ln -sf "../${R}_HS.amn" "${d}/${s}.amn"
      for ext in mmn eig; do ln -sf "../${R}_H.${ext}" "${d}/${s}.${ext}"; done
    else
      for ext in amn mmn eig; do ln -sf "../${R}_H.${ext}" "${d}/${s}.${ext}"; done
    fi
    echo "$(date '+%F %T') start ${R} ${v}"
    ( cd "$d" && srun --mpi=pmix --exclusive -n 1 wannier90.x "$s" < /dev/null > w90.log 2>&1 ) &
  done
done
wait
for R in "$@"; do
  for v in W1 W2 W3; do
    [ "$v" = W3 ] && [ ! -s "w90/${R}/H/W3/${R}_HW3.win" ] && continue
    f="w90/${R}/H/${v}/${R}_H${v}.wout"
    if [ -s "$f" ] && grep -q "All done" "$f"; then echo "[done] ${R} ${v}"; else echo "[FAIL] ${R} ${v}"; fi
  done
done
echo "$(date '+%F %T') finished"
