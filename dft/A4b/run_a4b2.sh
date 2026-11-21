#!/bin/bash
#SBATCH --job-name=scn-a4b2
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A4b-2 job 1 (make_a4b.py hse): HSE06 of one calibration ribbon, NC set, up to pw2wannier90 (recipe W2's set B).
# Usage: sbatch -J a4b2_<ribbon> run_a4b2.sh <ribbon>        (tags from tags_hse_<ribbon>.txt, in order)
#   hse_<r>     HSE06 SCF on the IBZ (7 k [100] / 9 k [110]), fixed occupations, q 6 / 8, nbnd NB_HSE
#   hpot_<r>    pp.x bare + Hartree potential of the HSE density (cube, gzipped) → vacuum level (H1)
#   hproj_<r>   projwfc.x on the HSE SCF → first non-model state (H1, H3)
#   hopen_<r>   open_grid.x → the full 12 / 16-k grid (prefix hse_<r>_open), as A2
#   hw90pp_<r>  the open_grid k list into w90/<r>/H/<r>_H.win, then wannier90.x -pp
#   hp2w_<r>    pw2wannier90 (amn, mmn, eig)
# Pools: 63 ranks, one pool per IBZ k (7 × 9 ranks / 9 × 7 ranks); 8 pools would leave [110]'s ninth k alone.
# Restartable: a finished tag is skipped. If hse_<r> has to run again, every later product of that ribbon is removed
# first. There is no EXX restart: an HSE SCF stopped by max_seconds starts again from scratch.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd "${A4B_DIR:-/path/to/scn_fet/A4b}"
R=${1:?usage: run_a4b2.sh <ribbon>}
TAGS=$(cat "tags_hse_${R}.txt")
mkdir -p outputs
case "$R" in r110_*) NK=9; NKF=16 ;; *) NK=7; NKF=12 ;; esac
D="w90/${R}/H"

finished() {   # finished <tag>
  local t=$1 o="outputs/$1.out"
  case "$t" in
    hpot_*)   [ -s "hpot_${R}.cube.gz" ] ;;
    hw90pp_*) [ -s "${D}/${R}_H.nnkp" ] && ! grep -q "@KLIST@" "${D}/${R}_H.win" ;;
    hp2w_*)   [ -s "${D}/${R}_H.amn" ] && [ -s "${D}/${R}_H.mmn" ] && [ -s "${D}/${R}_H.eig" ] && grep -q "JOB DONE" "$o" 2>/dev/null ;;
    *)        [ -f "$o" ] && grep -q "JOB DONE" "$o" ;;
  esac
}

klist() { awk '/List to be put in the .win file of wannier90/{f=1; next} f && NF==4 {printf "  %s %s %s\n", $1, $2, $3; n++; next} f && n>0 {exit}' "$1"; }

for tag in $TAGS; do
  if finished "$tag"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  inp="inputs/${tag}.in"
  out="outputs/${tag}.out"
  # </dev/null: srun would otherwise read stdin and swallow input (A0 job 141725)
  case "$tag" in
    hse_*)
      rm -rf "./tmp_${tag}" "hpot_${R}.dat" "hpot_${R}.cube" "hpot_${R}.cube.gz" "hpdos_${R}"* \
             outputs/{hpot,hproj,hopen,hw90pp,hp2w}_"${R}".out "${D}"/*.{nnkp,amn,mmn,eig} "${D}"/*.pw2wan.out
      mkdir -p "./tmp_${tag}"
      srun --mpi=pmix -n 63 pw.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1 ;;
    hpot_*)
      srun --mpi=pmix -n 16 pp.x -in "$inp" < /dev/null > "$out" 2>&1
      [ -s "hpot_${R}.cube" ] && gzip -f "hpot_${R}.cube" && rm -f "hpot_${R}.dat" ;;
    hproj_*)
      srun --mpi=pmix -n 63 projwfc.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1
      grep -q "JOB DONE" "$out" && rm -f "hpdos_${R}"* ;;
    hopen_*)   srun --mpi=pmix -n 1 open_grid.x -in "$inp" < /dev/null > "$out" 2>&1 ;;
    hw90pp_*)
      if grep -q "@KLIST@" "${D}/${R}_H.win"; then
        klist "outputs/hopen_${R}.out" > "${D}/klist.txt"
        n=$(wc -l < "${D}/klist.txt")
        echo "open_grid list: $n k-points"
        if [ "$n" -ne "$NKF" ]; then echo "[FAIL] $tag: k list has $n points, expected $NKF"; continue; fi
        sed -e '/@KLIST@/{r '"${D}/klist.txt" -e 'd}' "${D}/${R}_H.win" > "${D}/${R}_H.win.new" && mv "${D}/${R}_H.win.new" "${D}/${R}_H.win"
      fi
      ( cd "$D" && srun --mpi=pmix -n 1 wannier90.x -pp "${R}_H" < /dev/null > w90pp.log 2>&1 )
      cp "${D}/w90pp.log" "$out" ;;
    hp2w_*)
      ( cd "$D" && srun --mpi=pmix -n 16 pw2wannier90.x -in "${R}_H.pw2wan" < /dev/null > "${R}_H.pw2wan.out" 2>&1 )
      cp "${D}/${R}_H.pw2wan.out" "$out" ;;
    *) echo "[FAIL] $tag: unknown tag"; continue ;;
  esac
  if finished "$tag"; then echo "[done] $tag"; else echo "[FAIL] $tag"; fi
done
echo "$(date '+%F %T') all tags processed"
