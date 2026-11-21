#!/bin/bash
#SBATCH --job-name=scn-a3hse
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%x_%j.out
# A3-2 HSE06 (make_a3.py stage2): the 3-plane slab, NC set, up to pw2wannier90 for the H (Sc d + N p) and HS (+ Sc s)
# projection sets. Steps as dft/A4b/run_a4b2.sh, with dft/A5/hse/run_a5_hse.sh's guard and two slab changes:
#   - hse_slab3 counts as finished only with QE's converged "!!" energy; if it is not finished the job stops there;
#   - k: 45 IBZ points of the 16 × 16 × 1 grid on 63 ranks (9 pools × 7); open_grid.x must give 256 points;
#   - both .win files (H, HS) get the open_grid k list; pw2wannier90 writes amn + mmn + eig for H and amn for HS
#     (mmn / eig linked from H, as A4b-1's set S).
# Usage: sbatch -J a3hse run_a3_hse.sh        (tags from tags_hse_slab3.txt, in order)
# Restartable: a finished tag is skipped. If hse_slab3 has to run again, every later product is removed first.
# There is no EXX restart: an HSE SCF stopped by max_seconds starts again from scratch.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd "${A3_DIR:-/path/to/scn_fet/A3}"
R=slab3
TAGS=$(cat "tags_hse_${R}.txt")
mkdir -p outputs
NR=63; NK=9; NKF=256
D="w90/${R}/H"

finished() {   # finished <tag>
  local t=$1 o="outputs/$1.out"
  case "$t" in
    hse_*)    [ -f "$o" ] && grep -q "JOB DONE" "$o" && grep -q "^!!" "$o" ;;
    hpot_*)   [ -s "hpot_${R}.cube.gz" ] ;;
    hw90pp_*) [ -s "${D}/${R}_H.nnkp" ] && [ -s "${D}/${R}_HS.nnkp" ] \
              && ! grep -q "@KLIST@" "${D}/${R}_H.win" && ! grep -q "@KLIST@" "${D}/${R}_HS.win" ;;
    hp2w_*)   [ -s "${D}/${R}_H.amn" ] && [ -s "${D}/${R}_H.mmn" ] && [ -s "${D}/${R}_H.eig" ] && [ -s "${D}/${R}_HS.amn" ] \
              && grep -q "JOB DONE" "${D}/${R}_H.pw2wan.out" 2>/dev/null && grep -q "JOB DONE" "${D}/${R}_HS.pw2wan.out" 2>/dev/null ;;
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
      srun --mpi=pmix -n "$NR" pw.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1 ;;
    hpot_*)
      srun --mpi=pmix -n 16 pp.x -in "$inp" < /dev/null > "$out" 2>&1
      [ -s "hpot_${R}.cube" ] && gzip -f "hpot_${R}.cube" && rm -f "hpot_${R}.dat" ;;
    hproj_*)
      srun --mpi=pmix -n "$NR" projwfc.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1
      grep -q "JOB DONE" "$out" && rm -f "hpdos_${R}"* ;;
    hopen_*)   srun --mpi=pmix -n 1 open_grid.x -in "$inp" < /dev/null > "$out" 2>&1 ;;
    hw90pp_*)
      klist "outputs/hopen_${R}.out" > "${D}/klist.txt"
      n=$(wc -l < "${D}/klist.txt")
      echo "open_grid list: $n k-points"
      if [ "$n" -ne "$NKF" ]; then echo "[FAIL] $tag: k list has $n points, expected $NKF"; break; fi
      for s in H HS; do
        if grep -q "@KLIST@" "${D}/${R}_${s}.win"; then
          sed -e '/@KLIST@/{r '"${D}/klist.txt" -e 'd}' "${D}/${R}_${s}.win" > "${D}/${R}_${s}.win.new" && mv "${D}/${R}_${s}.win.new" "${D}/${R}_${s}.win"
        fi
        ( cd "$D" && srun --mpi=pmix -n 1 wannier90.x -pp "${R}_${s}" < /dev/null > "w90pp_${s}.log" 2>&1 )
      done
      cat "${D}"/w90pp_H.log "${D}"/w90pp_HS.log > "$out" ;;
    hp2w_*)
      for s in H HS; do
        ( cd "$D" && srun --mpi=pmix -n 16 pw2wannier90.x -in "${R}_${s}.pw2wan" < /dev/null > "${R}_${s}.pw2wan.out" 2>&1 )
      done
      ( cd "$D" && ln -sf "${R}_H.mmn" "${R}_HS.mmn" && ln -sf "${R}_H.eig" "${R}_HS.eig" )
      cat "${D}/${R}_H.pw2wan.out" "${D}/${R}_HS.pw2wan.out" > "$out" ;;
    *) echo "[FAIL] $tag: unknown tag"; continue ;;
  esac
  if finished "$tag"; then
    echo "[done] $tag"
  else
    echo "[FAIL] $tag"
    case "$tag" in hse_*) echo "no converged '!!' energy: stopping (resubmit to rerun from scratch)"; break ;; esac
  fi
done
echo "$(date '+%F %T') all tags processed"
