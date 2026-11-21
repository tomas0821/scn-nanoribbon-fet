#!/bin/bash
#SBATCH --job-name=scn-a5hse
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A5-3 (make_a5.py hse): the HSE06 gap single point of one hold-out (dft/A5/README.md "A5-3"). These are the first
# three steps of dft/A4b/run_a4b2.sh, unchanged, with two differences:
#   - the rank count follows the allocation: 63 on one node, 126 on two (sbatch -N 2 -n 128 …); one pool per IBZ k
#     either way (9 k [110] → 9 × 7 or 9 × 14 ranks; 7 k [100] → 7 × 9 or 7 × 18);
#   - hse_<r> counts as finished only with QE's converged "!!" energy (a max_seconds stop also prints JOB DONE), and
#     if it is not finished the job stops there, so hpot/hproj never run on an unconverged state.
# Usage: sbatch -J a5hse_<r> [-N 2 -n 128] run_a5_hse.sh <r>        (tags from tags_hse_<r>.txt, in order)
#   hse_<r>     HSE06 SCF on the IBZ, fixed occupations, q 8 / 6, nbnd nocc + 32
#   hpot_<r>    pp.x bare + Hartree potential of the HSE density (cube, gzipped) → vacuum level
#   hproj_<r>   projwfc.x on the HSE SCF → band-edge character and the first vacuum-like state
# Restartable: a finished tag is skipped. There is no EXX restart: an HSE SCF stopped by max_seconds starts again.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd "${A5_HSE_DIR:-/path/to/scn_fet/A5/hse}"
R=${1:?usage: run_a5_hse.sh <ribbon>}
TAGS=$(cat "tags_hse_${R}.txt")
mkdir -p outputs
case "$R" in r110_*) NK=9 ;; *) NK=7 ;; esac
NR=$(( (SLURM_NTASKS / 64) * 63 ))
[ "$NR" -ge 63 ] || { echo "[FAIL] need a full node (SLURM_NTASKS=${SLURM_NTASKS})"; exit 1; }
echo "$(date '+%F %T') ${R}: ${NR} ranks, ${NK} pools"

finished() {   # finished <tag>
  local t=$1 o="outputs/$1.out"
  case "$t" in
    hse_*)  [ -f "$o" ] && grep -q "JOB DONE" "$o" && grep -q "^!!" "$o" ;;
    hpot_*) [ -s "hpot_${R}.cube.gz" ] ;;
    *)      [ -f "$o" ] && grep -q "JOB DONE" "$o" ;;
  esac
}

for tag in $TAGS; do
  if finished "$tag"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  inp="inputs/${tag}.in"
  out="outputs/${tag}.out"
  # </dev/null: srun would otherwise read stdin and swallow input (A0 job 141725)
  case "$tag" in
    hse_*)
      rm -rf "./tmp_${tag}" "hpot_${R}.dat" "hpot_${R}.cube" "hpot_${R}.cube.gz" "hpdos_${R}"* \
             outputs/{hpot,hproj}_"${R}".out
      mkdir -p "./tmp_${tag}"
      srun --mpi=pmix -n "$NR" pw.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1 ;;
    hpot_*)
      srun --mpi=pmix -n 16 pp.x -in "$inp" < /dev/null > "$out" 2>&1
      [ -s "hpot_${R}.cube" ] && gzip -f "hpot_${R}.cube" && rm -f "hpot_${R}.dat" ;;
    hproj_*)
      srun --mpi=pmix -n "$NR" projwfc.x -nk "$NK" -in "$inp" < /dev/null > "$out" 2>&1
      grep -q "JOB DONE" "$out" && rm -f "hpdos_${R}"* ;;
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
