#!/bin/bash
#SBATCH --job-name=scn-a4b1
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%x_%j.out
# A4b-1 job 1 (make_a4b.py): PBE twin of one calibration ribbon, NC set, up to pw2wannier90.
# Usage: sbatch -J a4b1_<ribbon> run_a4b.sh <ribbon>        (tags from tags_<ribbon>.txt, in order)
# Restartable: a finished tag is skipped. If scf_<ribbon> has to run again, every later product of that ribbon is
# removed first, so nothing downstream is reused from an older density. tmp_scf_<ribbon> is kept (bands, projwfc,
# nscf and pw2wannier90 all read it); the per-k pdos files are deleted after projwfc (the stdout is what is parsed).
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd "${A4B_DIR:-/path/to/scn_fet/A4b}"
R=${1:?usage: run_a4b.sh <ribbon>}
TAGS=$(cat "tags_${R}.txt")
mkdir -p outputs
case "$R" in r110_*) NK_NSCF=8 ;; *) NK_NSCF=4 ;; esac     # 16 k / 12 k on 64 ranks

finished() {   # finished <tag>
  local t=$1 o="outputs/$1.out"
  case "$t" in
    pot_*)   [ -s "pot_${R}.cube.gz" ] ;;
    w90pp_*) [ -s "w90/${R}/${R}_B.nnkp" ] && [ -s "w90/${R}/${R}_S.nnkp" ] ;;
    p2w_*_B) [ -s "w90/${R}/${R}_B.amn" ] && [ -s "w90/${R}/${R}_B.mmn" ] && [ -s "w90/${R}/${R}_B.eig" ] && grep -q "JOB DONE" "$o" 2>/dev/null ;;
    p2w_*_S) [ -s "w90/${R}/${R}_S.amn" ] && [ -e "w90/${R}/${R}_S.mmn" ] && [ -e "w90/${R}/${R}_S.eig" ] && grep -q "JOB DONE" "$o" 2>/dev/null ;;
    *)       [ -f "$o" ] && grep -q "JOB DONE" "$o" ;;
  esac
}

for tag in $TAGS; do
  if finished "$tag"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  inp="inputs/${tag}.in"
  out="outputs/${tag}.out"
  # </dev/null: srun would otherwise read stdin and swallow input (A0 job 141725)
  case "$tag" in
    scf_*)
      rm -rf "./tmp_${tag}" "pot_${R}.dat" "pot_${R}.cube" "pot_${R}.cube.gz" "pdos_${R}"* \
             outputs/{pot,bands,projwfc,nscf,w90pp,p2w}_"${R}"*.out "w90/${R}"/*.{nnkp,amn,mmn,eig} "w90/${R}"/*.pw2wan.out
      mkdir -p "./tmp_${tag}"
      srun --mpi=pmix -n 64 pw.x -nk 8 -in "$inp" < /dev/null > "$out" 2>&1 ;;
    pot_*)
      srun --mpi=pmix -n 16 pp.x -in "$inp" < /dev/null > "$out" 2>&1
      [ -s "pot_${R}.cube" ] && gzip -f "pot_${R}.cube" && rm -f "pot_${R}.dat" ;;
    bands_*)   srun --mpi=pmix -n 64 pw.x -nk 8 -in "$inp" < /dev/null > "$out" 2>&1 ;;
    projwfc_*)
      srun --mpi=pmix -n 64 projwfc.x -nk 8 -in "$inp" < /dev/null > "$out" 2>&1
      grep -q "JOB DONE" "$out" && rm -f "pdos_${R}"* ;;
    nscf_*)    srun --mpi=pmix -n 64 pw.x -nk "$NK_NSCF" -in "$inp" < /dev/null > "$out" 2>&1 ;;
    w90pp_*)
      ( cd "w90/${R}" && for s in B S; do srun --mpi=pmix -n 1 wannier90.x -pp "${R}_${s}" < /dev/null > "w90pp_${s}.log" 2>&1; done )
      cat "w90/${R}"/w90pp_*.log > "$out" ;;
    p2w_*_B|p2w_*_S)
      s=${tag##*_}
      ( cd "w90/${R}" && srun --mpi=pmix -n 16 pw2wannier90.x -in "${R}_${s}.pw2wan" < /dev/null > "${R}_${s}.pw2wan.out" 2>&1 )
      cp "w90/${R}/${R}_${s}.pw2wan.out" "$out"
      if [ "$s" = S ]; then ( cd "w90/${R}" && ln -sf "${R}_B.mmn" "${R}_S.mmn" && ln -sf "${R}_B.eig" "${R}_S.eig" ); fi ;;
    *) echo "[FAIL] $tag: unknown tag"; continue ;;
  esac
  if finished "$tag"; then echo "[done] $tag"; else echo "[FAIL] $tag"; fi
done
echo "$(date '+%F %T') all tags processed"
