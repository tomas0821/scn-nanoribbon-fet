#!/bin/bash
#SBATCH --job-name=scn-a5
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A5: the hold-out ribbons (make_a5.py; README.md). Same logic as dft/A4a/run_a4a.sh.
# Usage: sbatch -J a5_<tag> run_a5.sh <tag>            one tag (A5-1: one relaxation per job, run in parallel)
#        sbatch run_a5.sh tags_<stage>.txt               a tag list (later stages)
# Restartable. A tag is done when its output has "End final coordinates" (relax_*) or "JOB DONE" (the rest): a relax
# stopped by max_seconds also prints JOB DONE, so that alone would wrongly mark it finished. An unfinished relax with
# its tmp dir is resumed with restart_mode='restart' (old output kept as .out.<n>). tmp_nm_* is kept for bands_* and
# projwfc_*; the other tmp dirs are removed after success (A4b runs the NC set, so PBE wavefunctions are not reused).
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd "${A5_DIR:-/path/to/scn_fet/A5}"
ARG=${1:?usage: run_a5.sh <tag | tag-list file>}
if [ -f "$ARG" ]; then TAGS=$(cat "$ARG"); else TAGS=$ARG; fi
mkdir -p outputs

finished() {   # finished <tag> <output>
  case "$1" in
    relax_*) grep -q "End final coordinates" "$2" ;;
    *)       grep -q "JOB DONE" "$2" ;;
  esac
}

for tag in $TAGS; do
  out="outputs/${tag}.out"
  inp="inputs/${tag}.in"
  if [ -f "$out" ] && finished "$tag" "$out"; then echo "[skip] $tag"; continue; fi
  case "$tag" in
    bands_*|projwfc_*) ;;                                   # reuse tmp_nm_<ribbon>
    relax_*)
      if [ -f "$out" ] && [ -d "./tmp_${tag}/${tag}.save" ]; then
        n=1; while [ -e "${out}.${n}" ]; do n=$((n + 1)); done
        mv "$out" "${out}.${n}"
        sed "s/calculation='relax',/calculation='relax', restart_mode='restart',/" "$inp" > "inputs/${tag}_restart.in"
        inp="inputs/${tag}_restart.in"
        echo "$(date '+%F %T') restart $tag (previous output -> ${out}.${n})"
      else
        rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}"
      fi ;;
    *) rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}" ;;
  esac
  echo "$(date '+%F %T') start $tag"
  # </dev/null: srun would otherwise read stdin and swallow input (A0 job 141725)
  case "$tag" in
    projwfc_*) srun --mpi=pmix -n 64 projwfc.x -nk 8 -in "$inp" < /dev/null > "$out" 2>&1 ;;
    *)         srun --mpi=pmix -n 64 pw.x -nk 8 -in "$inp" < /dev/null > "$out" 2>&1 ;;
  esac
  if finished "$tag" "$out"; then
    echo "[done] $tag"
    case "$tag" in nm_r*|bands_*|projwfc_*) ;; *) rm -rf "./tmp_${tag}" ;; esac
  else
    echo "[FAIL] $tag (not finished; resubmit to restart a relax)"
  fi
done
echo "$(date '+%F %T') all tags processed"
