#!/bin/bash
#SBATCH --job-name=scn-a0
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A0 cutoff/pseudopotential tests: all inputs in tags.txt order, in one node allocation (one queue wait).
# Restartable: a tag whose output already says "JOB DONE" is skipped, so just resubmit after a timeout.
# Each tmp_ dir is removed after its last user (a *_bands run reuses its scf's save).
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd /path/to/scn_fet/A0
while read -r tag; do
  [ -z "$tag" ] && continue
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "[skip] $tag"; continue; fi
  case "$tag" in
    bulk_*)  nk=8 ;;
    *)       nk=4 ;;
  esac
  scf="${tag%_bands}"
  [ "$tag" = "$scf" ] && { rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}"; }
  echo "$(date '+%F %T') start $tag (nk=$nk)"
  # </dev/null: srun would otherwise read the loop's stdin and swallow the rest of tags.txt
  # (job 141725 ran only the first tag for this reason)
  srun --mpi=pmix -n 64 pw.x -nk $nk -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  grep -q "JOB DONE" "$out" && echo "[done] $tag" || echo "[FAIL] $tag"
  # keep the scf save only while its bands run still needs it
  if ! grep -qx "${tag}_bands" tags.txt; then rm -rf "./tmp_${scf}"; fi
done < tags.txt
echo "$(date '+%F %T') all tags processed"
