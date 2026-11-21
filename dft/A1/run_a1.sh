#!/bin/bash
#SBATCH --job-name=scn-a1
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A1 stage 1: bulk ScN, NC 80 Ry, PBE + HSE06 electronic structure for Wannier90 (stage 2 runs after inspection).
# Restartable: a tag whose output says "JOB DONE" is skipped. The tmp_ directories are KEPT: stage 2
# (pw2wannier90) reads tmp_pbe_scf (after pbe_nscf), tmp_hse_ibz (prefix hse_ibz_open) and tmp_hse_nosym.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd /path/to/scn_fet/A1
mkdir -p outputs

nk_of() {
  case "$1" in
    pbe_nscf|hse_nosym) echo 16 ;;
    *)                  echo 8 ;;
  esac
}

while read -r tag; do
  [ -z "$tag" ] && continue
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  case "$tag" in
    pbe_scf|hse_ibz|hse_nosym|hse_q8|hse_fock320) rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}" ;;
  esac
  # </dev/null: srun would otherwise read the loop's stdin and swallow the rest of tags.txt (A0 job 141725)
  if [ "$tag" = "hse_ibz_open" ]; then
    srun --mpi=pmix -n 1 open_grid.x -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  else
    srun --mpi=pmix -n 64 pw.x -nk "$(nk_of "$tag")" -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  fi
  grep -q "JOB DONE" "$out" && echo "[done] $tag" || echo "[FAIL] $tag"
done < tags.txt
echo "$(date '+%F %T') all tags processed: $(grep -c '^\[done\]' "slurm_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out" 2>/dev/null) done this run"
