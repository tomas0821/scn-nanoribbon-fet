#!/bin/bash
#SBATCH --job-name=scn-a2
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
# A2 stage 1: M-ScN sheet, NC 80 Ry, PBE + HSE06 electronic structure (make_a2.py). Wannier90 (stage 2) runs after
# inspection of the eigenvalues and the vacuum level. Restartable: a tag whose output says "JOB DONE" is skipped.
# tmp_pbe_scf (after pbe_nscf_k16) and tmp_hse_ibz_k16 (prefix hse_ibz_k16_open) are KEPT for pw2wannier90.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd /path/to/scn_fet/A2
mkdir -p outputs

while read -r tag; do
  [ -z "$tag" ] && continue
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q -E "JOB DONE|PP: DONE|PP *:" "$out"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  case "$tag" in
    pbe_scf|pbe_c30|hse_k8_q4|hse_k8_q8|hse_ibz_k16) rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}" ;;
  esac
  # </dev/null: srun would otherwise read the loop's stdin and swallow the rest of tags.txt (A0 job 141725)
  case "$tag" in
    hse_ibz_k16_open) srun --mpi=pmix -n 1 open_grid.x -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1 ;;
    pbe_pot)          srun --mpi=pmix -n 16 pp.x -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1 ;;
    *)                srun --mpi=pmix -n 64 pw.x -nk 8 -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1 ;;
  esac
  if grep -q "JOB DONE" "$out" || { [ "$tag" = "pbe_pot" ] && [ -s pbe_scf_pot.cube ]; }; then echo "[done] $tag"; else echo "[FAIL] $tag"; fi
done < tags.txt
echo "$(date '+%F %T') all tags processed"
