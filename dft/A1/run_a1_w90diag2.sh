#!/bin/bash
#SBATCH --job-name=scn-a1d2
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_%x_%j.out
# A1 Wannier diagnostic round 2 (make_a1_diag2.py): pbe_mlwf (8³, num_iter 200), pbe_k12, pbe_k16 (projection-only).
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A1
echo "$(date '+%F %T') start"

# (a) maximally localised, same matrices as pbe
for ext in amn mmn eig; do cp w90/pbe/pbe.$ext w90/pbe_mlwf/pbe_mlwf.$ext; done
cd w90/pbe_mlwf
srun --mpi=pmix -n 1 wannier90.x pbe_mlwf < /dev/null > w90.log 2>&1
grep -q "All done" pbe_mlwf.wout && echo "[done] pbe_mlwf" || echo "[FAIL] pbe_mlwf"
cd ../..

# (b) denser Wannier grids, each NSCF in its own copy of the SCF save directory
for n in 12 16; do
  rm -rf "tmp_pbe_k$n"; mkdir -p "tmp_pbe_k$n"; cp -r tmp_pbe_scf/pbe_scf.save "tmp_pbe_k$n/"
  srun --mpi=pmix -n 64 pw.x -nk 16 -in "inputs/pbe_nscf_k$n.in" < /dev/null > "outputs/pbe_nscf_k$n.out" 2>&1
  grep -q "JOB DONE" "outputs/pbe_nscf_k$n.out" || { echo "[FAIL] pbe_nscf_k$n"; continue; }
  cd "w90/pbe_k$n"
  srun --mpi=pmix -n 1 wannier90.x -pp "pbe_k$n" < /dev/null > w90pp.log 2>&1
  srun --mpi=pmix -n 16 pw2wannier90.x -in "pbe_k$n.pw2wan" < /dev/null > "pbe_k$n.pw2wan.out" 2>&1
  srun --mpi=pmix -n 1 wannier90.x "pbe_k$n" < /dev/null > w90.log 2>&1
  grep -q "All done" "pbe_k$n.wout" && echo "[done] pbe_k$n" || echo "[FAIL] pbe_k$n"
  cd ../..
done
echo "$(date '+%F %T') finished"
