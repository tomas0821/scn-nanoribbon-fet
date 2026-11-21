#!/bin/bash
#SBATCH --job-name=scn-a1d2b
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_%x_%j.out
# Round-2 rerun of the Wannier steps for pbe_k12 / pbe_k16 (the NSCFs in job 141911 finished). pw2wannier90 on 16
# ranks: on 64 the 36-plane FFT grid left ranks without planes ("processes with no planes") and it aborted.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A1
for n in 12 16; do
  cd "w90/pbe_k$n"
  srun --mpi=pmix -n 16 pw2wannier90.x -in "pbe_k$n.pw2wan" < /dev/null > "pbe_k$n.pw2wan.out" 2>&1
  srun --mpi=pmix -n 1 wannier90.x "pbe_k$n" < /dev/null > w90.log 2>&1
  grep -q "All done" "pbe_k$n.wout" && echo "[done] pbe_k$n" || echo "[FAIL] pbe_k$n"
  cd ../..
done
echo "$(date '+%F %T') finished"
