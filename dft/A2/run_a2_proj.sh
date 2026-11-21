#!/bin/bash
#SBATCH --job-name=scn-a2proj
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=slurm_%x_%j.out
# A4b amendment 2026-11-14, item 6 (diagnostic): projwfc.x on the A2 sheet's HSE06 SCF (IBZ 16×16) and PBE NSCF
# (tmp_pbe_scf holds the 16×16 NSCF), to see whether the sheet has a vacuum-like state near its CBM.
# The per-k pdos files are deleted; the stdout is what is parsed (dft/A2/sheet_nonmodel.py).
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd "${A2_DIR:-/path/to/scn_fet/A2}"
mkdir -p outputs
for t in hproj_sheet pproj_sheet; do
  if grep -q "JOB DONE" "outputs/${t}.out" 2>/dev/null; then echo "[skip] $t"; continue; fi
  echo "$(date '+%F %T') start $t"
  srun --mpi=pmix -n 64 projwfc.x -nk 8 -in "inputs/${t}.in" < /dev/null > "outputs/${t}.out" 2>&1
  if grep -q "JOB DONE" "outputs/${t}.out"; then echo "[done] $t"; rm -f hpdos_sheet* ppdos_sheet*; else echo "[FAIL] $t"; fi
done
echo "$(date '+%F %T') finished"
