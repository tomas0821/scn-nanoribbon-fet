#!/bin/bash
#SBATCH --job-name=scn-a1d
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=slurm_%x_%j.out
# A1 Wannier diagnostic (criterion 4 failed): pbe_gi, pbe_ws, pbe9, pbe9_ws (see make_a1_diag.py). Wannier90 only,
# except one pw2wannier90 for the 9-WF projections. Each case writes <seed>_geninterp.dat on the QE band path.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A1/w90

w90run() {   # $1 = seed; runs in its directory
  cd "$1"
  srun --mpi=pmix -n 1 wannier90.x "$1" < /dev/null > w90.log 2>&1
  if grep -q "All done" "$1.wout" 2>/dev/null && [ -s "$1_geninterp.dat" ]; then echo "[done] $1"; else echo "[FAIL] $1"; fi
  cd ..
}
copy_mats() {   # $1 = from seed, $2 = to seed
  for ext in amn mmn eig; do cp "$1/$1.$ext" "$2/$2.$ext"; done
}

echo "$(date '+%F %T') start"
copy_mats pbe pbe_gi && w90run pbe_gi
copy_mats pbe pbe_ws && w90run pbe_ws
cd pbe9
srun --mpi=pmix -n 1 wannier90.x -pp pbe9 < /dev/null > w90pp.log 2>&1
srun --mpi=pmix -n 16 pw2wannier90.x -in pbe9.pw2wan < /dev/null > pbe9.pw2wan.out 2>&1
cd ..
w90run pbe9
copy_mats pbe9 pbe9_ws && w90run pbe9_ws
echo "$(date '+%F %T') finished"
