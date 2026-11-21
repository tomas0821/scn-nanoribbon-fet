#!/bin/bash
#SBATCH --job-name=scn-a2w
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%x_%j.out
# A2-2: HSE vacuum potential (pp.x) + Wannier90 for pbe_k16 and hse_k16 (make_a2_w90.py). pw2wannier90 on 16 ranks.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A2
echo "$(date '+%F %T') start hse_pot"
srun --mpi=pmix -n 16 pp.x -in inputs/hse_pot.in < /dev/null > outputs/hse_pot.out 2>&1
[ -s hse_pot.cube ] && echo "[done] hse_pot" || echo "[FAIL] hse_pot"
klist() { awk '/List to be put in the .win file of wannier90/{f=1; next} f && NF==4 {printf "  %s %s %s\n", $1, $2, $3; n++; next} f && n>0 {exit}' "$1"; }
if grep -q "@KLIST@" w90/hse_k16/hse_k16.win; then
  klist outputs/hse_ibz_k16_open.out > w90/hse_k16/klist.txt
  echo "open_grid list: $(wc -l < w90/hse_k16/klist.txt) k-points"
  [ "$(wc -l < w90/hse_k16/klist.txt)" -eq 256 ] || { echo "[FAIL] k list"; exit 1; }
  sed -e '/@KLIST@/{r w90/hse_k16/klist.txt' -e 'd}' w90/hse_k16/hse_k16.win > w90/hse_k16/hse_k16.win.new && mv w90/hse_k16/hse_k16.win.new w90/hse_k16/hse_k16.win
fi
for seed in pbe_k16 hse_k16; do
  cd "w90/$seed"
  echo "$(date '+%F %T') start w90_$seed"
  srun --mpi=pmix -n 1 wannier90.x -pp "$seed" < /dev/null > w90pp.log 2>&1
  srun --mpi=pmix -n 16 pw2wannier90.x -in "$seed.pw2wan" < /dev/null > "$seed.pw2wan.out" 2>&1
  srun --mpi=pmix -n 1 wannier90.x "$seed" < /dev/null > w90.log 2>&1
  grep -q "All done" "$seed.wout" && echo "[done] w90_$seed" || echo "[FAIL] w90_$seed"
  cd ../..
done
echo "$(date '+%F %T') finished"
