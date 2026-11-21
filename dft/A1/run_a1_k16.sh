#!/bin/bash
#SBATCH --job-name=scn-a1k16
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%x_%j.out
# A1 production HSE06 Wannier on 16³ (make_a1_k16.py). Restartable for the two QE steps.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A1
for tag in hse_ibz_k16 hse_ibz_k16_open; do
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  if [ "$tag" = "hse_ibz_k16_open" ]; then
    srun --mpi=pmix -n 1 open_grid.x -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  else
    rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}"
    srun --mpi=pmix -n 64 pw.x -nk 16 -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  fi
  grep -q "JOB DONE" "$out" && echo "[done] $tag" || { echo "[FAIL] $tag"; exit 1; }
done
klist() { awk '/List to be put in the .win file of wannier90/{f=1; next} f && NF==4 {printf "  %s %s %s\n", $1, $2, $3; n++; next} f && n>0 {exit}' "$1"; }
if grep -q "@KLIST@" w90/hse_k16/hse_k16.win; then
  klist outputs/hse_ibz_k16_open.out > w90/hse_k16/klist.txt
  echo "open_grid list: $(wc -l < w90/hse_k16/klist.txt) k-points"
  [ "$(wc -l < w90/hse_k16/klist.txt)" -eq 4096 ] || { echo "[FAIL] k list"; exit 1; }
  sed -e '/@KLIST@/{r w90/hse_k16/klist.txt' -e 'd}' w90/hse_k16/hse_k16.win > w90/hse_k16/hse_k16.win.new && mv w90/hse_k16/hse_k16.win.new w90/hse_k16/hse_k16.win
fi
cd w90/hse_k16
echo "$(date '+%F %T') start w90_hse_k16"
srun --mpi=pmix -n 1 wannier90.x -pp hse_k16 < /dev/null > w90pp.log 2>&1
srun --mpi=pmix -n 16 pw2wannier90.x -in hse_k16.pw2wan < /dev/null > hse_k16.pw2wan.out 2>&1
srun --mpi=pmix -n 1 wannier90.x hse_k16 < /dev/null > w90.log 2>&1
grep -q "All done" hse_k16.wout && echo "[done] w90_hse_k16" || echo "[FAIL] w90_hse_k16"
cd ../..
echo "$(date '+%F %T') finished"
