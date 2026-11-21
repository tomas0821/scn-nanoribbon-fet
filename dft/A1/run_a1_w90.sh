#!/bin/bash
#SBATCH --job-name=scn-a1w
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_%x_%j.out
# A1 stage 2: HSE06 q 8³ IBZ + open_grid (production; criterion 2 failed at q 4³), then Wannier90 for
# pbe (PBE twin), hse_n (q4, route N), hse_o (q4, route O) and hse_q8 (q8, route O). Restartable.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /path/to/scn_fet/A1
mkdir -p outputs

for tag in hse_ibz_q8 hse_ibz_q8_open; do
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "[skip] $tag"; continue; fi
  echo "$(date '+%F %T') start $tag"
  if [ "$tag" = "hse_ibz_q8_open" ]; then
    srun --mpi=pmix -n 1 open_grid.x -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  else
    rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}"
    srun --mpi=pmix -n 64 pw.x -nk 8 -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  fi
  grep -q "JOB DONE" "$out" && echo "[done] $tag" || { echo "[FAIL] $tag"; exit 1; }
done

# the q8 open_grid k list goes into hse_q8.win, in open_grid's order
klist() { awk '/List to be put in the .win file of wannier90/{f=1; next} f && NF==4 {printf "  %s %s %s\n", $1, $2, $3; n++; next} f && n>0 {exit}' "$1"; }
if grep -q "@KLIST@" w90/hse_q8/hse_q8.win; then
  klist outputs/hse_ibz_q8_open.out > w90/hse_q8/klist.txt
  klist outputs/hse_ibz_open.out > w90/hse_q8/klist_q4.txt
  echo "q8 open_grid list: $(wc -l < w90/hse_q8/klist.txt) k-points; same order as q4: $(cmp -s w90/hse_q8/klist.txt w90/hse_q8/klist_q4.txt && echo yes || echo NO)"
  [ "$(wc -l < w90/hse_q8/klist.txt)" -eq 512 ] || { echo "[FAIL] q8 k list"; exit 1; }
  sed -e '/@KLIST@/{r w90/hse_q8/klist.txt' -e 'd}' w90/hse_q8/hse_q8.win > w90/hse_q8/hse_q8.win.new && mv w90/hse_q8/hse_q8.win.new w90/hse_q8/hse_q8.win
fi

for seed in pbe hse_n hse_o hse_q8; do
  cd "w90/${seed}"
  if [ -f "${seed}.wout" ] && grep -q "All done" "${seed}.wout"; then echo "[skip] w90_${seed}"; cd ../..; continue; fi
  echo "$(date '+%F %T') start w90_${seed}"
  srun --mpi=pmix -n 1 wannier90.x -pp "${seed}" < /dev/null > w90pp.log 2>&1
  srun --mpi=pmix -n 16 pw2wannier90.x -in "${seed}.pw2wan" < /dev/null > "${seed}.pw2wan.out" 2>&1
  srun --mpi=pmix -n 1 wannier90.x "${seed}" < /dev/null > w90.log 2>&1
  grep -q "All done" "${seed}.wout" 2>/dev/null && echo "[done] w90_${seed}" || echo "[FAIL] w90_${seed}"
  cd ../..
done
echo "$(date '+%F %T') stage 2 finished"
