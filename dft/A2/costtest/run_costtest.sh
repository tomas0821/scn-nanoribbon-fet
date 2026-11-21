#!/bin/bash
#SBATCH --job-name=scn-hsecost
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --output=slurm_%x_%j.out
# HSE cost tests (README.md): T2 = wannier90.x on the A2 overlaps with lower outer windows (serial, run side by side),
# then T1 = four short HSE SCFs (fixed vs smearing, nbnd 24 vs 40) on the full node, one after another.
module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu          # pinned: a bare `module load qe` gives qe/7.5-gpu
cd /path/to/scn_fet/A2/costtest
mkdir -p outputs
echo "$(date '+%F %T') T2 start"
for d in w90/*_nb*/; do
  v=$(basename "$d"); seed=${v%_nb*}
  ( cd "$d" || exit 1
    for ext in amn mmn eig; do ln -sf "../../../w90/${seed}/${seed}.${ext}" .; done
    grep -q "All done" "${seed}.wout" 2>/dev/null && { echo "[skip] $v"; exit 0; }
    srun --exact -n 1 -c 1 --mpi=pmix wannier90.x "$seed" < /dev/null > w90.log 2>&1
    grep -q "All done" "${seed}.wout" && echo "[done] $v" || echo "[FAIL] $v" ) &
done
wait
echo "$(date '+%F %T') T2 finished; T1 start"
while read -r tag; do
  [ -z "$tag" ] && continue
  out="outputs/${tag}.out"
  if [ -f "$out" ] && grep -q "JOB DONE" "$out"; then echo "[skip] $tag"; continue; fi
  rm -rf "./tmp_${tag}"; mkdir -p "./tmp_${tag}"
  echo "$(date '+%F %T') start $tag"
  srun --mpi=pmix -n 64 pw.x -nk 8 -in "inputs/${tag}.in" < /dev/null > "$out" 2>&1
  grep -q "JOB DONE" "$out" && echo "[done] $tag" || echo "[FAIL] $tag"
  rm -rf "./tmp_${tag}"
done < tags_t1.txt
echo "$(date '+%F %T') finished"
