#!/bin/bash
# Array-C_G ToB runs (fits/PREREG_array.md): two processes, one per spacing, 2 threads and 2.5 GB cap each.
cd "$(dirname "$0")/../../.."   # repository root
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONPATH=.
for S in s2 sGAA; do
  ( ulimit -v 2500000; /usr/bin/time -f "$S: peak RSS %M kB, wall %E" nice -n 10 timeout 4h \
      ${PYTHON:-python} -P scripts/t1_cg_array.py tob $S > fits/t1/cg_array/tob_$S.log 2>&1;
    echo "exit $?" >> fits/t1/cg_array/tob_$S.log ) &
done
wait
