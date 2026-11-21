#!/bin/bash
# Headless OVITO 3.17 (venv at $OVITO_VENV, default ~/.venvs/ovito) for the structure renders, within the laptop's 2.5 GB virtual-memory cap.
# libOpenGL from $OVITO_VENV/extra; Vulkan from Mesa lavapipe. MALLOC_ARENA_MAX and LP_NUM_THREADS keep the
# virtual address space small enough for `ulimit -v 2500000` (without them the render hangs or Qt fails to load).
# Usage (from the project root): bash scripts/ovito_env.sh scripts/render_ribbons.py render
OVITO_VENV=${OVITO_VENV:-$HOME/.venvs/ovito}
L=$(dirname "$(find "$OVITO_VENV/extra" -name 'libOpenGL.so.0' | head -1)")
ulimit -v 2500000
exec env -u DISPLAY -u WAYLAND_DISPLAY LD_LIBRARY_PATH="$L" QT_QPA_PLATFORM=ovitoheadless \
  VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json MALLOC_ARENA_MAX=2 LP_NUM_THREADS=2 OVITO_THREAD_COUNT=2 \
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 nice -n 10 timeout 900 \
  "$OVITO_VENV/bin/python" -P "$@"
