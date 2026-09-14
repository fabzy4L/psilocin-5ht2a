#!/bin/bash
# Build GROMACS with GPU (CUDA) support inside a Linux/WSL2 environment,
# without needing sudo/apt (toolchain comes from conda, not the system).
#
# Verified end-to-end on: Windows 11 + WSL2 Ubuntu, RTX 5060 (Blackwell,
# CUDA arch compute_120), CUDA toolkit 12.9, GROMACS 2024.3. Full story of
# what this works around: ../docs/SESSION_HANDOFF_2026-09-10.md
#
# Usage:
#   bash build_gromacs.sh
# Override any default via env var, e.g.:
#   GMX_CUDA_ARCHITECTURES=89 INSTALL_PREFIX=/opt/gromacs bash build_gromacs.sh
#
# Requires: micromamba (or mambaforge/miniforge with `micromamba` on PATH)
# already installed. https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
set -e

GMX_VERSION="${GMX_VERSION:-2024.3}"
CUDA_TOOLKIT_VERSION="${CUDA_TOOLKIT_VERSION:-12.9}"
GMX_CUDA_ARCHITECTURES="${GMX_CUDA_ARCHITECTURES:-native}"
INSTALL_PREFIX="${INSTALL_PREFIX:-$HOME/gromacs}"
BUILD_DIR="${BUILD_DIR:-$HOME/build}"
ENV_NAME="${ENV_NAME:-gmx-build}"
RUN_CHECK="${RUN_CHECK:-0}"   # set to 1 to run `make check` (slow, optional)

command -v micromamba >/dev/null 2>&1 || {
  echo "ERROR: micromamba not found on PATH. Install it first:"
  echo "  https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
  exit 1
}

echo "== [1/5] Creating '$ENV_NAME' toolchain env: cmake, gcc/g++ 12, CUDA $CUDA_TOOLKIT_VERSION =="
echo "   (via conda, not apt -- avoids needing an interactive sudo password)"
micromamba create -y -n "$ENV_NAME" -c nvidia -c conda-forge \
  cmake make "gxx=12" "gcc=12" "cuda-toolkit=$CUDA_TOOLKIT_VERSION"

GMXENV="$HOME/micromamba/envs/$ENV_NAME"
if [ ! -d "$GMXENV" ]; then
  # fall back to whatever root prefix this micromamba install actually uses
  GMXENV="$(micromamba env list | awk -v n="$ENV_NAME" '$0 ~ n {print $NF; exit}')"
fi
[ -d "$GMXENV" ] || { echo "ERROR: could not locate the $ENV_NAME env directory"; exit 1; }

# Conda's split CUDA packaging needs these three additions beyond a normal
# system CUDA install: cicc lives under nvvm/bin (nvcc doesn't search it by
# default in this layout), and the headers/libs live under
# targets/x86_64-linux/ rather than directly under lib/include.
export PATH="$GMXENV/bin:$GMXENV/nvvm/bin:$PATH"
export CPATH="$GMXENV/targets/x86_64-linux/include:$CPATH"
export LIBRARY_PATH="$GMXENV/targets/x86_64-linux/lib:$LIBRARY_PATH"

echo "== [2/5] Fetching GROMACS $GMX_VERSION source =="
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
if [ ! -f "gromacs-$GMX_VERSION.tar.gz" ]; then
  wget -q "https://ftp.gromacs.org/gromacs/gromacs-$GMX_VERSION.tar.gz"
fi
[ -d "gromacs-$GMX_VERSION" ] || tar xf "gromacs-$GMX_VERSION.tar.gz"
cd "gromacs-$GMX_VERSION"
mkdir -p build && cd build

echo "== [3/5] Configuring =="
echo "   Using the classic Makefiles generator, not Ninja: GMX_BUILD_OWN_FFTW=ON"
echo "   fails outright under Ninja (GROMACS's bundled FFTW sub-build needs make)."
echo "   CUDA arch: $GMX_CUDA_ARCHITECTURES (native = auto-detect this machine's GPU)"
cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=CUDA \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DCMAKE_CUDA_ARCHITECTURES="$GMX_CUDA_ARCHITECTURES"

echo "== [4/5] Building (the long step -- 15-45+ min depending on cores/GPU) =="
NPROC="$(nproc)"
if ! make -j"$NPROC"; then
  echo ""
  echo "First link attempt failed -- retrying with an explicit cuFFT link fix."
  echo "(Known issue with conda-packaged CUDA toolkits: libgromacs.so builds"
  echo " fine with unresolved cufft* symbols -- shared libs tolerate that on"
  echo " Linux -- but the final gmx executable's link needs them resolved,"
  echo " and -lcufft doesn't get pulled in automatically. If this retry also"
  echo " fails, the error is something else -- check the log above, don't"
  echo " assume it's this same issue.)"
  CUFFT_FLAGS="-L$GMXENV/lib -Wl,-rpath,$GMXENV/lib -lcufft"
  cmake -DCMAKE_EXE_LINKER_FLAGS="$CUFFT_FLAGS" -DCMAKE_SHARED_LINKER_FLAGS="$CUFFT_FLAGS" .
  make -j"$NPROC"
fi

if [ "$RUN_CHECK" = "1" ]; then
  echo "== [4.5/5] Running GROMACS's regression test suite (RUN_CHECK=1) =="
  make check
fi

echo "== [5/5] Installing to $INSTALL_PREFIX =="
make install

echo ""
echo "== Verifying =="
# shellcheck disable=SC1091
source "$INSTALL_PREFIX/bin/GMXRC"
gmx --version | grep -iE 'GROMACS version|GPU support|CUDA driver|CUDA runtime'

echo ""
echo "Done. Add this to your shell profile (or source it manually each session):"
echo "  source $INSTALL_PREFIX/bin/GMXRC"
