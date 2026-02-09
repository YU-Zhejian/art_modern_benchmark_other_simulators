#!/usr/bin/env bash
# shellcheck disable=SC1091
# shellcheck disable=SC2317
set -uex

rm -fr src/art_modern/
env -C src git clone https://github.com/YU-Zhejian/art_modern.git -b master --depth 1

rm -fr opt/art_modern_gcc_build/
mkdir -p opt/art_modern_gcc_build/
env -C opt/art_modern_gcc_build/ cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=PCG \
    -DUSE_MALLOC=NOP \
    -G Ninja "$(pwd)"/src/art_modern
env -C opt/art_modern_gcc_build/ ninja

set +uex
. /opt/intel/oneapi/setvars.sh
set -uex

rm -fr opt/art_modern_build/
mkdir -p opt/art_modern_build/
env -C opt/art_modern_build/ cmake \
    -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_COMPILER=icpx \
    -DCEU_CM_SHOULD_USE_NATIVE=ON \
    -DCEU_CM_SHOULD_ENABLE_TEST=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DUSE_RANDOM_GENERATOR=ONEMKL \
    -DUSE_HTSLIB=hts \
    -DUSE_MALLOC=NOP \
    -DCMAKE_PREFIX_PATH="$(pwd)"/opt \
    -DC_INCLUDE_PATH="$(pwd)"/opt/include \
    -G Ninja "$(pwd)"/src/art_modern
env -C opt/art_modern_build/ ninja
