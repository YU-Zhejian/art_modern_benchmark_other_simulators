#!/usr/bin/env bash
# shellcheck disable=SC1091
# shellcheck disable=SC2317
# shellcheck disable=SC2046

set +ue
. /opt/intel/oneapi/setvars.sh
set -uex

CFLAGS=(-O3 -mtune=native -march=native -w)

# Build HTSLib
env -C src/htslib-1.21 \
    ./configure --prefix="$(pwd)"/opt \
    CC=icx \
    CFLAGS="${CFLAGS[*]}"
env -C src/htslib-1.21 make -j20
env -C src/htslib-1.21 make -j20 install

# Build GSL
env -C src/gsl-2.8 \
    ./configure --prefix="$(pwd)"/opt \
    --enable-shared=yes \
    --enable-static=yes \
    CC=icx \
    CFLAGS="${CFLAGS[*]}"
env -C src/gsl-2.8 make -j20
env -C src/gsl-2.8 make -j20 install

# Build Original ART using optimized GSL
icpx "${CFLAGS[@]}" \
    -lgsl -lgslcblas \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/art_original src/art_original/*.cpp

# Build wgsim
icpx "${CFLAGS[@]}" \
    -lhts -lz -lpthread -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/wgsim src/wgsim.c

# Build DWGSIM
icx "${CFLAGS[@]}" \
    -lhts -lm -lc \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -D_FILE_OFFSET_BITS=64 \
    -D_LARGEFILE64_SOURCE \
    -D_USE_KNETFILE \
    -DPACKAGE_VERSION='"0.1.15"' \
    -o bin/dwgsim src/dwgsim/*.c

# Build PIRS
icpx "${CFLAGS[@]}" -fopenmp -std=c++17 \
    -lz -lpthread \
    -DSFMT_MEXP=19937 -DHAVE_CONFIG_H -DPKGDATADIR='"/usr/local/share/pirs"' \
    -Isrc/pirs/SFMT-src-1.4 \
    -o bin/pirs \
    src/pirs/*.cpp src/pirs/SFMT-src-1.4/SFMT.c

# Build generate_large_contigs and generate_more_contigs
for f in generate_large_contigs generate_more_contigs; do
    icx "${CFLAGS[@]}" -std=c11 \
        $(pkgconf --cflags --libs mkl-sdl) \
        -o bin/"${f}" \
        "${f}".c
done
