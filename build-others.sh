#!/usr/bin/env bash
# shellcheck disable=SC1091
# shellcheck disable=SC2317
# shellcheck disable=SC2046

set -uex

if [ -n "${CFLAGS:-}" ]; then
    echo "CFLAGS is set to ${CFLAGS}"
else
    CFLAGS=(-O3 -mtune=native -march=native -w)
fi
if [ -n "${CC:-}" ]; then
    echo "CC is set to ${CC}"
else
    CC=icx
fi
if [ -n "${CXX:-}" ]; then
    echo "CXX is set to ${CXX}"
else
    CXX=icpx
fi
if [ "${CXX}" = "icpx" ]; then
set +ue
. /opt/intel/oneapi/setvars.sh
set -uex
fi


# Build HTSLib
env -C src/htslib-1.21 make distclean -j20 || true
env -C src/htslib-1.21 \
    ./configure --prefix="$(pwd)"/opt \
    CC="${CC}" \
    CFLAGS="${CFLAGS[*]}"
env -C src/htslib-1.21 make -j20
env -C src/htslib-1.21 make -j20 install

# Build GSL
env -C src/gsl-2.8 make distclean -j20 || true
env -C src/gsl-2.8 \
    ./configure --prefix="$(pwd)"/opt \
    --enable-shared=yes \
    --enable-static=yes \
    CC="${CC}" \
    CFLAGS="${CFLAGS[*]}"
env -C src/gsl-2.8 make -j20
env -C src/gsl-2.8 make -j20 install

# Build Original ART using optimized GSL
"${CXX}" "${CFLAGS[@]}" \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/art_original src/art_original/*.cpp \
    -lgsl -lgslcblas 

# Build wgsim
"${CXX}" "${CFLAGS[@]}" \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -o bin/wgsim src/wgsim.c \
    -lhts -lz -lpthread -lm -lc

# Build DWGSIM
"${CC}" "${CFLAGS[@]}" \
    -Lopt/lib/ \
    -Iopt/include \
    -Wl,-rpath,"$(pwd)"/opt/lib \
    -D_FILE_OFFSET_BITS=64 \
    -D_LARGEFILE64_SOURCE \
    -D_USE_KNETFILE \
    -DPACKAGE_VERSION='"0.1.15"' \
    -o bin/dwgsim src/dwgsim/*.c \
    -lhts -lm -lc

# Build PIRS
"${CXX}" "${CFLAGS[@]}" -fopenmp -std=c++17 \
    -DSFMT_MEXP=19937 -DHAVE_CONFIG_H -DPKGDATADIR='"/usr/local/share/pirs"' \
    -Isrc/pirs/SFMT-src-1.4 \
    -o bin/pirs \
    src/pirs/*.cpp src/pirs/SFMT-src-1.4/SFMT.c \
    -lz -lpthread

# Build generate_large_contigs and generate_more_contigs
for f in generate_large_contigs generate_more_contigs; do
    "${CC}" "${CFLAGS[@]}" -std=c11 \
        $(pkgconf --cflags mkl-sdl) \
        -o bin/"${f}" \
        "${f}".c \
        $(pkgconf --libs mkl-sdl)
done
