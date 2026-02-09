# shellcheck shell=bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR=/tmp/data_out
mkdir -p "${OUT_DIR}"
export ART_MODERN_THREADS=6
LD_LIBRARY_PATH="$(pwd)"/opt/lib:"${LD_LIBRARY_PATH:-}"
LD_RUN_PATH="$(pwd)"/opt/lib:"${LD_RUN_PATH:-}"
export LD_LIBRARY_PATH
export LD_RUN_PATH
export FRAG_DIST_MEAN=500
export FRAG_DIST_STD_DEV=20

function export_on_adjusted_rlen() {
    if [ "${RLEN}" -eq 100 ]; then
        export PIRS_B="data/e_coli_HiSeq2K_pirs_bcm.count.matrix"
        export ART_MTX_PREFIX="data/e_coli_HiSeq2K_art"
    else
        export PIRS_B="data/soybean_HiSeq2500_pirs_bcm.count.matrix"
        export ART_MTX_PREFIX="data/soybean_HiSeq2500_art"
    fi
}

function run() {
    echo "${1}" -- "${@:2}"
    /bin/time \
        -a \
        -o "${OUT_TSV}" \
        -f "${1}"'\t%e\t%S\t%U\t%M\t%F\t%R\t%w\t%c' \
        "${@:2}" &>"${OUT_DIR}"/"${1}".log || return 1
    rm -rf "${OUT_DIR:?}"/*
    mkdir -p "${OUT_DIR}"
}
