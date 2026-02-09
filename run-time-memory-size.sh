#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue

export OUT_DIR=/tmp/data_out
export ART_MODERN_THREADS=6
export LD_LIBRARY_PATH="$(pwd)"/opt/lib:"${LD_LIBRARY_PATH:-}"
export LD_RUN_PATH="$(pwd)"/opt/lib:"${LD_RUN_PATH:-}"
FRAG_DIST_MEAN=500
FRAG_DIST_STD_DEV=20

WGS_FCOV=4
TRANSCRIPTOME_FCOV=4
GENOME_SIZES=("1048576" "4194304" "16777216" "67108864" "268435456")
TRANSCRIPTOME_SIZES=("1024" "4096" "16384" "65536" "262144")
TRANSCRIPT_LEN=1024

mkdir -p "${OUT_DIR}"
OUT_TSV="time-size.tsv"

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

function export_on_adjusted_rlen() {
    if [ "${RLEN}" -eq 100 ]; then
        export PIRS_B="data/e_coli_HiSeq2K_pirs_bcm.count.matrix"
        export ART_MTX_PREFIX="data/e_coli_HiSeq2K_art"
    else
        export PIRS_B="data/soybean_HiSeq2500_pirs_bcm.count.matrix"
        export ART_MTX_PREFIX="data/soybean_HiSeq2500_art"
    fi
}

for GENOME_SIZE in "${GENOME_SIZES[@]}"; do
    if [ ! -f "data/gen_genome_${GENOME_SIZE}bp.fa" ]; then
        echo "Generating genome ${GENOME_SIZE}bp"
        python3 gen_genome.py "${GENOME_SIZE}" >data/gen_genome_"${GENOME_SIZE}"bp.fa
    fi
done
for TRANSCRIPTOME_SIZE in "${TRANSCRIPTOME_SIZES[@]}"; do
    if [ ! -f "data/gen_transcriptome_${TRANSCRIPTOME_SIZE}bp.fa" ]; then
        echo "Generating transcriptome ${TRANSCRIPTOME_SIZE}bp"
        python3 gen_transcriptome.py "${TRANSCRIPTOME_SIZE}" "${TRANSCRIPT_LEN}" >data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa
    fi
done

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >"${OUT_TSV}"
for i in {1..3}; do
    echo "Run ${i}"
    for RLEN in 100 300; do
        export_on_adjusted_rlen
        for GENOME_SIZE in "${GENOME_SIZES[@]}"; do
            run wgsim-genome-pe"${RLEN}"-size"${GENOME_SIZE}" bin/wgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -N "$(echo "${WGS_FCOV} * $(grep -v '^$' <data/gen_genome_"${GENOME_SIZE}"bp.fa | wc -c) / 2  / ${RLEN}" | bc -q)" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -r 0 \
                data/gen_genome_"${GENOME_SIZE}"bp.fa \
                "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_wgsim_1.fq "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_wgsim_2.fq

            run pirs-genome-pe"${RLEN}"-size"${GENOME_SIZE}" bin/pirs simulate \
                -A dist \
                -m "${FRAG_DIST_MEAN}" \
                -l "${RLEN}" \
                -x "${WGS_FCOV}" \
                -v "${FRAG_DIST_STD_DEV}" \
                -t "${ART_MODERN_THREADS}" \
                -B "${PIRS_B}" \
                --no-gc-bias \
                --no-indel-errors \
                -o "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_Illumina \
                data/gen_genome_"${GENOME_SIZE}"bp.fa

            run dwgsim-genome-pe"${RLEN}"-size"${GENOME_SIZE}" bin/dwgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -C "${WGS_FCOV}" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -o 2 \
                -r 0 \
                -y 0 \
                data/gen_genome_"${GENOME_SIZE}"bp.fa \
                "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_dwgsim

            run art_original-genome-pe"${RLEN}"-size"${GENOME_SIZE}" bin/art_original \
                --in data/gen_genome_"${GENOME_SIZE}"bp.fa \
                --out "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_art_ \
                --qprof1 "${ART_MTX_PREFIX}"_R1.txt \
                --qprof2 "${ART_MTX_PREFIX}"_R2.txt \
                -f "${WGS_FCOV}" \
                --len "${RLEN}" \
                --mflen "${FRAG_DIST_MEAN}" \
                --sdev "${FRAG_DIST_STD_DEV}" \
                --noALN \
                --paired

            run art_modern-genome-pe"${RLEN}"-size"${GENOME_SIZE}" opt/art_modern_build/art_modern \
                --mode wgs \
                --lc pe \
                --i-file data/gen_genome_"${GENOME_SIZE}"bp.fa \
                --i-fcov "${WGS_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
            run art_modern_gcc-genome-pe"${RLEN}"-size"${GENOME_SIZE}" opt/art_modern_gcc_build/art_modern \
                --mode wgs \
                --lc pe \
                --i-file data/gen_genome_"${GENOME_SIZE}"bp.fa \
                --i-fcov "${WGS_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/gen_genome_"${GENOME_SIZE}"bp_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
        done
        for TRANSCRIPTOME_SIZE in "${TRANSCRIPTOME_SIZES[@]}"; do
            run wgsim-transcriptome-pe"${RLEN}"-size"${TRANSCRIPTOME_SIZE}" bin/wgsim \
                -1 "${RLEN}" -2 "${RLEN}" \
                -N "$(echo "${TRANSCRIPTOME_FCOV} * $(grep -v '^$' <data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa | wc -c) / 2  / ${RLEN}" | bc -q)" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -r 0 \
                data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa \
                "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_wgsim_1.fq \
                "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_wgsim_2.fq

            #    run pirs-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" [...]

            run dwgsim-transcriptome-pe"${RLEN}"-size"${TRANSCRIPTOME_SIZE}" bin/dwgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -C "${TRANSCRIPTOME_FCOV}" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -o 2 \
                -r 0 \
                -y 0 \
                data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa \
                "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_dwgsim

            run art_original-transcriptome-pe"${RLEN}"-size"${TRANSCRIPTOME_SIZE}" bin/art_original \
                --in data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa \
                --out "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_art_ \
                --qprof1 "${ART_MTX_PREFIX}"_R1.txt \
                --qprof2 "${ART_MTX_PREFIX}"_R2.txt \
                -f "${TRANSCRIPTOME_FCOV}" \
                --len "${RLEN}" \
                --mflen "${FRAG_DIST_MEAN}" \
                --sdev "${FRAG_DIST_STD_DEV}" \
                --noALN \
                --paired

            run art_modern-transcriptome-pe"${RLEN}"-size"${TRANSCRIPTOME_SIZE}" opt/art_modern_build/art_modern \
                --mode trans --lc pe \
                --i-file data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --i-parser stream \
                --o-fastq "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
            run art_modern_gcc-transcriptome-pe"${RLEN}"-size"${TRANSCRIPTOME_SIZE}" opt/art_modern_gcc_build/art_modern  \
                --mode trans --lc pe \
                --i-file data/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --i-parser stream \
                --o-fastq "${OUT_DIR}"/gen_transcriptome_"${TRANSCRIPTOME_SIZE}"bp_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
        done
    done
done
rm -fr "${OUT_DIR}"

Rscript plot_time_memory-size.R
