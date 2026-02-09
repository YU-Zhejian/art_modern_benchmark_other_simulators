#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set -ue
. run-functions.sh

OUT_TSV="time-speedup.tsv"

WGS_FCOV=4
TRANSCRIPTOME_FCOV=4

printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >"${OUT_TSV}"
for i in {1..3}; do
    echo "Run ${i}"
    for RLEN in 100 300; do
        export_on_adjusted_rlen
        for NTHREADS in 1 2 4 8 16; do
            run art_modern-genome-pe"${RLEN}"-nthreads"${NTHREADS}" opt/art_modern_build/art_modern \
                --mode wgs \
                --lc pe \
                --i-file data/ce11.fa \
                --i-fcov "${WGS_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${NTHREADS}"
            run art_modern_gcc-genome-pe"${RLEN}"-nthreads"${NTHREADS}" opt/art_modern_gcc_build/art_modern \
                --mode wgs \
                --lc pe \
                --i-file data/ce11.fa \
                --i-fcov "${WGS_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/ce11_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${NTHREADS}"

            run art_modern-transcriptome-pe"${RLEN}"-nthreads"${NTHREADS}" opt/art_modern_build/art_modern \
                --mode trans --lc pe \
                --i-file data/hg38_long_mrna.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${NTHREADS}"
            run art_modern_gcc-transcriptome-pe"${RLEN}"-nthreads"${NTHREADS}" opt/art_modern_gcc_build/art_modern  \
                --mode trans --lc pe \
                --i-file data/hg38_long_mrna.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${NTHREADS}"
        done
    done
done
rm -fr "${OUT_DIR}"

Rscript plot_time_memory-speedup.R
