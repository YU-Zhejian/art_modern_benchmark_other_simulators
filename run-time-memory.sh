#!/usr/bin/env bash
# shellcheck disable=SC2317
# shellcheck disable=SC1091
set -ue
. run-functions.sh
OUT_TSV="time.tsv"


printf 'TEST_CASE\tWALL_CLOCK\tSYSTEM\tUSER\tRSS\tMAJ_PG_F\tMIN_PG_F\tVOL_CTX_S\tIV_CTX_S\n' >"${OUT_TSV}"
for i in {1..3}; do
    echo "Run ${i}"
    for RLEN in 100 300; do
        export_on_adjusted_rlen
        for WGS_FCOV in 1 2 4 8 16; do
            run wgsim-genome-pe"${RLEN}"-cov"${WGS_FCOV}" bin/wgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -N "$(echo "${WGS_FCOV} * $(grep -v '^$' <data/ce11.fa | wc -c) / 2  / ${RLEN}" | bc -q)" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -r 0 \
                data/ce11.fa \
                "${OUT_DIR}"/ce11_wgsim_1.fq "${OUT_DIR}"/ce11_wgsim_2.fq

            run pirs-genome-pe"${RLEN}"-cov"${WGS_FCOV}" bin/pirs simulate \
                -A dist \
                -m "${FRAG_DIST_MEAN}" \
                -l "${RLEN}" \
                -x "${WGS_FCOV}" \
                -v "${FRAG_DIST_STD_DEV}" \
                -t "${ART_MODERN_THREADS}" \
                -B "${PIRS_B}" \
                --no-gc-bias \
                --no-indel-errors \
                -o "${OUT_DIR}"/Illumina \
                data/ce11.fa

            run dwgsim-genome-pe"${RLEN}"-cov"${WGS_FCOV}" bin/dwgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -C "${WGS_FCOV}" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -o 2 \
                -r 0 \
                -y 0 \
                data/ce11.fa \
                "${OUT_DIR}"/ce11_dwgsim

            run art_original-genome-pe"${RLEN}"-cov"${WGS_FCOV}" bin/art_original \
                --in data/ce11.fa \
                --out "${OUT_DIR}"/ce11_art_ \
                --qprof1 "${ART_MTX_PREFIX}"_R1.txt \
                --qprof2 "${ART_MTX_PREFIX}"_R2.txt \
                -f "${WGS_FCOV}" \
                --len "${RLEN}" \
                --mflen "${FRAG_DIST_MEAN}" \
                --sdev "${FRAG_DIST_STD_DEV}" \
                --noALN \
                --paired

            run art_modern-genome-pe"${RLEN}"-cov"${WGS_FCOV}" opt/art_modern_build/art_modern \
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
                --parallel "${ART_MODERN_THREADS}"
            run art_modern_gcc-genome-pe"${RLEN}"-cov"${WGS_FCOV}" opt/art_modern_gcc_build/art_modern \
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
                --parallel "${ART_MODERN_THREADS}"
        done
        for TRANSCRIPTOME_FCOV in 1 2 4 8 16; do

            run wgsim-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" bin/wgsim \
                -1 "${RLEN}" -2 "${RLEN}" \
                -N "$(echo "${TRANSCRIPTOME_FCOV} * $(grep -v '^$' <data/hg38_long_mrna.fa | wc -c) / 2  / ${RLEN}" | bc -q)" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -r 0 \
                data/hg38_long_mrna.fa \
                "${OUT_DIR}"/hg38_long_mrna_wgsim_1.fq \
                "${OUT_DIR}"/hg38_long_mrna_wgsim_2.fq

            #    run pirs-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" [...]

            run dwgsim-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" bin/dwgsim \
                -1 "${RLEN}" \
                -2 "${RLEN}" \
                -C "${TRANSCRIPTOME_FCOV}" \
                -d "${FRAG_DIST_MEAN}" \
                -s "${FRAG_DIST_STD_DEV}" \
                -o 2 \
                -r 0 \
                -y 0 \
                data/hg38_long_mrna.fa \
                "${OUT_DIR}"/hg38_long_mrna_dwgsim

            run art_original-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" bin/art_original \
                --in data/hg38_long_mrna.fa \
                --out "${OUT_DIR}"/hg38_long_mrna_art_ \
                --qprof1 "${ART_MTX_PREFIX}"_R1.txt \
                --qprof2 "${ART_MTX_PREFIX}"_R2.txt \
                -f "${TRANSCRIPTOME_FCOV}" \
                --len "${RLEN}" \
                --mflen "${FRAG_DIST_MEAN}" \
                --sdev "${FRAG_DIST_STD_DEV}" \
                --noALN \
                --paired

            run art_modern-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" opt/art_modern_build/art_modern \
                --mode trans --lc pe \
                --i-file data/hg38_long_mrna.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --i-parser stream \
                --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
            run art_modern_gcc-transcriptome-pe"${RLEN}"-cov"${TRANSCRIPTOME_FCOV}" opt/art_modern_gcc_build/art_modern  \
                --mode trans --lc pe \
                --i-file data/hg38_long_mrna.fa \
                --i-fcov "${TRANSCRIPTOME_FCOV}" \
                --read_len "${RLEN}" \
                --i-parser stream \
                --o-fastq "${OUT_DIR}"/hg38_long_mrna_art_modern_wgs_memory.fastq \
                --qual_file_1 "${ART_MTX_PREFIX}"_R1.txt \
                --qual_file_2 "${ART_MTX_PREFIX}"_R2.txt \
                --pe_frag_dist_mean "${FRAG_DIST_MEAN}" \
                --pe_frag_dist_std_dev "${FRAG_DIST_STD_DEV}" \
                --parallel "${ART_MODERN_THREADS}"
        done
    done

done
rm -fr "${OUT_DIR}"

Rscript plot_time_memory.R
