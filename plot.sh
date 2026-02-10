#!/usr/bin/env bash
set -ue
for fn in plot_time_memory*.R; do
  Rscript "${fn}"
done
# Convert all PDFs in fig/ to PNG using ghostscript
# for fn in fig/*.pdf; do
#     gs -dSAFER -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile="${fn%.pdf}.png" "${fn}"
# done

