replacement_list <- list(
  "CPU_TIME" = "CPU Time (s)",
  "WALL_CLOCK" = "Clock Time (s)",
  "RSS" = "Residential Memory (MB)",
  "MAJ_PG_F" = "Major Page Faults",
  "MIN_PG_F" = "Minor Page Faults",
  "VOL_CTX_S" = "Voluntary Context Switches",
  "IV_CTX_S" = "Involuntary Context Switches",
  "pirs" = "pIRS",
  "dwgsim" = "DWGSIM",
  "wgsim" = "wgsim",
  "art_original" = "Original ART",
  "art_modern" = "art_modern (Intel)",
  "art_modern_gcc" = "art_modern (GCC)"
)
software_levels <- c(
  "wgsim",
  "DWGSIM",
  "pIRS",
  "art_modern (Intel)",
  "art_modern (GCC)",
  "Original ART"
)
software_colors <- c(
  "wgsim" = "#1b9e77",
  "DWGSIM" = "#1fd902",
  "pIRS" = "#00caf7",
  "art_modern (Intel)" = "#e7298a",
  "art_modern (GCC)" = "#a61e47",
  "Original ART" = "#e6ab02"
)