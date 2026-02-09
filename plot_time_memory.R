library(ggplot2)
library(dplyr)

df <- readr::read_tsv("time.tsv") %>%
  dplyr::mutate(CPU_TIME = SYSTEM + USER) %>%
  dplyr::mutate(
    DATA = ifelse(
      stringr::str_count(TEST_CASE, "genome") != 0,
      "GENOME",
      "TRANSCRIPTOME"
    ),
    SOFTWARE = stringr::str_extract(TEST_CASE, "^[^-]+"),
    RLEN = ifelse(stringr::str_count(TEST_CASE, "300") != 0, "300", "100"),
    COVERAGE = stringr::str_extract(TEST_CASE, "cov([0-9]+)", group = 1)
  )

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
  "art_modern_prev_ver" = "art_modern (prev)",
  "art_modern_gcc" = "art_modern (GCC)",
  "art_modern_jemalloc" = "art_modern (master/jemalloc)",
  "art_modern_asio" = "art_modern (master/ASIO)",
  "art_modern_stlmap" = "art_modern (master/STLMAP)",
  "art_modern_mimalloc" = "art_modern (master/mi-malloc)"
)
software_levels <- c(
  "wgsim",
  "DWGSIM",
  "pIRS",
  "art_modern (Intel)",
  "art_modern (prev)",
  "art_modern (GCC)",
  "art_modern (master/jemalloc)",
  "art_modern (master/mi-malloc)",
  "art_modern (master/ASIO)",
  "art_modern (master/STLMAP)",
  "Original ART"
)
df_p <- df %>%
  dplyr::select(CPU_TIME, WALL_CLOCK, RSS, DATA, SOFTWARE, RLEN, COVERAGE) %>%
  dplyr::mutate(RSS = RSS / 1024) %>%
  tidyr::pivot_longer(
    cols = c("CPU_TIME", "WALL_CLOCK", "RSS"),
    names_to = "ASPECTS",
    values_to = "VALUES"
  ) %>%
  dplyr::mutate(
    ASPECTS = sapply(ASPECTS, function(x) {
      replacement_list[[x]]
    }),
    SOFTWARE = sapply(SOFTWARE, function(x) {
      replacement_list[[x]]
    })
  ) %>%
  dplyr::group_by(ASPECTS, DATA, SOFTWARE, COVERAGE, RLEN) %>%
  dplyr::summarise(VALUES_MEAN = mean(VALUES), VALUES_SD = sd(VALUES)) %>%
  dplyr::ungroup()

df_p %>%
  dplyr::filter(SOFTWARE == "art_modern (GCC)" | SOFTWARE == "Original ART") %>%
  dplyr::filter(ASPECTS == "Clock Time (s)") %>%
  dplyr::select(-VALUES_SD, -ASPECTS) %>%
  dplyr::group_by(DATA, SOFTWARE, COVERAGE, RLEN) %>%
  tidyr::pivot_wider(
    names_from = SOFTWARE,
    values_from = VALUES_MEAN
  ) %>%
  dplyr::mutate(SPEEDUP = `Original ART` / `art_modern (GCC)`) %>%
  dplyr::ungroup() %>%
  print()

for (rlen in c("100", "300")) {
  p <- df_p %>%
    dplyr::filter(RLEN == rlen) %>%
    ggplot() +
    geom_line(aes(
      y = VALUES_MEAN,
      x = factor(COVERAGE, levels = c("1", "2", "4", "8", "16")),
      color = factor(
        SOFTWARE,
        levels = software_levels
      ),
      group = SOFTWARE
    )) +
    geom_errorbar(
      aes(
        ymin = VALUES_MEAN - VALUES_SD,
        ymax = VALUES_MEAN + VALUES_SD,
        x = factor(COVERAGE, levels = c("1", "2", "4", "8", "16")),
        color = factor(
          SOFTWARE,
          levels = software_levels
        ),
        group = SOFTWARE
      ),
      width = .2
    ) +
    scale_x_discrete(
      "Coverage"
    ) +
    scale_y_continuous(
      "Values",
      trans = "log10",
      labels = scales::label_number(),
      # breaks = scales::breaks_pretty(n = 5),
      limits = c(0.1, NA)
    ) +
    scale_fill_discrete("Read length") +
    scale_color_discrete("Software") +
    facet_grid(ASPECTS ~ DATA, scales = "free_y") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(paste0("fig/time_memory_", rlen, ".pdf"), p, width = 10, height = 10)
}

sessionInfo()
