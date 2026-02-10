library(ggplot2)
library(dplyr)

df <- readr::read_tsv("time-speedup.tsv") %>%
  dplyr::mutate(CPU_TIME = SYSTEM + USER) %>%
  dplyr::mutate(
    DATA = ifelse(
      stringr::str_count(TEST_CASE, "genome") != 0,
      "GENOME",
      "TRANSCRIPTOME"
    ),
    SOFTWARE = stringr::str_extract(TEST_CASE, "^[^-]+"),
    RLEN = ifelse(stringr::str_count(TEST_CASE, "300") != 0, "300", "100"),
    NTHREADS = stringr::str_extract(TEST_CASE, "nthreads([0-9]+)", group = 1)
  )

source("plot_legends.R")
df_p <- df %>%
  dplyr::select(CPU_TIME, WALL_CLOCK, RSS, DATA, SOFTWARE, RLEN, NTHREADS) %>%
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
  dplyr::group_by(ASPECTS, DATA, SOFTWARE, NTHREADS, RLEN) %>%
  dplyr::summarise(VALUES_MEAN = mean(VALUES), VALUES_SD = sd(VALUES)) %>%
  dplyr::ungroup()

for (rlen in c("100", "300")) {
  p <- df_p %>%
    dplyr::filter(RLEN == rlen) %>%
    ggplot() +
    geom_line(aes(
      y = VALUES_MEAN,
      x = factor(NTHREADS, levels = c("1", "2", "4", "8", "16")),
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
        x = factor(NTHREADS, levels = c("1", "2", "4", "8", "16")),
        color = factor(
          SOFTWARE,
          levels = software_levels
        ),
        group = SOFTWARE
      ),
      width = .2
    ) +
    scale_x_discrete(
      "# parallel threads"
    ) +
    scale_y_continuous(
      "",
      trans = "log10",
      labels = scales::label_number(),
      limits = c(1, NA),
      breaks = scales::breaks_log(n = 7),
      expand = expansion(mult = c(0.05, 0.2))
    ) +
    scale_fill_discrete("Read length") +
    scale_color_manual("Software", values = software_colors) +
    facet_grid(ASPECTS ~ DATA, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text=element_text(size=14),
      axis.title=element_text(size=16),
      legend.text=element_text(size=14),
      legend.title=element_text(size=16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(), 
      strip.text=element_text(color="black", size=16)
    )
  ggsave(paste0("fig/time_memory-speedup_", rlen, ".pdf"), p, width = 10, height = 10)
}

sessionInfo()
