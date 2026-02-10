library(ggplot2)
library(dplyr)

size_replace_list <- list(
  "1048576" = "1M",
  "4194304" = "4M",
  "16777216" = "16M",
  "67108864" = "64M",
  "268435456" = "256M"
)
size_order <- c("1M", "4M", "16M", "64M", "256M")

df <- readr::read_tsv("time-size.tsv") %>%
  dplyr::mutate(CPU_TIME = SYSTEM + USER) %>%
  dplyr::mutate(
    DATA = ifelse(
      stringr::str_count(TEST_CASE, "genome") != 0,
      "GENOME",
      "TRANSCRIPTOME"
    ),
    SOFTWARE = stringr::str_extract(TEST_CASE, "^[^-]+"),
    RLEN = ifelse(stringr::str_count(TEST_CASE, "300") != 0, "300", "100"),
    SIZE = as.numeric(stringr::str_extract(
      TEST_CASE,
      "size([0-9]+)",
      group = 1
    ))
  ) %>%
  dplyr::mutate(
    SIZE = as.character(SIZE * ifelse(DATA == "GENOME", 1, 1024))
  )
print(unique(df$SIZE))

source("plot_legends.R")

df_p <- df %>%
  dplyr::select(CPU_TIME, WALL_CLOCK, RSS, DATA, SOFTWARE, RLEN, SIZE) %>%
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
    }),
    SIZE = sapply(SIZE, function(x) {
      size_replace_list[[x]]
    })
  ) %>%
  dplyr::group_by(ASPECTS, DATA, SOFTWARE, SIZE, RLEN) %>%
  dplyr::summarise(VALUES_MEAN = mean(VALUES), VALUES_SD = sd(VALUES)) %>%
  dplyr::ungroup()

for (rlen in c("100", "300")) {
  p <- df_p %>%
    dplyr::filter(RLEN == rlen) %>%
    ggplot() +
    geom_line(aes(
      y = VALUES_MEAN,
      x = factor(SIZE, levels = size_order),
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
        x = factor(SIZE, levels = size_order),
        color = factor(
          SOFTWARE,
          levels = software_levels
        ),
        group = SOFTWARE
      ),
      width = .2
    ) +
    scale_x_discrete(
      "# bases in the genome/transcriptome",
    ) +
    scale_y_continuous(
      "",
      trans = "log10",
      labels = scales::label_number(),
      limits = c(0.03, NA),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    scale_fill_discrete("Read length") +
    scale_color_manual("Software", values = software_colors) +
    facet_grid(ASPECTS ~ DATA, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", size = 16)
    )

  ggsave(
    paste0("fig/time_memory-size_", rlen, ".pdf"),
    p,
    width = 10,
    height = 10
  )
}

sessionInfo()
