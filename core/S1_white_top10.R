############################################################
# S1_white_top10.R
# ----------------------------------------------------------
# Supplementary code for extracting top-ranked white map
# candidates based on triad-derived priority scores.
#
# Input:
#   output/white_map_candidates_*.csv
#
# Output:
#   fig_suppl/S1_white_top10_bar.pdf
#   fig_suppl/S1_white_top10.csv
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

# -------------------------
# Paths (relative)
# -------------------------
DIR_OUT <- "output"
DIR_FIG <- "fig_suppl"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load white map candidates
# -------------------------
f <- list.files(DIR_OUT, pattern = "^white_map_candidates_.*\\.csv$", full.names = TRUE)
stopifnot(length(f) > 0)
f <- f[which.max(file.info(f)$mtime)]

WM <- read_csv(f, show_col_types = FALSE)

need <- c("A","C","AC_known","best_score_q","ABBC_strength","priority")
stopifnot(all(need %in% names(WM)))

# -------------------------
# Select A–C unknown candidates
# -------------------------
WM_u <- WM %>% filter(!AC_known)

if (nrow(WM_u) == 0) {
  warning("No A–C unknown pairs found; using all pairs.")
  WM_u <- WM
}

# -------------------------
# Top10 extraction
# -------------------------
TOP_N <- 10

top10 <- WM_u %>%
  arrange(desc(priority), desc(best_score_q), desc(ABBC_strength)) %>%
  slice_head(n = TOP_N) %>%
  mutate(pair = paste0(A, " → ", C))

# Save table
write_csv(
  top10 %>% select(A, C, priority, best_score_q, ABBC_strength),
  file.path(DIR_FIG, "S1_white_top10.csv")
)

# -------------------------
# Bar plot
# -------------------------
p_bar <- ggplot(
  top10,
  aes(x = priority, y = fct_reorder(pair, priority))
) +
  geom_col(width = 0.65, fill = "#4C72B0") +
  labs(
    title = "Top-ranked white map candidates",
    subtitle = "A–C unknown pairs ranked by triad-based priority",
    x = "Priority score",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(
  file.path(DIR_FIG, "S1_white_top10_bar.pdf"),
  p_bar,
  width = 7,
  height = max(3, 0.4 * nrow(top10) + 1)
)

message("=== DONE (S1_white_top10) ===")
