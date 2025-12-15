############################################################
# S1_white_map_core.R
# ----------------------------------------------------------
# Supplementary code for the core definition of the White map.
# Identifies A–C pairs without direct co-occurrence and ranks
# them using triad-based priority scores.
#
# Input:
#   output/abc_rankings_*.csv
#   data_proc/hits_matrix_*.csv
#
# Output:
#   fig_suppl/S1_white_map_core.pdf
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# -------------------------
# Paths (relative)
# -------------------------
DIR_OUT  <- "output"
DIR_PROC <- "data_proc"
DIR_FIG  <- "fig_suppl"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load inputs
# -------------------------
# Latest ABC rankings
f_abc <- list.files(DIR_OUT, pattern="^abc_rankings_.*\\.csv$", full.names=TRUE)
stopifnot(length(f_abc) > 0)
f_abc <- f_abc[which.max(file.info(f_abc)$mtime)]
ABC <- read_csv(f_abc, show_col_types = FALSE)

# Latest hits matrix
f_hits <- list.files(DIR_PROC, pattern="^hits_matrix_.*\\.csv$", full.names=TRUE)
stopifnot(length(f_hits) > 0)
f_hits <- f_hits[which.max(file.info(f_hits)$mtime)]
H <- read_csv(f_hits, show_col_types = FALSE)

# -------------------------
# Helper
# -------------------------
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))
get_hit <- function(df, term){
  col <- paste0("hit__", term)
  if (!col %in% names(df)) return(integer(nrow(df)))
  binv(df[[col]])
}

# -------------------------
# A–C known / unknown
# -------------------------
A_terms <- unique(ABC$A)
C_terms <- unique(ABC$C)

AC_mat <- expand_grid(A = A_terms, C = C_terms) %>%
  rowwise() %>%
  mutate(
    AC_n11 = sum(get_hit(H, A) & get_hit(H, C)),
    AC_known = AC_n11 > 0
  ) %>%
  ungroup()

# -------------------------
# Triad-based summary
# -------------------------
ABC_sum <- ABC %>%
  group_by(A, C) %>%
  summarise(
    best_score_q = max(score_q, na.rm = TRUE),
    AB_max = max(AB_n11, na.rm = TRUE),
    BC_max = max(BC_n11, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ABBC_strength = pmax(AB_max, 0) + pmax(BC_max, 0),
    priority = best_score_q + 0.1 * ABBC_strength
  )

# -------------------------
# White map core table
# -------------------------
WM <- AC_mat %>%
  left_join(ABC_sum, by = c("A","C")) %>%
  mutate(
    best_score_q   = coalesce(best_score_q, 0),
    ABBC_strength  = coalesce(ABBC_strength, 0),
    priority       = coalesce(priority, 0)
  ) %>%
  arrange(desc(!AC_known), desc(priority))

# -------------------------
# Plot (tile map)
# -------------------------
WM_plot <- WM %>%
  mutate(status = ifelse(AC_known, "known A–C", "A–C unknown"))

A_ord <- WM_plot %>%
  group_by(A) %>%
  summarise(m = mean(priority, na.rm = TRUE)) %>%
  arrange(desc(m)) %>%
  pull(A)

C_ord <- WM_plot %>%
  count(C) %>%
  arrange(desc(n)) %>%
  pull(C)

p_white <- ggplot(
  WM_plot,
  aes(x = factor(C, levels = C_ord),
      y = factor(A, levels = A_ord),
      fill = best_score_q)
) +
  geom_tile(color = "grey85") +
  geom_point(aes(size = ABBC_strength, shape = status), alpha = 0.85) +
  scale_shape_manual(values = c("A–C unknown" = 16, "known A–C" = 1)) +
  scale_fill_gradient(low = "#deebf7", high = "#08519c") +
  labs(
    title = "White map (core definition)",
    subtitle = "A–C unknown pairs ranked by triad-based priority",
    x = "Outcome (C)",
    y = "Candidate factor (A)",
    fill = "best_score_q",
    size = "ABBC strength",
    shape = "A–C status"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(DIR_FIG, "S1_white_map_core.pdf"),
  p_white,
  width = 10,
  height = 0.25 * length(unique(WM_plot$A)) + 3
)

message("=== DONE (S1_white_map_core) ===")
