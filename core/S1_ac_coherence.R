############################################################
# S1_AC_coherence.R
# ----------------------------------------------------------
# Supplementary code for A–C coherence (NPMI / lift)
# Reproduces Fig3 using precomputed co-occurrence statistics.
#
# Input:
#   output/cooc_npmi_lift_*.csv
#
# Output:
#   fig_suppl/S1_AC_coherence_heatmap.pdf
#   fig_suppl/S1_AC_coherence_scatter.pdf
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# -------------------------
# Paths (relative)
# -------------------------
DIR_OUT   <- "output"
DIR_FIG   <- "fig_suppl"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load data
# -------------------------
# Use the latest cooc_npmi_lift file
f <- list.files(DIR_OUT, pattern = "^cooc_npmi_lift_.*\\.csv$", full.names = TRUE)
stopifnot(length(f) > 0)
f <- f[which.max(file.info(f)$mtime)]

D <- read_csv(f, show_col_types = FALSE)

stopifnot(all(c("A","C","npmi","lift","pAC") %in% names(D)))

# Derived metric
D <- D %>%
  mutate(log2_lift = log2(pmax(lift, 1e-12)))

# -------------------------
# Filters (as in main analysis)
# -------------------------
TH_NPMI <- 0.05
TH_L2L  <- log2(1.2)

D_f <- D %>%
  filter(is.finite(npmi), is.finite(log2_lift))

# -------------------------
# Heatmap (two panels)
# -------------------------
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

H1 <- D_f %>% mutate(metric = "NPMI", value = clip(npmi, -0.2, 1.0))
H2 <- D_f %>% mutate(metric = "log2(Lift)", value = clip(log2_lift, -2, 2))
H  <- bind_rows(H1, H2) %>%
  mutate(metric = factor(metric, levels = c("NPMI","log2(Lift)")))

# Order for readability
A_ord <- H %>% group_by(A) %>% summarise(m = mean(value, na.rm=TRUE)) %>% arrange(m) %>% pull(A)
C_ord <- H %>% count(C) %>% arrange(desc(n)) %>% pull(C)

H <- H %>%
  mutate(A = factor(A, levels = A_ord),
         C = factor(C, levels = C_ord))

p_heat <- ggplot(H, aes(x = A, y = C, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0) +
  labs(title = "A–C coherence", subtitle = basename(f), x = "A term", y = "Outcome (C)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(DIR_FIG, "S1_AC_coherence_heatmap.pdf"),
       p_heat, width = 10, height = 8)

# -------------------------
# Scatter
# -------------------------
p_sc <- ggplot(D_f, aes(x = npmi, y = log2_lift)) +
  geom_vline(xintercept = TH_NPMI, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = TH_L2L,  linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = pAC, color = C), alpha = 0.8, show.legend = TRUE) +
  scale_size_continuous(name = "p(A,C)", range = c(1.5, 5)) +
  labs(title = "A–C coherence scatter",
       x = "NPMI", y = "log2(Lift)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(DIR_FIG, "S1_AC_coherence_scatter.pdf"),
       p_sc, width = 10, height = 7)

message("=== DONE (S1_AC_coherence) ===")
