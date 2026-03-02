# =========================================================
# 10_Fig3_AC_Coherence_Final.R
# Figure 3 (legacy-style, restored)
#   - Fig 3A: A↔C Coherence Heatmap (NPMI & log2(Lift)) [2 panels]
#   - Fig 3B: A↔C Coherence Scatter (NPMI vs log2(Lift)), faceted by Outcome(C)
# Design follows the previous "Figure 3 (final tuning)" style.
# =========================================================

## ---- setup ----
# Use the project-level setup script (public/reproducible).
.find_setup <- function() {
  cand <- c(
    file.path(getwd(), "00_setup_Final.R"),
    Sys.getenv("RAILD_ROOT", unset=""),
    ""
  )
  if (cand[2] != "") cand[2] <- file.path(cand[2], "00_setup_Final.R")
  for (p in cand) if (nzchar(p) && file.exists(p)) return(p)
  stop("Cannot find 00_setup_Final.R. Open the project root (recommended) or set env var RAILD_ROOT to the project folder.")
}
source(.find_setup())


suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(forcats)
  library(scales); library(glue)
  library(systemfonts); library(showtext)
})

# ---- Font ----
cands  <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic","IPAexGothic","Noto Sans CJK JP")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

# ---- Input (tagged) ----
f_npmi <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_npmi)) stop("cooc_npmi_lift not found: ", f_npmi)

D0 <- readr::read_csv(f_npmi, show_col_types = FALSE)

need_cols <- c("A","C","pA","pC","pAC","lift","npmi")
mis <- setdiff(need_cols, names(D0))
if (length(mis)) stop("Required columns are missing in the CSV: ", paste(mis, collapse=", "))

# =========================================================
# Parameters (same spirit as legacy)
# =========================================================
MIN_pAC   <- 0       # e.g., set 0.002 for a relaxed filter (if needed)
TOP_PER_C <- 50      # How many top A-terms to show per C (rank by NPMI, tie-break by log2_lift)
LABEL_PER_C <- 8     # Number of labels per C in the scatter plot

# thresholds (same as legacy)
th_npmi <- 0.05
th_l2l  <- log2(1.20)  # lift>=1.2

# clip ranges (same as legacy)
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# =========================================================
# Preprocess
# =========================================================
D <- D0 %>%
  mutate(
    A = as.character(A),
    C = as.character(C),
    log2_lift = log2(ifelse(lift > 0, lift, NA_real_))
  ) %>%
  filter(!is.na(npmi), !is.na(log2_lift), pAC >= MIN_pAC)

# Select top terms per C (rank by NPMI, tie-break by log2_lift)
D <- D %>%
  group_by(C) %>%
  arrange(desc(npmi), desc(log2_lift), .by_group = TRUE) %>%
  mutate(rank_c = row_number()) %>%
  ungroup() %>%
  filter(rank_c <= TOP_PER_C)

# Global A ordering (by mean rank across C: key legacy behavior)
ord_A <- D %>%
  group_by(A) %>%
  summarise(mr = mean(rank_c, na.rm = TRUE), .groups="drop") %>%
  arrange(mr) %>%
  pull(A)
D <- D %>% mutate(A = factor(A, levels = ord_A))

# C ordering (by frequency, descending)
ord_C <- D %>% count(C, name="n") %>% arrange(desc(n)) %>% pull(C)
D <- D %>% mutate(C = factor(C, levels = ord_C))

# =========================================================
# Figure 3A: Heatmap (NPMI & log2(Lift))  [2 panels]
# =========================================================
D_heat_npmi <- D %>%
  mutate(
    npmi_clip = clip(npmi, -0.2, 1.0),
    metric = factor("NPMI", levels=c("NPMI","log2(Lift)")),
    value = npmi_clip
  )

D_heat_l2l <- D %>%
  mutate(
    l2l_clip = clip(log2_lift, -2, 2),
    metric = factor("log2(Lift)", levels=c("NPMI","log2(Lift)")),
    value = l2l_clip
  )

HH <- bind_rows(D_heat_npmi, D_heat_l2l)

p_heat <- ggplot(HH, aes(x = A, y = C, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = "#4575b4", mid = "#ffffbf", high = "#d73027",
    midpoint = 0, name = "Value"
  ) +
  labs(
    title = "Figure 3A. A↔C Coherence Heatmap",
    subtitle = paste0("NPMI & log2(Lift) | source: ", basename(f_npmi)),
    x = "A term", y = "Outcome (C)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = jpfont),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(face="bold")
  )

f3a_pdf <- file.path(DIR_FIG2, sprintf("Fig3A_AC_coherence_heatmap_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f3a_png <- file.path(DIR_FIG2, sprintf("Fig3A_AC_coherence_heatmap_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggsave(f3a_pdf, p_heat, device = cairo_pdf_device, width = 10, height = 8)
ggsave(f3a_png, p_heat, width = 10, height = 8, dpi = 300)
log_msg("WROTE:", f3a_pdf)

# =========================================================
# Figure 3B: Scatter (npmi vs log2(lift)), faceted by Outcome(C)
# =========================================================
lab_df <- D %>%
  group_by(C) %>%
  arrange(desc(npmi + 0.5*log2_lift), .by_group = TRUE) %>%
  slice_head(n = LABEL_PER_C) %>%
  ungroup()

p_sc <- ggplot(D, aes(x = npmi, y = log2_lift)) +
  geom_hline(yintercept = th_l2l, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = th_npmi, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = pAC, color = C), alpha = 0.8, show.legend = TRUE) +
  ggplot2::geom_text(
    data = lab_df,
    aes(label = as.character(A)),
    size = 3.2, family = jpfont,
    nudge_y = 0.03, check_overlap = TRUE
  ) +
  facet_wrap(~ C, scales = "free") +
  scale_size_continuous(
    name  = "p(A,C) (co-occur rate)",
    range = c(1, 5)
  ) +
  # Colors differ by outcome; set guide="none" to hide the legend if desired
  scale_color_discrete(name = "Outcome (C)") +
  labs(
    title = "Figure 3B. A↔C Coherence Scatter",
    subtitle = glue("Guides: NPMI ≥ {th_npmi}, Lift ≥ 1.2 | source: {basename(f_npmi)}"),
    x = "NPMI",
    y = "log2(Lift)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = jpfont),
    panel.grid = element_blank(),
    legend.position = "right"
  )

f3b_pdf <- file.path(DIR_FIG2, sprintf("Fig3B_AC_coherence_scatter_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f3b_png <- file.path(DIR_FIG2, sprintf("Fig3B_AC_coherence_scatter_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggsave(f3b_pdf, p_sc, device = cairo_pdf_device, width = 10, height = 7)
ggsave(f3b_png, p_sc, width = 10, height = 7, dpi = 300)
log_msg("WROTE:", f3b_pdf)

log_msg("=== DONE Fig3A/Fig3B (legacy restored) ===")


# =========================================================
# Figure 3B (subset): Outcomes limited to 4 (AE-ILD, hospitalization, mortality, progression)
# =========================================================

C_FOCUS4 <- c("AE-ILD", "hospitalization", "mortality", "progression")

D4 <- D %>%
  mutate(C = as.character(C)) %>%
  filter(C %in% C_FOCUS4) %>%
  mutate(C = factor(C, levels = C_FOCUS4))

# Labeling: only a few top points per C (readability first)
LABEL_PER_C4 <- 10
lab_df4 <- D4 %>%
  group_by(C) %>%
  arrange(desc(npmi + 0.5*log2_lift), .by_group = TRUE) %>%
  slice_head(n = LABEL_PER_C4) %>%
  ungroup()

p_sc4 <- ggplot(D4, aes(x = npmi, y = log2_lift)) +
  geom_hline(yintercept = th_l2l, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = th_npmi, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = pAC, color = C), alpha = 0.8, show.legend = TRUE) +
  ggplot2::geom_text(
    data = lab_df4,
    aes(label = as.character(A)),
    size = 3.2, family = jpfont,
    nudge_y = 0.03, check_overlap = TRUE
  ) +
  facet_wrap(~ C, scales = "free") +
  scale_size_continuous(
    name  = "p(A,C) (co-occur rate)",
    range = c(1, 5)
  ) +
  scale_color_discrete(name = "Outcome (C)") +
  labs(
    title = "Figure 3B (focused). A↔C Coherence Scatter (4 outcomes)",
    subtitle = glue("Outcomes: AE-ILD / hospitalization / mortality / progression | NPMI ≥ {th_npmi}, Lift ≥ 1.2"),
    x = "NPMI",
    y = "log2(Lift)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = jpfont),
    panel.grid = element_blank(),
    legend.position = "right"
  )

f3b4_pdf <- file.path(DIR_FIG2, sprintf("Fig3B_AC_coherence_scatter_focus4_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f3b4_png <- file.path(DIR_FIG2, sprintf("Fig3B_AC_coherence_scatter_focus4_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggsave(f3b4_pdf, p_sc4, device = cairo_pdf_device, width = 10, height = 7)
ggsave(f3b4_png, p_sc4, width = 10, height = 7, dpi = 300)
log_msg("WROTE:", f3b4_pdf)


# =========================================================
# Figure 3A (focused): Heatmap limited to 4 outcomes
#   AE-ILD / hospitalization / mortality / progression
# =========================================================

C_FOCUS4 <- c("AE-ILD", "hospitalization", "mortality", "progression")

# Extract only 4 outcomes from D (the processed data used for Fig3A/3B)
D4 <- D %>%
  mutate(C = as.character(C)) %>%
  filter(C %in% C_FOCUS4) %>%
  mutate(C = factor(C, levels = C_FOCUS4))

# Recompute A ordering within the 4-outcome subset (stable appearance)
ord_A4 <- D4 %>%
  group_by(A) %>%
  summarise(mr = mean(rank_c, na.rm = TRUE), .groups="drop") %>%
  arrange(mr) %>%
  pull(A)

D4 <- D4 %>% mutate(A = factor(A, levels = ord_A4))

# --- Heatmap: NPMI ---
D4_npmi <- D4 %>%
  mutate(
    npmi_clip = clip(npmi, -0.2, 1.0),
    metric = factor("NPMI", levels=c("NPMI","log2(Lift)")),
    value = npmi_clip
  )

# --- Heatmap: log2(Lift) ---
D4_l2l <- D4 %>%
  mutate(
    l2l_clip = clip(log2_lift, -2, 2),
    metric = factor("log2(Lift)", levels=c("NPMI","log2(Lift)")),
    value = l2l_clip
  )

HH4 <- bind_rows(D4_npmi, D4_l2l)

p_heat4 <- ggplot(HH4, aes(x = A, y = C, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = "#4575b4", mid = "#ffffbf", high = "#d73027",
    midpoint = 0, name = "Value"
  ) +
  labs(
    title = "Figure 3A (focused). A↔C Coherence Heatmap (4 outcomes)",
    subtitle = paste0("Outcomes: ", paste(C_FOCUS4, collapse=", "),
                      " | source: ", basename(f_npmi)),
    x = "A term", y = "Outcome (C)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = jpfont),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(face="bold")
  )

f3a4_pdf <- file.path(DIR_FIG2, sprintf("Fig3A_AC_coherence_heatmap_focus4_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f3a4_png <- file.path(DIR_FIG2, sprintf("Fig3A_AC_coherence_heatmap_focus4_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggsave(f3a4_pdf, p_heat4, device = cairo_pdf_device, width = 10, height = 6)
ggsave(f3a4_png, p_heat4, width = 10, height = 6, dpi = 300)
log_msg("WROTE:", f3a4_pdf)


# =========================================================
# [ADD-ON] Manuscript Figure 3: combine focus4 Fig3A + Fig3B into ONE page
#   - DO NOT modify the existing code above; just append this block.
#   - Uses p_heat4 (focused heatmap) and p_sc4 (focused scatter) already created above.
# =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# --- Safety checks ---
if (!exists("p_heat4")) stop("p_heat4 not found. Run the script above that creates p_heat4 first.")
if (!exists("p_sc4"))   stop("p_sc4 not found. Run the script above that creates p_sc4 first.")
if (!exists("DIR_FIG2")) stop("DIR_FIG2 not found. Ensure 00_setup_Final.R defines DIR_FIG2.")
if (!exists("CORPUS_TAG")) CORPUS_TAG <- "corpus"
if (!exists("DIC_TAG"))    DIC_TAG    <- "dic"

# --- 1) Fix overlapping A-term labels in Fig3A (focused heatmap) ---
p_heat4_pub <- p_heat4 +
  labs(
    title = "Figure 3A. A↔C coherence heatmap (4 outcomes)",
    subtitle = NULL
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),

    # key: reduce overlap for many A terms
    axis.text.x  = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),

    strip.text   = element_text(face = "bold", size = 11),
    panel.grid   = element_blank(),
    # extra bottom margin so 90-degree labels are not clipped
    plot.margin  = margin(t = 6, r = 8, b = 26, l = 8)
  )

# If guide_axis is available (ggplot2 >= 3.4.0), dodge labels into 2 rows for readability
if ("guide_axis" %in% getNamespaceExports("ggplot2")) {
  p_heat4_pub <- p_heat4_pub +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}

# --- 2) Make Fig3B (focused scatter) a bit more "paper-like" ---
p_sc4_pub <- p_sc4 +
  labs(
    title = "Figure 3B. A↔C coherence scatter (4 outcomes)",
    subtitle = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    strip.text   = element_text(face = "bold", size = 11),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10),
    panel.grid   = element_blank(),
    plot.margin  = margin(t = 6, r = 8, b = 6, l = 8)
  )

# --- 3) Combine into a single Figure 3 page (A on top, B below) ---
out_pdf <- file.path(DIR_FIG2, sprintf("Figure3_AB_focus4_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
out_png <- file.path(DIR_FIG2, sprintf("Figure3_AB_focus4_%s__%s.png", CORPUS_TAG, DIC_TAG))

# Prefer patchwork; fallback to cowplot if patchwork is unavailable
if (requireNamespace("patchwork", quietly = TRUE)) {

  suppressPackageStartupMessages(library(patchwork))

  p_comb <- (p_heat4_pub / p_sc4_pub) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.99)
    )

  ggsave(out_pdf, p_comb, device = cairo_pdf_device, width = 10, height = 12)
  ggsave(out_png, p_comb, width = 10, height = 12, dpi = 300)

} else if (requireNamespace("cowplot", quietly = TRUE)) {

  suppressPackageStartupMessages(library(cowplot))

  pA <- cowplot::ggdraw(p_heat4_pub) + cowplot::draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 16)
  pB <- cowplot::ggdraw(p_sc4_pub)   + cowplot::draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 16)

  p_comb <- cowplot::plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 1.15))

  ggsave(out_pdf, p_comb, device = cairo_pdf_device, width = 10, height = 12)
  ggsave(out_png, p_comb, width = 10, height = 12, dpi = 300)

} else {
  stop("Neither 'patchwork' nor 'cowplot' is available. Please install one: install.packages('patchwork') or install.packages('cowplot').")
}

log_msg("WROTE combined Figure 3:", out_pdf)
log_msg("WROTE combined Figure 3:", out_png)

