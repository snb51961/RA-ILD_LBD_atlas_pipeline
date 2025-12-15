############################################################
# 05) Signed Effects vs AE-ILD (Fig5A + Fig5B)
# ----------------------------------------------------------
# Inputs:
#   output/signed_effects_summary_withCI_{STAMP}.csv
#   output/abc_rankings_{STAMP}.csv
#
# Outputs:
#   fig_pub/signed_effects_AE-ILD_{STAMP}.pdf        (Fig5A)
#   fig_pub/signed_effects_AE-ILD_bubble_{STAMP}.pdf (old)
#   fig_pub/signed_effects_AE-ILD_triadLBD_{STAMP}.pdf (Fig5B)
#
# Fully faithful to the implementation in:
#   "Figure調整用追加コード.rtf"
############################################################

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(systemfonts)
  library(showtext)
})

############################################################
# 0) PATHS & FONT
############################################################

ROOT     <- here::here()
DIR_OUT  <- file.path(ROOT, "output")
DIR_FIG  <- file.path(ROOT, "fig_pub")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# Japanese font
cands <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic",
           "IPAexGothic","Noto Sans CJK JP")
avail <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

############################################################
# 1) Load latest signed-effects summary (with CI)
############################################################

files_se <- list.files(
  DIR_OUT,
  pattern="^signed_effects_summary_withCI_\\d{8}_\\d{4}\\.csv$",
  full.names=TRUE
)
stopifnot(length(files_se) > 0)

f_csv <- files_se[ which.max(file.info(files_se)$mtime) ]

# Extract timestamp from filename
STAMP <- sub("^.*signed_effects_summary_withCI_(\\d{8}_\\d{4})\\.csv$",
             "\\1", basename(f_csv))

DF <- read_csv(f_csv, show_col_types = FALSE)

############################################################
# 2) Filter AE-ILD + Prepare DataFrame
############################################################

target_C <- "AE-ILD"

dat <- DF %>%
  filter(C == target_C) %>%
  mutate(
    articles     = coalesce(articles, 0),
    balance      = coalesce(balance, 0),
    balance_low  = coalesce(balance_low, NA_real_),
    balance_high = coalesce(balance_high, NA_real_),
    pos_ratio    = coalesce(pos_ratio, NA_real_)
  )

stopifnot(nrow(dat) > 0)

# Order A: articles desc → balance desc
dat <- dat %>%
  arrange(desc(articles), desc(balance)) %>%
  mutate(A_f = forcats::fct_inorder(A))

# Color scale (exact implementation)
col_fun <- scales::col_numeric(
  palette = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
              "#FDDBC7","#F4A582","#D6604D","#B2182B"),
  domain  = c(-1,1),
  na.color = "grey80"
)

############################################################
# 3) Fig5A — Forest Plot
############################################################

p_forest <- ggplot(dat, aes(x = balance, y = A_f)) +
  geom_errorbarh(
    aes(xmin = balance_low, xmax = balance_high),
    height = 0, alpha = 0.7, na.rm = TRUE
  ) +
  geom_point(
    aes(size = articles, fill = balance),
    shape = 21, stroke = 0.3, color = "black"
  ) +
  scale_fill_gradientn(
    colours = col_fun(seq(-1,1,length.out=9)),
    limits  = c(-1,1),
    oob     = scales::squish
  ) +
  scale_size_continuous(
    name   = "Articles",
    range  = c(2.3,7),
    breaks = scales::pretty_breaks(3)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  coord_cartesian(xlim = c(-1,1)) +
  labs(
    title    = "Signed effect vs AE-ILD (forest)",
    subtitle = paste0("Balance with Wilson CI | ", target_C, " | ", STAMP),
    x        = "Balance (risk_down ← 0 → risk_up)",
    y        = NULL,
    fill     = "Balance"
  ) +
  theme_minimal(base_family = jpfont, base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

outfile_forest <- file.path(DIR_FIG, paste0("signed_effects_", target_C, "_", STAMP, ".pdf"))
ggsave(outfile_forest, p_forest, device=cairo_pdf_device, width=5, height=5)
message("Saved Fig5A: ", outfile_forest)

############################################################
# 4) Old bubble (pos_ratio) — Supplement version
############################################################

dat_bub <- dat %>% filter(!is.na(pos_ratio))

p_bubble <- ggplot(dat_bub, aes(x = pos_ratio, y = A_f)) +
  geom_point(
    aes(size = articles, fill = balance),
    shape = 21, color = "black", stroke = 0.3, alpha = 0.9
  ) +
  scale_fill_gradientn(
    colours = col_fun(seq(-1,1,length.out=9)),
    limits  = c(-1,1),
    oob     = scales::squish
  ) +
  scale_size_continuous(
    name   = "Articles",
    range  = c(2.3,7),
    breaks = scales::pretty_breaks(3)
  ) +
  scale_x_continuous(
    limits = c(0,1),
    labels = scales::percent_format(accuracy=1),
    breaks = seq(0,1,0.2)
  ) +
  labs(
    title    = "Signed effect vs AE-ILD (bubble)",
    subtitle = paste0("x = Positive-articles ratio | ", target_C, " | ", STAMP),
    x        = "Positive articles ratio",
    y        = NULL,
    fill     = "Balance"
  ) +
  theme_minimal(base_family = jpfont, base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

outfile_bubble <- file.path(DIR_FIG, paste0("signed_effects_", target_C, "_bubble_", STAMP, ".pdf"))
ggsave(outfile_bubble, p_bubble, device=cairo_pdf_device, width=5, height=5)
message("Saved old bubble: ", outfile_bubble)

############################################################
# 5) Fig5B — Triad-based LBD support bubble
############################################################

# Load ABC rankings
files_abc <- list.files(DIR_OUT, pattern="^abc_rankings_\\d{8}_\\d{4}\\.csv$", full.names=TRUE)
stopifnot(length(files_abc) > 0)
f_abc <- files_abc[ which.max(file.info(files_abc)$mtime) ]
abc   <- read_csv(f_abc, show_col_types = FALSE)

abc_ae <- abc %>% filter(C == target_C)

triad_A <- abc_ae %>%
  group_by(A) %>%
  summarise(
    triad_strength = sum(score_q, na.rm=TRUE),
    max_score_q    = max(score_q, na.rm=TRUE),
    n_triads       = n(),
    .groups="drop"
  )

dat_triad <- dat %>%
  left_join(triad_A, by="A") %>%
  filter(!is.na(triad_strength)) %>%
  mutate(triad_strength_log = log10(triad_strength + 1e-6))

p_triad <- ggplot(dat_triad, aes(x = triad_strength_log, y = A_f)) +
  geom_point(
    aes(size = articles, fill = balance),
    shape=21, color="black", stroke=0.3, alpha=0.9
  ) +
  scale_fill_gradientn(
    colours = col_fun(seq(-1,1,length.out=9)),
    limits=c(-1,1),
    oob=scales::squish
  ) +
  scale_size_continuous(
    name="Articles",
    range=c(2.3,7),
    breaks=scales::pretty_breaks(3)
  ) +
  labs(
    title="Triad-based LBD support vs AE-ILD (bubble)",
    subtitle=paste0(
      "x = log10 Σ_B score_q(A,B,", target_C, ") | ", STAMP
    ),
    x="Triad-based LBD support (log10 Σ score_q)",
    y=NULL,
    fill="Balance"
  ) +
  theme_minimal(base_family=jpfont, base_size=12) +
  theme(
    legend.position="right",
    panel.grid.minor=element_blank()
  )

outfile_triad <- file.path(DIR_FIG, paste0("signed_effects_", target_C, "_triadLBD_", STAMP, ".pdf"))
ggsave(outfile_triad, p_triad, device=cairo_pdf_device, width=5, height=5)
message("Saved Fig5B: ", outfile_triad)

############################################################
# END
############################################################
message("=== DONE (05) ===")
