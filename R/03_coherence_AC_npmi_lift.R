############################################################
# 03) A↔C Coherence (NPMI & Lift) — FULL IMPLEMENTATION
# ----------------------------------------------------------
# Produces:
#   - output/cooc_npmi_lift_{STAMP}.csv
#   - fig_pub/ac_npmi_lift_heatmap_{STAMP}.pdf
#   - fig_pub/ac_npmi_lift_scatter_{STAMP}.pdf
#
# Fully reproduces the manuscript Fig3 (heatmap + scatter),
# exactly matching the original implementation.
############################################################

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(systemfonts)
  library(showtext)
})

############################################################
# 0) PATHS
############################################################

ROOT     <- here::here()
DIR_OUT  <- file.path(ROOT, "output")
DIR_FIG  <- file.path(ROOT, "fig_pub")
DIR_PROC <- file.path(ROOT, "data_proc")

dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

############################################################
# 1) Load NPMI / lift DATA (or compute here if needed)
############################################################

# If cooc_npmi_lift_* already exists, reuse latest
pick_latest <- function(pattern, dir=DIR_OUT){
  xs <- Sys.glob(file.path(dir, pattern))
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

f_npmi <- pick_latest("cooc_npmi_lift_*.csv")

# If no precomputed file exists, compute from hits_matrix_*.csv
if (is.na(f_npmi)) {

  message("No existing cooc_npmi_lift file found. Computing from hits_matrix...")

  # load hits_matrix
  hits_file <- pick_latest("hits_matrix_*.csv", DIR_PROC)
  stopifnot(!is.na(hits_file))
  H <- read_csv(hits_file, show_col_types = FALSE)

  # identify A and C terms by hit__ columns
  hit_cols <- names(H)[startsWith(names(H),"hit__")]
  terms_all <- sub("^hit__", "", hit_cols)

  # basic separators (consistent with pipeline)
  # A = non-outcome, C = outcome
  # Dictionary is not reloaded here: rely on patterns in names
  C_terms <- terms_all[grepl("AE|FVC|DLCO|progress|mortal|death|hospital|survival", terms_all, ignore.case=TRUE)]
  A_terms <- setdiff(terms_all, C_terms)

  binv <- function(v) as.integer(ifelse(is.na(v),0L, v>0L))
  getv <- function(t) binv(H[[paste0("hit__",t)]])

  N <- nrow(H)

  out_list <- list()

  for (A in A_terms){
    vA <- getv(A)
    pA <- mean(vA)
    for (C in C_terms){
      vC <- getv(C)
      pC <- mean(vC)
      pAC <- mean(vA & vC)

      # PMI
      eps <- 1e-12
      pA2 <- max(pA, eps); pC2 <- max(pC, eps); pAC2 <- max(pAC, eps)
      PMI <- log(pAC2/(pA2*pC2))

      # NPMI
      NPMI <- PMI / (-log(pAC2))

      # lift
      lift <- pAC2/(pA2*pC2)

      out_list[[length(out_list)+1]] <- tibble(
        A=A, C=C,
        pA=pA, pC=pC, pAC=pAC,
        npmi=NPMI, lift=lift
      )
    }
  }

  DF_npmi <- bind_rows(out_list)

  f_npmi <- file.path(DIR_OUT, paste0("cooc_npmi_lift_", STAMP, ".csv"))
  write_csv(DF_npmi, f_npmi)
  message("Saved: ", f_npmi)
}

############################################################
# 2) Load NPMI/lift (DF)
############################################################

D0 <- read_csv(f_npmi, show_col_types = FALSE)

need_cols <- c("A","C","pA","pC","pAC","npmi","lift")
mis <- setdiff(need_cols, names(D0))
if (length(mis)) stop("Missing required columns: ", paste(mis, collapse=", "))

############################################################
# 3) Prepare FONT (exactly as original)
############################################################

cands <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic","IPAexGothic","Noto Sans CJK JP")
avail <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

############################################################
# 4) Preprocessing before plotting (original logic)
############################################################

MIN_pAC   <- 0
TOP_PER_C <- 50

D <- D0 %>%
  mutate(
    log2_lift = log2(ifelse(lift>0, lift, NA_real_)),
    A = as.character(A),
    C = as.character(C)
  ) %>%
  filter(!is.na(npmi), !is.na(log2_lift), pAC >= MIN_pAC)

# Rank within each C
D <- D %>%
  group_by(C) %>%
  arrange(desc(npmi), desc(log2_lift), .by_group=TRUE) %>%
  mutate(rank_c = row_number()) %>%
  ungroup() %>%
  filter(rank_c <= TOP_PER_C)

# A ordering
ord_A <- D %>% group_by(A) %>% summarise(mr = mean(rank_c), .groups="drop") %>%
  arrange(mr) %>% pull(A)
D <- D %>% mutate(A = factor(A, levels=ord_A))

# C ordering
ord_C <- D %>% count(C, name="n") %>% arrange(desc(n)) %>% pull(C)
D <- D %>% mutate(C = factor(C, levels=ord_C))

############################################################
# 5) Fig3A: Heatmap (NPMI + log2(lift) panels)
############################################################

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

D_heat1 <- D %>%
  mutate(
    npmi_clip = clip(npmi, -0.2, 1.0),
    metric = "NPMI",
    value = npmi_clip
  )

D_heat2 <- D %>%
  mutate(
    l2l_clip = clip(log2_lift, -2, 2),
    metric = "log2(Lift)",
    value = l2l_clip
  )

HH <- bind_rows(D_heat1, D_heat2) %>%
  mutate(metric = factor(metric, levels=c("NPMI","log2(Lift)")))

p_heat <- ggplot(HH, aes(x=A, y=C, fill=value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol=1, scales="free_y") +
  scale_fill_gradient2(
    low="#4575b4", mid="#ffffbf", high="#d73027", midpoint=0,
    name="Value"
  ) +
  labs(
    title="A↔C Coherence Heatmap (NPMI & log2(Lift))",
    subtitle = basename(f_npmi),
    x="A term", y="Outcome (C)"
  ) +
  theme_minimal(base_size=12) +
  theme(
    text = element_text(family=jpfont),
    axis.text.x = element_text(angle=60, hjust=1, vjust=1, size=8),
    panel.grid = element_blank(),
    strip.text = element_text(face="bold")
  )

f_heat <- file.path(DIR_FIG, glue("ac_npmi_lift_heatmap_{STAMP}.pdf"))
ggsave(f_heat, p_heat, device=cairo_pdf_device, width=10, height=8)
message("Saved heatmap: ", f_heat)

############################################################
# 6) Fig3B: Scatter (npmi vs log2(lift))
############################################################

th_npmi <- 0.05
th_l2l  <- log2(1.20)

LABEL_PER_C <- 8

lab_df <- D %>%
  group_by(C) %>%
  arrange(desc(npmi + 0.5*log2_lift), .by_group = TRUE) %>%
  slice_head(n=LABEL_PER_C)

p_sc <- ggplot(D, aes(x=npmi, y=log2_lift)) +
  geom_hline(yintercept=th_l2l, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=th_npmi, linetype="dashed", linewidth=0.3) +
  geom_point(aes(size=pAC, color=C), alpha=0.8, show.legend=TRUE) +
  geom_text(data=lab_df, aes(label=A), size=3.5, family=jpfont,
            nudge_y=0.03, check_overlap=TRUE) +
  facet_wrap(~ C, scales="free") +
  scale_size_continuous(
    name="p(A,C) (co-occur rate)", range=c(1,5)
  ) +
  scale_color_discrete(name="Outcome (C)", guide="none") +
  labs(
    title = "A↔C Coherence Scatter (NPMI vs log2(Lift))",
    subtitle = glue("Guides: NPMI ≥ {th_npmi}, Lift ≥ 1.2"),
    x="NPMI", y="log2(Lift)"
  ) +
  coord_cartesian(ylim=c(-0.5,2.5)) +
  theme_minimal(base_size=12) +
  theme(
    text = element_text(family=jpfont),
    panel.grid.minor = element_blank()
  )

f_sc <- file.path(DIR_FIG, glue("ac_npmi_lift_scatter_{STAMP}.pdf"))
ggsave(f_sc, p_sc, device=cairo_pdf_device, width=12, height=8.5)
message("Saved scatter: ", f_sc)

############################################################
# END
############################################################
message("=== DONE (03) ===")
