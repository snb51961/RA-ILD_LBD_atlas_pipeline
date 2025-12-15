############################################################
# 06) White Map + B-breakdown + Fig5C (AE-ILD)
# ----------------------------------------------------------
# Fully replicates:
#   - white_map_candidates_*.csv
#   - white_map tile (PDF)
#   - top_white_candidates_*.csv + barplot (PDF)
#   - B-breakdown (triad-based) per A→C (PDFs)
#   - B-breakdown fallback2 (hits_matrix-based) per A→C (PDFs)
#   - Fig5C: smoking / male_sex / MUC5B → AE-ILD (3-panel PDF)
#
# ROOT = here::here() for GitHub reproducibility
############################################################

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(patchwork)
  library(systemfonts)
  library(showtext)
})

# ---------------------------------------------------------
# 0) PATHS & FONTS
# ---------------------------------------------------------
ROOT     <- here::here()
DIR_PROC <- file.path(ROOT, "data_proc")
DIR_OUT  <- file.path(ROOT, "output")
DIR_FIG  <- file.path(ROOT, "fig")
DIR_FIGB <- file.path(DIR_FIG, "whitemap_B")

dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIGB, recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

cands <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic",
           "IPAexGothic","Noto Sans CJK JP")
avail <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

pick_latest <- function(pattern, dir){
  xs <- list.files(dir, pattern, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[order(file.info(xs)$mtime, decreasing = TRUE)][1]
}

# ---------------------------------------------------------
# 1) Load inputs
# ---------------------------------------------------------
f_hits <- pick_latest("^hits_matrix_.*\\.csv$", DIR_PROC)
stopifnot(!is.na(f_hits))
hits <- readr::read_csv(f_hits, show_col_types = FALSE)

f_abc <- pick_latest("^abc_rankings_.*\\.csv$", DIR_OUT)
stopifnot(!is.na(f_abc))
abc <- readr::read_csv(f_abc, show_col_types = FALSE)

hit_cols  <- names(hits)[startsWith(names(hits), "hit__")]
terms_all <- sub("^hit__", "", hit_cols)

binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))
get_hit <- function(df, t){
  col <- paste0("hit__", t)
  if (!(col %in% names(df))) return(integer(nrow(df)))
  binv(df[[col]])
}

# ---------------------------------------------------------
# 2) White map candidates
# ---------------------------------------------------------
C_terms <- unique(abc$C)

# ★重要：A は triad の A に揃え、hits_matrix にも存在するもののみ
A_terms <- intersect(unique(abc$A), terms_all)

ac_known <- lapply(A_terms, function(a){
  tibble::tibble(
    A = a,
    !!!setNames(lapply(C_terms, function(c){
      sum(get_hit(hits, a) & get_hit(hits, c))
    }), paste0("AC_n11__", C_terms))
  )
}) |> bind_rows()

ac_long <- tidyr::pivot_longer(
  ac_known,
  starts_with("AC_n11__"),
  names_to = "C",
  values_to = "AC_n11"
) %>%
  mutate(
    C = sub("^AC_n11__", "", C),
    AC_known = AC_n11 > 0
  )

abc_summ <- abc %>%
  group_by(A, C) %>%
  summarise(
    best_score_q = max(score_q, na.rm = TRUE),
    AB_evi_sum   = sum(AB_n11, na.rm = TRUE),
    BC_evi_sum   = sum(BC_n11, na.rm = TRUE),
    AB_evi_max   = max(AB_n11, na.rm = TRUE),
    BC_evi_max   = max(BC_n11, na.rm = TRUE),
    n_B          = n(),
    .groups = "drop"
  )

white_map <- ac_long %>%
  left_join(abc_summ, by = c("A","C")) %>%
  mutate(
    ABBC_strength = pmax(AB_evi_max, 0) + pmax(BC_evi_max, 0),
    priority      = dplyr::coalesce(best_score_q, 0) + 0.1 * ABBC_strength
  ) %>%
  arrange(desc(!AC_known), desc(priority))

out_cands <- file.path(DIR_OUT, paste0("white_map_candidates_", STAMP, ".csv"))
readr::write_csv(white_map, out_cands)
message("Saved: ", basename(out_cands))

# ---------------------------------------------------------
# 3) White map tile plot (PDF)
# ---------------------------------------------------------
plot_df <- white_map %>%
  mutate(status = ifelse(AC_known, "known A–C", "WHITE (A–C unknown)"))

p_white <- ggplot(plot_df, aes(x=C, y=A, fill=best_score_q)) +
  geom_tile(color="grey85") +
  geom_point(aes(size=ABBC_strength, shape=status), alpha=0.85) +
  scale_shape_manual(values=c("WHITE (A–C unknown)"=16, "known A–C"=1)) +
  guides(fill=guide_colorbar(title="best score_q"),
         size=guide_legend(title="AB/BC strength"),
         shape=guide_legend(title="A–C status")) +
  labs(title="White map (RA-ILD): A–C 未既知と AB/BC 根拠の俯瞰",
       x="Outcome (C)", y="Candidate factor (A)") +
  theme_minimal(base_size=11) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(family = jpfont))

ggsave(
  file.path(DIR_FIG, paste0("white_map_AC_tiles_", STAMP, ".pdf")),
  p_white,
  width = 10,
  height = 0.25*length(unique(plot_df$A)) + 3,
  device = cairo_pdf_device
)

# ---------------------------------------------------------
# 4) Top white candidates + bar plot
# ---------------------------------------------------------
cand_unknown <- white_map %>% filter(!AC_known)
if (nrow(cand_unknown) == 0) cand_unknown <- white_map

cand_ranked <- cand_unknown %>%
  mutate(ABBC_sum = AB_evi_sum + BC_evi_sum) %>%
  arrange(desc(priority), desc(best_score_q),
          desc(ABBC_strength), desc(ABBC_sum), desc(n_B))

topN  <- min(10, nrow(cand_ranked))
top10 <- cand_ranked %>% slice_head(n = topN)

out_top10 <- file.path(DIR_OUT, paste0("top_white_candidates_", STAMP, ".csv"))
readr::write_csv(
  top10 %>% select(A, C, priority, best_score_q, ABBC_strength,
                   AB_evi_max, BC_evi_max, AB_evi_sum, BC_evi_sum, n_B),
  out_top10
)
message("Saved: ", basename(out_top10))

p_top10 <- ggplot(top10, aes(x = priority,
                             y = reorder(paste0(A, " → ", C), priority))) +
  geom_col(width = 0.65) +
  labs(title = "Top White Map Candidates (A–C unknown)",
       x = "Priority (score blend)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = jpfont))

ggsave(
  file.path(DIR_FIG, paste0("top_white_candidates_bar_", STAMP, ".pdf")),
  p_top10,
  width = 7,
  height = max(3, 0.4 * nrow(top10) + 1),
  device = cairo_pdf_device
)

# ---------------------------------------------------------
# 5) B-breakdown (triad-based; mortality excluded from plots)
# ---------------------------------------------------------
pairs <- top10 %>% distinct(A, C)

detail <- abc %>%
  inner_join(pairs, by = c("A","C")) %>%
  group_by(A, C, B) %>%
  summarise(
    n_ev_AB      = sum(AB_n11, na.rm = TRUE),
    n_ev_BC      = sum(BC_n11, na.rm = TRUE),
    best_AB      = max(AB_n11, na.rm = TRUE),
    best_BC      = max(BC_n11, na.rm = TRUE),
    best_score_q = max(score_q, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  mutate(ABBC_strength = pmax(best_AB,0) + pmax(best_BC,0)) %>%
  arrange(A, C, desc(best_score_q), desc(ABBC_strength), desc(n_ev_AB + n_ev_BC))

k <- 6
detail_top <- detail %>% group_by(A,C) %>% slice_head(n = k) %>% ungroup()

out_Btab <- file.path(DIR_OUT, paste0("whitemap_B_breakdown_top", k, "_", STAMP, ".csv"))
readr::write_csv(detail_top, out_Btab)
message("Saved: ", basename(out_Btab))

for (i in seq_len(nrow(pairs))) {
  A0 <- pairs$A[i]; C0 <- pairs$C[i]
  D  <- detail_top %>%
    filter(A == A0, C == C0) %>%
    filter(B != "mortality") %>%
    arrange(best_score_q)
  if (nrow(D) == 0) next

  gp <- ggplot(D, aes(x = best_score_q, y = factor(B, levels = D$B))) +
    geom_segment(aes(x = 0, xend = best_score_q, yend = B), linewidth = 0.6) +
    geom_point(aes(size = ABBC_strength), alpha = 0.9) +
    labs(title = paste0("B breakdown: ", A0, " → ", C0),
         x = "best_score_q", y = "B (intermediate)") +
    theme_minimal(base_size = 11) +
    theme(text = element_text(family = jpfont))

  ggsave(
    file.path(DIR_FIGB,
              paste0("B_breakdown_",
                     gsub("[^A-Za-z0-9]+","_", A0),
                     "_to_",
                     gsub("[^A-Za-z0-9]+","_", C0),
                     "_", STAMP, ".pdf")),
    gp,
    width  = 5,
    height = max(3, 0.5*nrow(D) + 1),
    device = cairo_pdf_device
  )
}

# ---------------------------------------------------------
# 6) fallback2 (hits_matrix-based; mortality excluded from plots)
# ---------------------------------------------------------
for (i in seq_len(nrow(pairs))) {
  A0 <- pairs$A[i]; C0 <- pairs$C[i]
  message("Running fallback2 for: ", A0, " → ", C0)

  col_A <- paste0("hit__", A0)
  col_C <- paste0("hit__", C0)
  if (!all(c(col_A, col_C) %in% names(hits))) next

  vA <- binv(hits[[col_A]])
  vC <- binv(hits[[col_C]])

  B_pool <- unique(abc$B)
  B_pool <- B_pool[nzchar(B_pool)]
  B_cols <- paste0("hit__", B_pool)
  B_cols <- B_cols[B_cols %in% names(hits)]
  if (!length(B_cols)) {
    B_cols <- setdiff(hit_cols, c(col_A, col_C))
  }

  AB_counts <- sapply(B_cols, function(bc) sum(binv(hits[[bc]]) & vA))
  BC_counts <- sapply(B_cols, function(bc) sum(binv(hits[[bc]]) & vC))

  DF <- tibble::tibble(
    B      = sub("^hit__", "", B_cols),
    AB_n11 = as.integer(AB_counts),
    BC_n11 = as.integer(BC_counts)
  ) %>%
    mutate(ABBC_strength = AB_n11 + BC_n11) %>%
    filter(ABBC_strength > 0) %>%
    filter(B != "mortality") %>%
    arrange(desc(ABBC_strength), desc(AB_n11), desc(BC_n11)) %>%
    slice_head(n = k)

  save_path <- file.path(DIR_FIGB,
    paste0("B_breakdown_", gsub("[^A-Za-z0-9]+","_", A0), "_to_",
           gsub("[^A-Za-z0-9]+","_", C0), "_", STAMP, "_fallback2.pdf"))

  if (nrow(DF) > 0) {
    gp <- ggplot(DF, aes(x = ABBC_strength, y = factor(B, levels = DF$B))) +
      geom_segment(aes(x = 0, xend = ABBC_strength, yend = B), linewidth = 0.6) +
      geom_point(aes(size = ABBC_strength), alpha = 0.9) +
      labs(title = paste0("B breakdown: ", A0, " → ", C0, " (fallback2)"),
           x = "AB_n11 + BC_n11 (co-occurrence)", y = "B (intermediate)") +
      theme_minimal(base_size = 11) +
      theme(text = element_text(family = jpfont))

    ggsave(save_path, gp,
           width = 5, height = max(3, 0.5*nrow(DF) + 1),
           device = cairo_pdf_device)
  } else {
    gp_empty <- ggplot() +
      annotate("text", x=0, y=0,
               label=paste0("No B found even from hits_matrix for ", A0, " → ", C0),
               size=4) +
      theme_void()
    ggsave(sub("\\.pdf$", "_EMPTY.pdf", save_path),
           gp_empty, width=5, height=3, device=cairo_pdf_device)
  }
}

# ---------------------------------------------------------
# 7) Fig5C: selected A → AE-ILD (fallback2) 3-panel
# ---------------------------------------------------------
A_targets <- c("smoking","male_sex","MUC5B")
C0 <- "AE-ILD"

get_B_fallback2 <- function(A0, C0, k = 6){
  col_A <- paste0("hit__", A0)
  col_C <- paste0("hit__", C0)
  if (!all(c(col_A, col_C) %in% names(hits))) return(NULL)

  vA <- binv(hits[[col_A]])
  vC <- binv(hits[[col_C]])

  B_pool <- unique(abc$B)
  B_pool <- B_pool[nzchar(B_pool)]
  B_cols <- paste0("hit__", B_pool)
  B_cols <- B_cols[B_cols %in% names(hits)]
  if (!length(B_cols)) {
    B_cols <- setdiff(hit_cols, c(col_A, col_C))
  }

  AB_counts <- sapply(B_cols, function(bc) sum(binv(hits[[bc]]) & vA))
  BC_counts <- sapply(B_cols, function(bc) sum(binv(hits[[bc]]) & vC))

  DF <- tibble::tibble(
    A      = A0,
    C      = C0,
    B      = sub("^hit__", "", B_cols),
    AB_n11 = as.integer(AB_counts),
    BC_n11 = as.integer(BC_counts)
  ) %>%
    mutate(ABBC_strength = AB_n11 + BC_n11) %>%
    filter(ABBC_strength > 0) %>%
    filter(B != "mortality") %>%
    arrange(desc(ABBC_strength), desc(AB_n11), desc(BC_n11)) %>%
    slice_head(n = k)

  if (!nrow(DF)) return(NULL)
  DF
}

df_list <- lapply(A_targets, get_B_fallback2, C0 = C0, k = 6)
df_combined <- dplyr::bind_rows(df_list)

if (!nrow(df_combined)) {
  message("No B-breakdown data for selected A → AE-ILD; Fig5C not created.")
} else {
  df_combined <- df_combined %>%
    group_by(A) %>%
    arrange(ABBC_strength, .by_group = TRUE) %>%
    mutate(B_f = factor(B, levels = unique(B))) %>%
    ungroup()

  p_5C <- ggplot(df_combined, aes(x = ABBC_strength, y = B_f)) +
    geom_segment(aes(x = 0, xend = ABBC_strength, yend = B_f,
                     colour = ABBC_strength),
                 linewidth = 0.6) +
    geom_point(aes(size = ABBC_strength, fill = ABBC_strength),
               shape = 21, color = "black", alpha = 0.9) +
    scale_size_continuous(
      name  = "Bubble size\n#(A∧B) + #(B∧C)\n(C = AE-ILD)",
      range = c(2.5, 6)
    ) +
    scale_fill_gradient(
      name = "ABBC_strength",
      low  = "#deebf7",
      high = "#08519c"
    ) +
    scale_colour_gradient(
      low = "#deebf7",
      high = "#08519c",
      guide = "none"
    ) +
    facet_wrap(~ A, nrow = 1, scales = "free_y") +
    labs(
      title = "B breakdown for selected A → AE-ILD (fallback2)",
      x     = "Co-occurrence count\n[#(A∧B) + #(B∧C)]  (C = AE-ILD)",
      y     = "Intermediate B"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      text            = element_text(family = jpfont),
      strip.text      = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.title.y    = element_text(margin = margin(r = 4)),
      axis.title.x    = element_text(margin = margin(t = 4))
    )

  out_5C <- file.path(DIR_FIGB, paste0("B_breakdown_Fig5C_selected_AE-ILD_", STAMP, ".pdf"))
  ggsave(out_5C, p_5C, device = cairo_pdf_device, width = 10, height = 3.5)
  message("Saved Fig5C candidate: ", out_5C)
}

message("=== DONE (06) ===")
