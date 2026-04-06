# =========================================================
# 11_SuppFigS3_ABC_Final.R
# Supplementary Figure S3 (ABC block) — FINAL
#  - SuppFigS3A: Top triads barplot (A → B/B′ → C_main)
#  - SuppFigS3B: ABC network
#      * AB/BC edges separated by color + linetype (AB=solid, BC=longdash)
#      * edge width ∝ log1p(evidence count) to avoid domination by one edge
#      * node shapes denote roles (A/B/B′/C)
#  - Inputs: tagged outputs from 04_abc_rankings.R (+ optional abc_edges_top*)
#  - Dictionary: dic/ra_ild_dictionary_analysis_v1_genetics.csv (fixed)
# =========================================================

# ---- Load project setup (public-friendly) ----
# This script expects 00_setup_Final.R to define DIR_* paths and helper functions.
setup_candidates <- c(
  file.path(getwd(), "00_setup_Final.R"),
  file.path(Sys.getenv("RAILD_ROOT"), "00_setup_Final.R")
)
setup_path <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_path) || !nzchar(setup_path)) {
  stop("Cannot find 00_setup_Final.R. Run from the project root or set env var RAILD_ROOT to the project folder.")
}
source(setup_path)


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(systemfonts)
  library(showtext)
})

# ------------------ Font ------------------
cands  <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic","IPAexGothic","Noto Sans CJK JP")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

# ------------------ User controls ------------------
TOP_TRIADS  <- 60
FDR_THRES   <- 0.20
MIN_N11     <- 1

# Optional: cap the number of triads per A (set NA to disable).
CAP_PER_A <- NA_integer_   # e.g., 2 means "max 2 triads per same A"

# network layout: "fr" (recommended) or "kk"
NET_LAYOUT <- "fr"

# ------------------ Inputs (tagged) ------------------
f_rank <- file.path(DIR_TABLE, sprintf("abc_rankings_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_rank)) stop("abc_rankings not found: ", f_rank)
rank <- readr::read_csv(f_rank, show_col_types = FALSE)
log_msg("09 loaded rankings:", f_rank, " n=", nrow(rank))

# Prefer edges file from 04 if present
cand_edges <- Sys.glob(file.path(DIR_TABLE, sprintf("abc_edges_top*_%s__%s.csv", CORPUS_TAG, DIC_TAG)))
f_edges <- if (length(cand_edges)) cand_edges[which.max(file.info(cand_edges)$mtime)] else NA_character_

# Dictionary (fixed)
dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
dic <- readr::read_csv(dic_path, show_col_types = FALSE) %>%
  mutate(term = as.character(term),
         class = tolower(trimws(as.character(class))))

# ------------------ Role mapping ------------------
ROLE_A_CLASSES   <- c("drug","bio","gene","exposure","host","event","nonpharm","vaccine","population","system","strategy")
ROLE_B_CLASSES   <- c("biomarker","molecular","microbiome","cell","activity","pft","phenotype","trend","complication")
ROLE_BPR_CLASSES <- c("pattern","radiology","qct","airway")   # B′
ROLE_C_CLASSES   <- c("outcome")

class_lookup <- function(term_vec){
  out <- dic$class[match(term_vec, dic$term)]
  ifelse(is.na(out), NA_character_, out)
}
role_of_class <- function(cls){
  dplyr::case_when(
    cls %in% ROLE_C_CLASSES   ~ "C",
    cls %in% ROLE_BPR_CLASSES ~ "Bprime",
    cls %in% ROLE_B_CLASSES   ~ "B",
    cls %in% ROLE_A_CLASSES   ~ "A",
    TRUE ~ "Other"
  )
}
role_lookup <- function(term_vec){
  role_of_class(class_lookup(term_vec))
}

# ------------------ Pick top triads ------------------
score_col <- if ("score_q" %in% names(rank)) "score_q" else if ("score" %in% names(rank)) "score" else NULL
if (is.null(score_col)) stop("No score column in rankings.")
if (!all(c("A","B","C") %in% names(rank))) stop("rankings must have A/B/C columns.")

R <- rank
has_q <- all(c("AB_q","BC_q") %in% names(R))

R_try <- R
if (has_q) R_try <- R_try %>% filter(AB_q < FDR_THRES, BC_q < FDR_THRES)
if ("AB_n11" %in% names(R_try)) R_try <- R_try %>% filter(AB_n11 >= MIN_N11)
if ("BC_n11" %in% names(R_try)) R_try <- R_try %>% filter(BC_n11 >= MIN_N11)

R_try <- R_try %>% arrange(desc(.data[[score_col]]))
if (nrow(R_try) == 0) {
  log_msg("09 WARNING: triads empty after FDR/n11 filters; relaxing to score-only.")
  R_try <- R %>% arrange(desc(.data[[score_col]]))
}

Rf <- R_try %>% slice_head(n = min(TOP_TRIADS, nrow(R_try)))
stopifnot(nrow(Rf) > 0)

# Optionally cap per A
if (!is.na(CAP_PER_A)) {
  Rf <- Rf %>% group_by(A) %>% slice_head(n = CAP_PER_A) %>% ungroup() %>%
    arrange(desc(.data[[score_col]]))
}

C_main <- Rf$C[1]
log_msg("09 selected triads:", nrow(Rf), " | C_main=", C_main)

# =========================================================
# SuppFigS3A: Top triads barplot
# =========================================================
R4A <- Rf %>%
  mutate(
    roleB = role_lookup(B),
    B_disp = ifelse(roleB == "Bprime", paste0(B, " [B′]"), B),
    triad_label = paste0(A, " → ", B_disp, " → ", C)
  )

p4a <- ggplot(R4A, aes(x = .data[[score_col]], y = fct_rev(fct_inorder(triad_label)))) +
  geom_col() +
  labs(
    title = "Supplementary Figure S3A. Top A–B–C triads",
    subtitle = paste0("C = ", C_main, " | source: ", basename(f_rank),
                      if (has_q) paste0(" | FDR<", FDR_THRES) else ""),
    x = score_col,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(text = element_text(family = jpfont))

f4a_png <- file.path(DIR_FIG2, sprintf("SuppFigS3A_top_triads_%s__%s.png", CORPUS_TAG, DIC_TAG))
f4a_pdf <- file.path(DIR_FIG2, sprintf("SuppFigS3A_top_triads_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
ggsave(f4a_png, p4a, width = 11, height = max(6, 0.22*nrow(R4A) + 2), dpi = 300)
ggsave(f4a_pdf, p4a, width = 11, height = max(6, 0.22*nrow(R4A) + 2), device = cairo_pdf_device)
log_msg("WROTE:", f4a_pdf)

# =========================================================
# SuppFigS3B: ABC network
#  - AB/BC separated by color + linetype
#  - width uses log1p(count) for stability
# =========================================================

edges <- if (!is.na(f_edges) && file.exists(f_edges)) {
  log_msg("09 using edges:", f_edges)
  readr::read_csv(f_edges, show_col_types = FALSE)
} else {
  log_msg("09 WARNING: edges file not found; using rankings-derived fallback (n11).")
  bind_rows(
    R4A %>% transmute(from=A, to=B, w=as.integer(ifelse("AB_n11" %in% names(R4A), AB_n11, 1)), kind="AB"),
    R4A %>% transmute(from=B, to=C, w=as.integer(ifelse("BC_n11" %in% names(R4A), BC_n11, 1)), kind="BC")
  ) %>% group_by(from,to,kind) %>% summarise(w=max(w), .groups="drop")
}

# restrict to nodes involved in selected triads
keep_nodes <- unique(c(R4A$A, R4A$B, R4A$C))
edges <- edges %>% filter(from %in% keep_nodes, to %in% keep_nodes)

# node table (IMPORTANT: vertex name column must be 'name')
nodes <- tibble::tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(
    class = class_lookup(name),
    role  = role_of_class(class),
    role  = ifelse(name == C_main, "C", role),
    role2 = ifelse(is.na(role) | role == "Other", "Other", role)
  )

g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
E(g)$kind   <- edges$kind
E(g)$weight <- as.numeric(edges$w)
E(g)$w_plot <- log1p(E(g)$weight)   # stabilize edge width
V(g)$deg    <- igraph::degree(g, mode="all")

# colors / shapes
col_edge_AB <- "#e74c3c"  # muted red
col_edge_BC <- "#2471a3"  # muted blue

col_A      <- "#1b4f72"
col_B      <- "#117a65"
col_Bprime <- "#d35400"
col_C      <- "#b03a2e"
col_Other  <- "#7d7d7d"

shape_map <- c(A=22, B=21, Bprime=24, C=23, Other=25) # give B' a distinct shape
fill_map  <- c(A=col_A, B=col_B, Bprime=col_Bprime, C=col_C, Other=col_Other)

set.seed(42)
p4b <- ggraph::ggraph(g, layout = NET_LAYOUT) +
  ggraph::geom_edge_link(
    aes(edge_colour = kind, edge_width = w_plot, linetype = kind),
    alpha = 0.60,
    lineend = "round"
  ) +
  ggraph::geom_node_point(
    aes(shape = role2, size = deg, fill = role2),
    colour = "grey20", stroke = 0.7
  ) +
  ggraph::geom_node_text(
    aes(label = name),
    family = jpfont, size = 3, repel = TRUE
  ) +
  ggraph::scale_edge_colour_manual(
    values = c(AB = col_edge_AB, BC = col_edge_BC),
    labels = c(AB="AB (A → B/B′)", BC=paste0("BC (B/B′ → ", C_main, ")")),
    name = "Edge type"
  ) +
  # ★線種：論文記述に合わせる（AB solid / BC longdash）
  ggraph::scale_edge_linetype_manual(
    values = c(AB="solid", BC="longdash"),
    name = "Edge type"
  ) +
  ggraph::scale_edge_width_continuous(
    range = c(0.4, 2.5),
    name  = "Evidence (log1p count)"
  ) +
  scale_shape_manual(values = shape_map, name = "Node role") +
  scale_fill_manual(values = fill_map, guide = "none") +
  scale_size_continuous(range = c(3, 9), guide = "none") +
  labs(
    title = "Supplementary Figure S3B. ABC network (Top triads)",
    subtitle = paste0("AB: red solid | BC: blue dashed | width ∝ supporting records (log1p). ",
                      "rankings: ", basename(f_rank))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = jpfont),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

f4b_png <- file.path(DIR_FIG2, sprintf("SuppFigS3B_abc_network_%s__%s.png", CORPUS_TAG, DIC_TAG))
f4b_pdf <- file.path(DIR_FIG2, sprintf("SuppFigS3B_abc_network_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
ggsave(f4b_png, p4b, width = 10, height = 6, dpi = 300)
ggsave(f4b_pdf, p4b, width = 10, height = 6, device = cairo_pdf_device)

log_msg("WROTE:", f4b_pdf)
log_msg("=== DONE 09_SuppFigS3_ABC_FINAL_v2 ===")



# =========================================================
# ADD-ON (GitHub-ready): Create final MANUSCRIPT Supplementary Figure S3 (A+B) in ONE page
#   - SuppFigS3A: two-column barplot (Left: top, Right: next)
#   - SuppFigS3B: network (already created as p4b above)
#   - Output: SupplementaryFigureS3_AB_*.pdf/.png  (this is the manuscript figure to use)
# =========================================================

suppressPackageStartupMessages({
  library(patchwork)
})

# ---- SuppFigS3A (two-column, compact, muted academic color) ----
BAR_COL <- "#4E79A7"  # muted blue

R4A_2col <- R4A %>%
  mutate(.rank = row_number(),
         .half = ceiling(n() / 2),
         .panel = ifelse(.rank <= .half, "Left (Top)", "Right (Next)"),
         .panel = factor(.panel, levels = c("Left (Top)", "Right (Next)"))) %>%
  group_by(.panel) %>%
  mutate(.lab = forcats::fct_inorder(triad_label)) %>%
  ungroup()

p4a_2col <- ggplot(R4A_2col, aes(x = .data[[score_col]], y = forcats::fct_rev(.lab))) +
  geom_col(fill = BAR_COL) +
  facet_wrap(~ .panel, nrow = 1, scales = "free_y") +
  labs(
    title = "Supplementary Figure S3A. Top A–B–C triads",
    x = score_col,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = jpfont),
    plot.title   = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    strip.text   = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 8)
  )

# ---- SuppFigS3B (reuse p4b, slightly tighten subtitle to save space) ----
p4b_pub <- p4b +
  labs(subtitle = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.margin = margin(t = 6, r = 8, b = 6, l = 8)
  )

# ---- Combine (A thinner, B larger) ----
p4_comb <- (p4a_2col / p4b_pub) +
  plot_layout(heights = c(0.75, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)
  )

f4ab_pdf <- file.path(DIR_FIG2, sprintf("SupplementaryFigureS3_AB_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f4ab_png <- file.path(DIR_FIG2, sprintf("SupplementaryFigureS3_AB_%s__%s.png", CORPUS_TAG, DIC_TAG))

ggsave(f4ab_pdf, p4_comb, width = 12, height = 10, device = cairo_pdf_device)
ggsave(f4ab_png, p4_comb, width = 12, height = 10, dpi = 300)

log_msg("WROTE FINAL MANUSCRIPT SUPPLEMENTARY FIGURE S3:", f4ab_pdf)
