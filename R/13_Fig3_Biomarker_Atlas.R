# =========================================================
# 13_Fig3_Biomarker_Atlas.R
# Figure 3 = Biomarker evidence atlas (RA-ILD) — fully compatible + dense Panel B (fixed outcomes)
#
# Fixes vs denseB.R:
#  1) Panel B no longer uses a pre-filtered A candidate set. It builds the network from the full
#     biomarker-filtered table (S0) so edges exist for all outcomes.
#  2) Outcome nodes are ALWAYS added to the graph vertices even if an outcome has few edges.
#  3) Density is controlled ONLY by TOP_PER_A / TOP_PER_C / TOP_EDGES / MIN_EDGE_EVIDENCE.
# =========================================================

# --- Load shared setup (robust root handling) ---
.find_setup <- function(fname = "00_setup_Final.R") {
  # 1) current working directory
  p1 <- file.path(getwd(), fname)
  if (file.exists(p1)) return(p1)
  # 2) search up to 5 parent directories
  cur <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:5) {
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    p <- file.path(parent, fname)
    if (file.exists(p)) return(p)
    cur <- parent
  }
  # 3) environment variable override
  root_env <- Sys.getenv("RAILD_ROOT", unset = NA_character_)
  if (!is.na(root_env) && nzchar(root_env)) {
    p3 <- file.path(root_env, fname)
    if (file.exists(p3)) return(p3)
  }
  stop("Cannot find ", fname, ". Run from the project root or set RAILD_ROOT to the project folder.")
}

source(.find_setup())

# --- Minimal run log for reproducibility ---
stamp_run <- format(Sys.time(), "%Y%m%d-%H%M%S")
log_file <- file.path(DIR_LOG, paste0("run13_fig3_biomarker_atlas_", stamp_run, ".txt"))
try({
  dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)
  con <- file(log_file, open = "wt", encoding = "UTF-8")
  writeLines(c(
    paste0("script: ", "13_Fig3_Biomarker_Atlas.R"),
    paste0("run_at: ", Sys.time()),
    paste0("CORPUS_TAG: ", CORPUS_TAG),
    paste0("DIC_TAG: ", DIC_TAG),
    paste0("DIC_FILE: ", DIC_FILE),
    paste0("DIR_TABLE: ", DIR_TABLE),
    paste0("DIR_FIG3: ", if (exists("DIR_FIG3")) DIR_FIG3 else if (exists("DIR_FIG2")) DIR_FIG2 else if (exists("DIR_FIG")) DIR_FIG else "NA"),
    ""
  ), con)
  # Capture key inputs if present
  key_inputs <- c(
    file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    file.path(DIR_TABLE, sprintf("abc_rankings_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  )
  for (p in key_inputs) {
    if (file.exists(p)) {
      info <- file.info(p)
      writeLines(paste0("input: ", p, " | size=", info$size, " | mtime=", info$mtime), con)
    } else {
      writeLines(paste0("input_missing: ", p), con)
    }
  }
  writeLines("", con)
  writeLines("sessionInfo():", con)
  capture.output(sessionInfo(), file = con)
  close(con)
}, silent = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(igraph)
  library(ggraph)
library(ggrepel)
})



canon_outcome <- function(x){
  xl <- stringr::str_to_lower(stringr::str_replace_all(x, "–", "-"))
  dplyr::case_when(
    stringr::str_detect(xl, "ae-ild|aeild|acute exacerb|exacerb") ~ "AE-ILD",
    stringr::str_detect(xl, "progress|decline|fvc|dlco|worsen")   ~ "progression",
    stringr::str_detect(xl, "mortality|death|surviv")             ~ "mortality",
    TRUE ~ NA_character_
  )
}

# -------------------------
# Parameters (edit here)
# -------------------------
OUTCOMES <- c("AE-ILD","progression","mortality")
FOCUS_OUTCOME <- "progression"

# Core B-side classes (legacy)
B_CLASSES <- c("biomarker","molecular","cell","microbiome")

# Meta-B terms
META_TERMS <- c("Cytokine_inflammation","Fibrosis_pathway","Oxidative_stress_pathway","Immune_checkpoint")

# Always keep
MUST_KEEP <- c("KL6","SP_D","ACPA_RF","CRP_ESR","MMP7","MPO","NLR","CAR","CXCL9","IL17A","APRIL")

# Panel A
N_A_SHOW <- 45

# Panel B density controls
MIN_EDGE_EVIDENCE <- 0.30   # 0.30 drops single-hit edges (pos+neg=1)
TOP_PER_A <- 4              # per biomarker/family
TOP_PER_C <- 20             # per outcome
TOP_EDGES <- 300           # global cap (set NA to disable)
KK_SEED <- 42

# Panel C
MAX_LABELS_C <- 30

# -------------------------
# Helpers
# -------------------------
is_family <- function(x) str_detect(x, "_family$")
var_to_parent <- function(x, DIC=NULL) sub("_var$", "", x)

pick_dic <- function(){
  if (exists("DIC_FILE") && is.character(DIC_FILE) && length(DIC_FILE)==1 && file.exists(DIC_FILE)) return(DIC_FILE)
  p1 <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
  if (file.exists(p1)) return(p1)
  cand <- Sys.glob(file.path(DIR_DIC, paste0("*", DIC_TAG, "*.csv")))
  if (length(cand)) {
    cand <- cand[order(file.info(cand)$mtime, decreasing=TRUE)]
    return(cand[1])
  }
  cand <- Sys.glob(file.path(DIR_DIC, "*.csv"))
  if (!length(cand)) return(NA_character_)
  cand[which.max(file.info(cand)$mtime)]
}

build_bio_terms <- function(DIC, b_classes=B_CLASSES){
  term <- as.character(DIC$term)
  class <- tolower(as.character(DIC$class))
  if ("role" %in% names(DIC)) {
    role <- toupper(as.character(DIC$role))
    keep <- (role %in% c("B","BPRIME")) & (class %in% b_classes)
  } else {
    keep <- class %in% b_classes
  }
  unique(term[keep])
}

# -------------------------
# Load dictionary
# -------------------------
dic_path <- pick_dic()
if (is.na(dic_path) || !file.exists(dic_path)) stop("Dictionary CSV not found in DIR_DIC.")

DIC <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","class") %in% names(DIC))) stop("Dictionary must include columns term and class.")

message("FIG3 (denseB_v2) dictionary: ", dic_path,
        "  n=", nrow(DIC), "  has role? ", ("role" %in% names(DIC)))

BIO_TERMS <- build_bio_terms(DIC, b_classes=B_CLASSES)
BIO_TERMS <- unique(c(BIO_TERMS, META_TERMS, MUST_KEEP))

# -------------------------
# Load signed summary (06 output)
# -------------------------
f_sum <- file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_sum)) stop("Missing: ", f_sum, " (run 06_signed_effects.R first)")

S <- readr::read_csv(f_sum, show_col_types = FALSE)
need <- c("A","C","pos_articles","neg_articles","articles","balance")
miss <- setdiff(need, names(S))
if (length(miss)) stop("summary_withCI missing columns: ", paste(miss, collapse=", "))

S0 <- S %>%
  mutate(
    A        = as.character(A),
    C_raw    = as.character(C),
    C        = canon_outcome(C_raw),              # normalize first
    A_parent = var_to_parent(A, DIC),
    articles     = as.numeric(articles),
    pos_articles = as.numeric(pos_articles),
    neg_articles = as.numeric(neg_articles),
    balance      = as.numeric(balance),
    evidence     = log10(pmax(pos_articles + neg_articles, 1)),
    score        = (abs(balance) + 0.20) * (evidence + 0.50)
  ) %>%
  filter(!is.na(C)) %>%                            # keep only outcomes that can be normalized
  mutate(C = factor(C, levels = OUTCOMES)) %>%     # align to OUTCOMES ordering
  filter(is_family(A_parent) | (A_parent %in% BIO_TERMS))


if (nrow(S0) == 0) stop("After biomarker-atlas filtering, no terms remained.")

message("C raw unique: ", paste(unique(S$C), collapse=", "))
message("C canon table:")
print(table(S0$C))

# -------------------------
# Panel A (compact)
# -------------------------
A_rank <- S0 %>%
  group_by(A_parent) %>%
  summarise(scoreA = max(score, na.rm=TRUE), evA = sum(evidence, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(scoreA), desc(evA))

keepA_show <- unique(c(intersect(MUST_KEEP, A_rank$A_parent),
                       intersect(META_TERMS, A_rank$A_parent),
                       head(A_rank$A_parent, N_A_SHOW)))

S_A <- S0 %>% filter(A_parent %in% keepA_show)

A_order <- S_A %>%
  group_by(A_parent) %>%
  summarise(ev = sum(evidence, na.rm=TRUE), .groups="drop") %>%
  arrange(ev) %>%
  pull(A_parent)

S_A <- S_A %>% mutate(A_f = factor(A_parent, levels = A_order))

pA <- ggplot(S_A, aes(x = C, y = A_f)) +
  geom_point(aes(size = evidence, fill = balance),
             shape = 21, colour="grey25", stroke=0.3, alpha=0.95) +
  scale_fill_gradient2(low="steelblue", mid="white", high="firebrick",
                       midpoint=0, limits=c(-1,1), oob=scales::squish) +
  scale_size_continuous(range = c(1.2, 8.0)) +
  labs(
    title = "Biomarker evidence atlas – direction vs evidence (RA-ILD)",
    x = "RA-ILD outcomes (C)",
    y = "Biomarker / Family (B-role)",
    fill = "Directional balance\n(-1=protective, +1=risk)",
    size = "Evidence volume\nlog10(pos+neg)"
  ) +
  theme_minimal(base_size=10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=7),
    legend.position = "right",
    plot.title = element_text(face="bold", size=12)
  )

# -------------------------
# Panel B (dense): build from full S0, then legacy selection with balancing
# -------------------------
AtlasB <- S0 %>%
  filter(evidence >= MIN_EDGE_EVIDENCE) %>%
  mutate(A = A_parent)

# ensure each outcome has candidates
if (any(table(AtlasB$C) == 0)) {
  message("Warning: some outcomes have 0 candidate edges after MIN_EDGE_EVIDENCE filter.")
}

G_A <- AtlasB %>% group_by(A) %>% slice_max(order_by = score, n = TOP_PER_A, with_ties = FALSE) %>% ungroup()
G_C <- AtlasB %>% group_by(C) %>% slice_max(order_by = score, n = TOP_PER_C, with_ties = FALSE) %>% ungroup()

Gdat <- bind_rows(G_A, G_C) %>%
  distinct(A, C, .keep_all = TRUE) %>%
  group_by(C) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  slice_head(n = TOP_PER_C) %>%
  ungroup()

if (!is.na(TOP_EDGES)) {
  Gdat <- Gdat %>% arrange(desc(score)) %>% slice_head(n = TOP_EDGES)
}

E <- Gdat %>%
  transmute(
    from    = as.character(A),
    to      = as.character(C),
    balance = as.numeric(balance),
    weight  = pmax(as.numeric(evidence), 1e-3)
  )

# ---- sanitize edges/nodes (avoid unnamed nodes) ----
E <- E %>%
  filter(!is.na(from), from != "", !is.na(to), to != "")


# ALWAYS include outcome vertices (even if sparse)
nodes <- tibble(name = unique(c(E$from, E$to, OUTCOMES))) %>%
  mutate(
    type  = if_else(name %in% OUTCOMES, "outcome", "biomarker"),
    label = if_else(is.na(name) | name == "", "(unnamed)", name)
  )


g <- igraph::graph_from_data_frame(E, directed = FALSE, vertices = nodes)


# ---- Always have vertex labels (safer than relying on nodes tibble)
igraph::V(g)$label <- igraph::V(g)$name

# ---- KK layout
set.seed(KK_SEED)
dist_w <- 1 / igraph::E(g)$weight
coords <- igraph::layout_with_kk(g, weights = dist_w)

# ---- Fix Outcome node positions (anchor near center)
vnames <- igraph::V(g)$name
coords2 <- coords

out_idx <- which(vnames %in% OUTCOMES)
if (length(out_idx) > 0) {
  theta <- seq(0, 2*pi, length.out = length(out_idx) + 1)[-(length(out_idx) + 1)]
  coords2[out_idx, 1] <- cos(theta) * 0.20
  coords2[out_idx, 2] <- sin(theta) * 0.20
}

# ---- Place isolated nodes (degree==0) on outer ring so they remain visible
deg <- igraph::degree(g)
iso_idx <- which(deg == 0 & !(vnames %in% OUTCOMES))
if (length(iso_idx) > 0) {
  theta <- seq(0, 2*pi, length.out = length(iso_idx) + 1)[-(length(iso_idx) + 1)]
  coords2[iso_idx, 1] <- cos(theta) * 1.25
  coords2[iso_idx, 2] <- sin(theta) * 1.25
}

pB <- ggraph(g, layout = "manual", x = coords2[,1], y = coords2[,2]) +
  geom_edge_link(aes(edge_width = weight, edge_colour = balance),
                 alpha = 0.85, lineend = "round",
                 end_cap = ggraph::circle(3, "mm")) +
  geom_node_point(aes(shape = type), size = 3.0, colour="grey15") +
  # Keep ggrepel; outcomes are on fixed coordinates so labels won't drift much.
  geom_node_text(aes(label = label), size = 3.0, repel = TRUE) +
  ggraph::scale_edge_colour_gradient2(
    low="steelblue", mid="grey90", high="firebrick",
    midpoint=0, limits=c(-1,1), oob=scales::squish
  ) +
  ggraph::scale_edge_width(range = c(0.22, 2.6)) +
  scale_shape_manual(values = c(outcome = 24, biomarker = 16)) +
  labs(
    title = "Biomarker/family–outcome network (RA-ILD)",
    edge_colour = "Directional balance",
    edge_width  = "Evidence"
  ) +
  theme_void(base_size=10) +
  theme(
    legend.position = "right",
    plot.title = element_text(face="bold", size=12)
  )

# -------------------------
# Panel C
# -------------------------
if (!(FOCUS_OUTCOME %in% OUTCOMES)) FOCUS_OUTCOME <- OUTCOMES[1]

Cdat <- S0 %>%
  filter(C == FOCUS_OUTCOME, A_parent %in% keepA_show) %>%
  mutate(
    lab = if_else(
      rank(-(abs(balance) * (evidence + 0.2)), ties.method="first") <= MAX_LABELS_C,
      as.character(A_parent),
      ""
    )
  )

pC <- ggplot(Cdat, aes(x = balance, y = evidence)) +
  geom_vline(xintercept = 0, linewidth=0.6, colour="grey50") +
  geom_point(aes(size = articles), colour="grey25", alpha=0.9) +
  ggrepel::geom_text_repel(
    aes(label = lab),
    size = 3,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_size_area(max_size = 8) +
  labs(
    title = paste0("Balance–evidence plot (→ ", FOCUS_OUTCOME, ")"),
    x = "Directional balance (risk-decreasing ← 0 → risk-increasing)",
    y = "Evidence volume (log10[pos+neg])",
    size = "Articles"
  ) +
  theme_minimal(base_size=10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold", size=12),
    legend.position = "right"
  )


# -------------------------
# Assemble and save
# -------------------------
fig3 <- (pA | pB) / pC +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(heights = c(1.05, 1))

FIG_OUT_DIR <- if (exists("DIR_FIG3")) DIR_FIG3 else if (exists("DIR_FIG2")) DIR_FIG2 else if (exists("DIR_FIG")) DIR_FIG else stop("No figure output directory found (DIR_FIG3/DIR_FIG2/DIR_FIG).")

stub <- file.path(FIG_OUT_DIR, sprintf("Figure3_%s__%s_BIOMARKER_ATLAS_COMPAT_denseB_v2", CORPUS_TAG, DIC_TAG))
ggsave(paste0(stub, ".pdf"), fig3, width=12, height=8.2, device=grDevices::cairo_pdf)
ggsave(paste0(stub, ".png"), fig3, width=12, height=8.2, dpi=350)

log_msg("WROTE: ", paste0(stub, ".pdf"))
log_msg("WROTE: ", paste0(stub, ".png"))


atlas_out <- S0 %>%
  transmute(
    A = A_parent,
    C = as.character(C),
    n_inc = pos_articles,
    n_dec = neg_articles,
    articles = articles,
    balance = balance,
    evidence = evidence
  ) %>%
  distinct()

write_csv(
  atlas_out,
  file.path(DIR_TABLE, sprintf("biomarker_outcome_atlas_%s__%s.csv", CORPUS_TAG, DIC_TAG))
)