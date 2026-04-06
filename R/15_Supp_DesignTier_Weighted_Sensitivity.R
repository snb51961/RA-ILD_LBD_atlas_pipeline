
# =========================================================
# 15_Supp_DesignTier_Weighted_Sensitivity_v2.R
# Supplementary Tables S2–S5 for the RA-ILD LBD project
#
# PURPOSE
#   Assess whether key signals remain stable after modest study-design
#   weighting at the PMID level, and write manuscript-ready supplementary
#   tables aligned with the table legends.
#
# WRITES
#   - Internal working tables (raw / analysis-friendly)
#   - Manuscript-ready supplementary tables:
#       S2  Design-tier distribution and weights
#       S3  Weighted versus unweighted signed-effect summaries
#       S4  Design-tier weighted biomarker support
#       S5  Weighted versus unweighted non-review drug–outcome co-mentions
#
# IMPORTANT
#   - Main analyses remain unchanged. This script is for sensitivity only.
#   - This is NOT meta-analysis and does NOT use effect-size weighting.
#   - Signed-effects weighting is applied only at article-label aggregation,
#     not by re-running sentence-level inference.
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("readr","dplyr","stringr","purrr","tibble"))

# -------------------------
# User-editable design weights
# Keep these modest.
# -------------------------
W_TIER1_TRIAL        <- 1.50
W_TIER1_STRUCTURED   <- 1.25
W_TIER2_OBSERVATION  <- 1.15
W_TIER3_OTHER        <- 1.00

STAMP_DATE <- format(Sys.Date(), "%Y%m%d")

# -------------------------
# Helpers
# -------------------------
find_latest_file <- function(dir, regex){
  if (is.na(dir) || !dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

find_first_existing <- function(paths){
  ok <- paths[!is.na(paths) & file.exists(paths)]
  if (length(ok)) ok[1] else NA_character_
}

collapse_na <- function(x) ifelse(is.na(x), "", x)

assign_design_tier <- function(pubtype, title = "", abstract = ""){
  pt <- tolower(collapse_na(pubtype))
  tx <- tolower(paste(collapse_na(title), collapse_na(abstract)))

  if (stringr::str_detect(pt, "randomized controlled trial|controlled clinical trial|clinical trial")) {
    return("tier1_trial")
  }

  if (stringr::str_detect(pt, "multicenter study|multicentre study|validation study|evaluation study")) {
    return("tier1_structured")
  }

  if (stringr::str_detect(pt, "observational study|comparative study|case-control")) {
    return("tier2_observational")
  }

  if (stringr::str_detect(
    tx,
    "\\bprospective\\b|\\bretrospective\\b|\\bcohort\\b|\\bcase-control\\b|\\bcross-sectional\\b"
  )) {
    return("tier2_observational")
  }

  "tier3_other_original"
}

design_weight_from_tier <- function(tier){
  dplyr::case_when(
    tier == "tier1_trial"         ~ W_TIER1_TRIAL,
    tier == "tier1_structured"    ~ W_TIER1_STRUCTURED,
    tier == "tier2_observational" ~ W_TIER2_OBSERVATION,
    TRUE                          ~ W_TIER3_OTHER
  )
}

safe_numeric <- function(x){
  suppressWarnings(as.numeric(x))
}

rank_desc <- function(x){
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  dplyr::dense_rank(dplyr::desc(x))
}

pct <- function(x, denom){
  ifelse(is.na(denom) | denom == 0, NA_real_, 100 * x / denom)
}

round_if_not_na <- function(x, digits = 3){
  if (is.numeric(x)) round(x, digits = digits) else x
}

design_tier_order <- function(x){
  dplyr::case_when(
    x == "tier1_trial"         ~ 1L,
    x == "tier1_structured"    ~ 2L,
    x == "tier2_observational" ~ 3L,
    TRUE                       ~ 4L
  )
}

design_tier_label <- function(x){
  dplyr::case_when(
    x == "tier1_trial"         ~ "Tier 1: clinical trial",
    x == "tier1_structured"    ~ "Tier 1: structured original study",
    x == "tier2_observational" ~ "Tier 2: observational study",
    TRUE                       ~ "Tier 3: other original study"
  )
}

design_tier_description <- function(x){
  dplyr::case_when(
    x == "tier1_trial"         ~ "Randomized or controlled clinical trial",
    x == "tier1_structured"    ~ "Multicenter, validation, or evaluation original study",
    x == "tier2_observational" ~ "Observational, cohort, case-control, or cross-sectional study",
    TRUE                       ~ "Other original article not captured above"
  )
}

direction_label <- function(x){
  dplyr::case_when(
    x == "risk_up"            ~ "Risk up",
    x == "risk_down"          ~ "Risk down",
    x == "no_effect_or_mixed" ~ "No effect / mixed",
    TRUE                      ~ as.character(x)
  )
}

get_first_existing_col <- function(df, candidates, default = NA_real_){
  for (nm in candidates) {
    if (nm %in% names(df)) return(df[[nm]])
  }
  rep(default, nrow(df))
}

write_table_csv <- function(df, path){
  readr::write_csv(df, path, na = "")
  invisible(path)
}

# -------------------------
# Locate inputs
# -------------------------
f_articles <- find_first_existing(c(
  file.path(DIR_PROC, sprintf("articles_main_original_%s.csv", STAMP_DATE)),
  file.path(DIR_PROC, sprintf("articles_%s.csv", STAMP_DATE)),
  find_latest_file(DIR_PROC, "^articles_main_original_[0-9]{8}\\.csv$"),
  find_latest_file(DIR_PROC, "^articles_[0-9]{8}\\.csv$")
))

if (is.na(f_articles)) stop("Could not find articles_main_original_*.csv or articles_*.csv in DIR_PROC.")

f_signed_art <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("signed_effects_article_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^signed_effects_article_.*\\.csv$")
))

f_signed_sum <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^signed_effects_summary_withCI_.*\\.csv$")
))

f_biom <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("biomarker_metrics_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^biomarker_metrics_.*\\.csv$")
))

f_hits <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^hits_matrix_.*\\.csv$"),
  find_latest_file(DIR_PROC, "^hits_matrix_[0-9]{8}\\.csv$")
))

f_drug_full <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_full_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^nonreview_cooc_drug_outcomes_full_.*\\.csv$")
))

message("Using articles file: ", f_articles)
if (!is.na(f_signed_art)) message("Using signed_effects_article: ", f_signed_art)
if (!is.na(f_signed_sum)) message("Using signed_effects_summary: ", f_signed_sum)
if (!is.na(f_biom)) message("Using biomarker_metrics: ", f_biom)
if (!is.na(f_hits)) message("Using hits_matrix: ", f_hits)
if (!is.na(f_drug_full)) message("Using nonreview_cooc_drug_outcomes_full: ", f_drug_full)

# -------------------------
# 1) PMID-level design tiers
# -------------------------
ART <- readr::read_csv(f_articles, show_col_types = FALSE) |>
  dplyr::mutate(
    pmid = as.character(pmid),
    design_tier = purrr::pmap_chr(list(pubtype, title, abstract), assign_design_tier),
    design_weight = design_weight_from_tier(design_tier)
  )

# Internal / working output: PMID-level assignment
f_tiers_by_pmid <- file.path(DIR_TABLE, sprintf("design_tier_assignments_by_pmid_%s.csv", STAMP_DATE))
write_table_csv(
  ART |>
    dplyr::select(pmid, year, pubtype, design_tier, design_weight),
  f_tiers_by_pmid
)

# Supplementary Table S2: manuscript-ready
S2 <- ART |>
  dplyr::count(design_tier, design_weight, name = "original_articles_n") |>
  dplyr::mutate(
    tier_order = design_tier_order(design_tier),
    design_tier_label = design_tier_label(design_tier),
    design_tier_description = design_tier_description(design_tier),
    share_original_articles_pct = pct(original_articles_n, sum(original_articles_n))
  ) |>
  dplyr::arrange(tier_order, dplyr::desc(original_articles_n)) |>
  dplyr::transmute(
    `Design tier` = design_tier_label,
    `Description` = design_tier_description,
    `Article-level weight` = round_if_not_na(design_weight, 2),
    `Original articles, n` = original_articles_n,
    `Share of original articles, %` = round_if_not_na(share_original_articles_pct, 1)
  )

f_s2 <- file.path(DIR_TABLE, sprintf("Supplementary_Table_S2_Design_tier_distribution_and_weights_%s.csv", STAMP_DATE))
write_table_csv(S2, f_s2)

# -------------------------
# 2) Weighted signed-effects summary
#    Supplementary Table S3
# -------------------------
if (!is.na(f_signed_art)) {
  SA <- readr::read_csv(f_signed_art, show_col_types = FALSE) |>
    dplyr::mutate(pmid = as.character(pmid)) |>
    dplyr::left_join(
      ART |>
        dplyr::select(pmid, design_tier, design_weight),
      by = "pmid"
    ) |>
    dplyr::mutate(
      design_tier = ifelse(is.na(design_tier), "tier3_other_original", design_tier),
      design_weight = ifelse(is.na(design_weight), W_TIER3_OTHER, design_weight)
    )

  weighted_signed <- SA |>
    dplyr::group_by(A, C) |>
    dplyr::summarise(
      articles_unweighted = dplyr::n(),
      articles_weighted   = sum(design_weight, na.rm = TRUE),
      pos_unweighted      = sum(label == "risk_up", na.rm = TRUE),
      neg_unweighted      = sum(label == "risk_down", na.rm = TRUE),
      null_unweighted     = sum(label == "no_effect_or_mixed", na.rm = TRUE),
      pos_weighted        = sum(design_weight[label == "risk_up"], na.rm = TRUE),
      neg_weighted        = sum(design_weight[label == "risk_down"], na.rm = TRUE),
      null_weighted       = sum(design_weight[label == "no_effect_or_mixed"], na.rm = TRUE),
      weighted_balance    = (pos_weighted - neg_weighted) / pmax(1e-9, pos_weighted + neg_weighted),
      unweighted_balance  = (pos_unweighted - neg_unweighted) / pmax(1, pos_unweighted + neg_unweighted),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dir_unweighted = dplyr::case_when(
        unweighted_balance > 0 ~ "risk_up",
        unweighted_balance < 0 ~ "risk_down",
        TRUE ~ "no_effect_or_mixed"
      ),
      dir_weighted = dplyr::case_when(
        weighted_balance > 0 ~ "risk_up",
        weighted_balance < 0 ~ "risk_down",
        TRUE ~ "no_effect_or_mixed"
      ),
      direction_stable = dir_unweighted == dir_weighted
    ) |>
    dplyr::arrange(dplyr::desc(articles_weighted), dplyr::desc(abs(weighted_balance)))

  # Internal working output
  f_signed_w <- file.path(DIR_TABLE, sprintf("signed_effects_summary_design_weighted_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(weighted_signed, f_signed_w)

  compare_signed <- weighted_signed |>
    dplyr::mutate(
      unweighted_rank = rank_desc(articles_unweighted),
      weighted_rank   = rank_desc(articles_weighted),
      delta_rank      = weighted_rank - unweighted_rank
    )

  # Optional enrich with existing unweighted summary (e.g., CI columns) if available
  if (!is.na(f_signed_sum)) {
    SU <- readr::read_csv(f_signed_sum, show_col_types = FALSE)
    keep_cols <- intersect(
      c("A", "C", "balance_low", "balance_high", "directional_n"),
      names(SU)
    )
    if (length(keep_cols) > 2) {
      compare_signed <- compare_signed |>
        dplyr::left_join(
          SU |>
            dplyr::select(dplyr::all_of(keep_cols)),
          by = c("A", "C")
        )
    }
  }

  # Internal comparison output
  f_signed_cmp <- file.path(DIR_TABLE, sprintf("signed_effects_weighted_vs_unweighted_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(compare_signed, f_signed_cmp)

  # Supplementary Table S3: manuscript-ready
  S3 <- compare_signed |>
    dplyr::arrange(weighted_rank, unweighted_rank, dplyr::desc(articles_weighted), A, C) |>
    dplyr::transmute(
      `A term` = A,
      `Outcome` = C,
      `Unweighted articles, n` = articles_unweighted,
      `Weighted articles` = round_if_not_na(articles_weighted, 2),
      `Unweighted risk-up articles, n` = pos_unweighted,
      `Unweighted risk-down articles, n` = neg_unweighted,
      `Unweighted no-effect/mixed articles, n` = null_unweighted,
      `Weighted risk-up articles` = round_if_not_na(pos_weighted, 2),
      `Weighted risk-down articles` = round_if_not_na(neg_weighted, 2),
      `Weighted no-effect/mixed articles` = round_if_not_na(null_weighted, 2),
      `Unweighted directional balance` = round_if_not_na(unweighted_balance, 3),
      `Weighted directional balance` = round_if_not_na(weighted_balance, 3),
      `Unweighted direction label` = direction_label(dir_unweighted),
      `Weighted direction label` = direction_label(dir_weighted),
      `Unweighted rank` = unweighted_rank,
      `Weighted rank` = weighted_rank,
      `Rank change (weighted - unweighted)` = delta_rank,
      `Direction label unchanged after weighting` = ifelse(direction_stable, "Yes", "No")
    )

  f_s3 <- file.path(DIR_TABLE, sprintf("Supplementary_Table_S3_Weighted_vs_unweighted_signed_effect_summaries_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(S3, f_s3)

} else {
  message("SKIP signed-effects weighting: signed_effects_article_*.csv not found.")
}

# -------------------------
# 3) Weighted biomarker summary
#    Supplementary Table S4
# -------------------------
if (!is.na(f_biom)) {
  BM <- readr::read_csv(f_biom, show_col_types = FALSE) |>
    dplyr::mutate(pmid = as.character(pmid)) |>
    dplyr::left_join(
      ART |>
        dplyr::select(pmid, design_tier, design_weight),
      by = "pmid"
    ) |>
    dplyr::mutate(
      design_tier = ifelse(is.na(design_tier), "tier3_other_original", design_tier),
      design_weight = ifelse(is.na(design_weight), W_TIER3_OTHER, design_weight),
      auc = safe_numeric(auc),
      sens = safe_numeric(sens),
      spec = safe_numeric(spec)
    )

  BM_pmid <- BM |>
    dplyr::group_by(marker, C, pmid) |>
    dplyr::summarise(
      year = if (all(is.na(year))) NA_real_ else suppressWarnings(min(year, na.rm = TRUE)),
      design_weight = dplyr::first(design_weight),
      auc_best  = if (all(is.na(auc))) NA_real_ else max(auc, na.rm = TRUE),
      sens_best = if (all(is.na(sens))) NA_real_ else max(sens, na.rm = TRUE),
      spec_best = if (all(is.na(spec))) NA_real_ else max(spec, na.rm = TRUE),
      .groups = "drop"
    )

  biom_weighted <- BM_pmid |>
    dplyr::group_by(marker, C) |>
    dplyr::summarise(
      pmids_unweighted = dplyr::n(),
      pmids_weighted   = sum(design_weight, na.rm = TRUE),
      auc_n            = sum(!is.na(auc_best)),
      auc_weighted_n   = sum(design_weight[!is.na(auc_best)], na.rm = TRUE),
      auc_mean         = if (all(is.na(auc_best))) NA_real_ else mean(auc_best, na.rm = TRUE),
      auc_weighted_mean = if (all(is.na(auc_best))) NA_real_
                          else stats::weighted.mean(
                            auc_best[!is.na(auc_best)],
                            w = design_weight[!is.na(auc_best)],
                            na.rm = TRUE
                          ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      unweighted_rank = rank_desc(pmids_unweighted),
      weighted_rank   = rank_desc(pmids_weighted),
      delta_rank      = weighted_rank - unweighted_rank
    ) |>
    dplyr::arrange(weighted_rank, unweighted_rank, marker, C)

  # Internal working output
  f_biom_w <- file.path(DIR_TABLE, sprintf("biomarker_metrics_design_weighted_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(biom_weighted, f_biom_w)

  # Supplementary Table S4: manuscript-ready
  S4 <- biom_weighted |>
    dplyr::transmute(
      `Biomarker` = marker,
      `Outcome` = C,
      `Unweighted PMID support, n` = pmids_unweighted,
      `Weighted PMID support` = round_if_not_na(pmids_weighted, 2),
      `Unweighted rank` = unweighted_rank,
      `Weighted rank` = weighted_rank,
      `Rank change (weighted - unweighted)` = delta_rank,
      `AUC summaries available, n` = auc_n,
      `Weighted AUC support` = round_if_not_na(auc_weighted_n, 2),
      `Mean AUC (unweighted)` = round_if_not_na(auc_mean, 3),
      `Mean AUC (weighted)` = round_if_not_na(auc_weighted_mean, 3)
    )

  f_s4 <- file.path(DIR_TABLE, sprintf("Supplementary_Table_S4_Design_tier_weighted_biomarker_support_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(S4, f_s4)

} else {
  message("SKIP biomarker weighting: biomarker_metrics_*.csv not found.")
}

# -------------------------
# 4) Weighted nonreview drug-outcome co-mentions
#    Supplementary Table S5
# -------------------------
if (!is.na(f_hits) && !is.na(f_drug_full)) {
  H <- readr::read_csv(f_hits, show_col_types = FALSE) |>
    dplyr::mutate(pmid = as.character(pmid)) |>
    dplyr::left_join(
      ART |>
        dplyr::select(pmid, design_tier, design_weight),
      by = "pmid"
    ) |>
    dplyr::mutate(
      design_weight = ifelse(is.na(design_weight), W_TIER3_OTHER, design_weight)
    )

  OUT_FULL <- readr::read_csv(f_drug_full, show_col_types = FALSE)

  get_hit <- function(df, term){
    nm <- paste0("hit__", term)
    if (nm %in% names(df)) as.integer(df[[nm]] %in% c(1L, 1, TRUE)) else integer(nrow(df))
  }

  weighted_drug <- purrr::map_dfr(seq_len(nrow(OUT_FULL)), function(i){
    a <- OUT_FULL$A[i]
    c <- OUT_FULL$C[i]
    vA  <- get_hit(H, a)
    vC  <- get_hit(H, c)
    w   <- H$design_weight
    wnA  <- sum(w * vA, na.rm = TRUE)
    wnC  <- sum(w * vC, na.rm = TRUE)
    wnAC <- sum(w * as.integer(vA == 1L & vC == 1L), na.rm = TRUE)
    wratio <- if (wnA > 0) wnAC / wnA else NA_real_

    tibble::tibble(
      A = a,
      C = c,
      weighted_nA = wnA,
      weighted_nC = wnC,
      weighted_nAC = wnAC,
      weighted_ratio = wratio
    )
  })

  # Internal weighted-only output
  f_drug_w <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_design_weighted_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(weighted_drug, f_drug_w)

  nA_unweighted <- safe_numeric(get_first_existing_col(OUT_FULL, c("nA")))
  nC_unweighted <- safe_numeric(get_first_existing_col(OUT_FULL, c("nC")))
  nAC_unweighted <- safe_numeric(get_first_existing_col(OUT_FULL, c("nAC")))
  ratio_unweighted <- safe_numeric(get_first_existing_col(OUT_FULL, c("ratio")))

  weighted_drug_cmp <- OUT_FULL |>
    dplyr::left_join(weighted_drug, by = c("A", "C")) |>
    dplyr::mutate(
      nA = nA_unweighted,
      nC = nC_unweighted,
      nAC = nAC_unweighted,
      ratio = ratio_unweighted,
      unweighted_rank = rank_desc(ratio),
      weighted_rank   = rank_desc(weighted_ratio),
      delta_rank      = weighted_rank - unweighted_rank
    ) |>
    dplyr::arrange(weighted_rank, unweighted_rank, A, C)

  # Internal comparison output
  f_drug_cmp <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_weighted_vs_unweighted_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(weighted_drug_cmp, f_drug_cmp)

  # Supplementary Table S5: manuscript-ready
  S5 <- weighted_drug_cmp |>
    dplyr::transmute(
      `Drug / exposure term (A)` = A,
      `Outcome` = C,
      `Unweighted nA` = nA,
      `Unweighted nC` = nC,
      `Unweighted nAC` = nAC,
      `Unweighted ratio` = round_if_not_na(ratio, 3),
      `Weighted nA` = round_if_not_na(weighted_nA, 2),
      `Weighted nC` = round_if_not_na(weighted_nC, 2),
      `Weighted nAC` = round_if_not_na(weighted_nAC, 2),
      `Weighted ratio` = round_if_not_na(weighted_ratio, 3),
      `Unweighted rank` = unweighted_rank,
      `Weighted rank` = weighted_rank,
      `Rank change (weighted - unweighted)` = delta_rank
    )

  f_s5 <- file.path(DIR_TABLE, sprintf("Supplementary_Table_S5_Weighted_vs_unweighted_nonreview_drug_outcome_co_mentions_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_table_csv(S5, f_s5)

} else {
  message("SKIP drug-outcome weighting: hits_matrix_*.csv and/or nonreview_cooc_drug_outcomes_full_*.csv not found.")
}

# -------------------------
# 5) Table legends / notes
# -------------------------
table_legends <- c(
  "Supplementary Table S2. Design-tier distribution and weights.",
  "Counts of original articles assigned to each pragmatic design tier used in the weighted sensitivity analysis, together with the corresponding article-level weights.",
  "",
  "Supplementary Table S3. Weighted versus unweighted signed-effect summaries.",
  "Side-by-side comparison of unweighted and design-tier weighted article-level signed-effect summaries for A-outcome pairs, including weighted directional balance, weighted direction label, and rank change.",
  "",
  "Supplementary Table S4. Design-tier weighted biomarker support.",
  "Weighted biomarker-outcome support table showing unweighted and weighted PMID support, weighted ranking, and weighted AUC summaries where available.",
  "",
  "Supplementary Table S5. Weighted versus unweighted non-review drug-outcome co-mention summaries.",
  "Comparison of unweighted and design-tier weighted non-review drug-outcome co-mention summaries, including weighted nA, nC, nAC, weighted ratio, and rank change.",
  "",
  "Notes",
  "1) This is a sensitivity analysis only. Main manuscript results remain unchanged.",
  "2) Article-level design weights are modest and intentionally conservative.",
  "3) Signed-effects weighting is applied at article-label aggregation, not by re-running sentence inference.",
  "4) Biomarker weighting reflects weighted PMID support, not effect-size meta-analysis.",
  "5) Drug-outcome weighting requires hits_matrix because aggregated nA/nC/nAC tables alone cannot be reweighted exactly."
)

f_legends <- file.path(DIR_TABLE, sprintf("Supplementary_Table_S2_to_S5_legends_%s.txt", STAMP_DATE))
writeLines(table_legends, con = f_legends)

message("DONE.")
message("Manuscript-ready supplementary tables written to DIR_TABLE:")
message("  - ", f_s2)
if (exists("f_s3")) message("  - ", f_s3)
if (exists("f_s4")) message("  - ", f_s4)
if (exists("f_s5")) message("  - ", f_s5)
message("Legends written to: ", f_legends)
