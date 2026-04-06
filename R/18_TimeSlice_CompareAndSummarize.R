# =========================================================
# 18_TimeSlice_CompareAndSummarize_v4.1.R
# Temporal holdout comparison / validation-like robustness
# for the RA-ILD atlas project
#
# v4.1 design
#   - Keeps "full atlas first, temporal holdout second"
#   - Focuses on discovery-defined core architecture
#   - A: hotspot retention dot plot (wider / more readable)
#   - B: AE-ILD bridge-intermediate retention with discovery importance encoded as point size
#   - C: signed-effect concordance using a stricter primary directional set
#   - D: biomarker-core persistence with values shown for all cells
#   - S6 includes permutation baseline for hotspot strict retention
# =========================================================

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","tibble","purrr","ggplot2","patchwork","forcats","scales","tidyr","grid"))

.this_file <- function(default_name) {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(getwd(), default_name), winslash = "/", mustWork = FALSE)
}
SCRIPT_FILE <- .this_file("18_TimeSlice_CompareAndSummarize_v4.1.R")
SCRIPT_DIR  <- dirname(SCRIPT_FILE)
TH_ROOT     <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = FALSE)
TH_OUTPUT   <- file.path(TH_ROOT, "output")
TH_LOG      <- file.path(TH_OUTPUT, "log")
dir.create(TH_LOG, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(TH_LOG, paste0("18_compare_v4_1_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  message(msg)
  write(msg, file = log_file, append = TRUE)
}

manifest_path <- file.path(TH_OUTPUT, "temporal_holdout_manifest.csv")
if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path, "\nRun 17_TimeSlice_BuildAndAnalyze_v2.R first.")
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

if (!all(c("slice","npmi_path","abc_path","sign_path","biom_path") %in% names(manifest))) {
  stop("Manifest is missing required columns.")
}

OUTCOME_GROUPS <- c("AE-ILD", "progression", "mortality")
HOTSPOT_TOP_K <- 12L
BRIDGE_TOP_K  <- 12L
SIGN_MIN_ARTICLES <- 3L
SIGN_MIN_DIRECTIONAL_N <- 3L
SIGN_MIN_ABS_BALANCE <- 0.33
CORE_BIOMARKERS <- c("KL6","SP_D","ACPA_RF","MMP7","APRIL","Fibrosis_pathway","Cytokine_inflammation")
FULL_HOTSPOT_CORE_DESIRED <- c("male_sex","smoking","disease_duration","UIP","honeycombing","traction_bronchiectasis")
FULL_BIOMARKER_CORE_DESIRED <- c("KL6","SP_D","ACPA_RF","MMP7","APRIL")
N_PERM <- 2000L
set.seed(42)

read_slice_tbl <- function(which_slice, suffix){
  row <- manifest |> dplyr::filter(slice == which_slice)
  if (nrow(row) != 1) stop("Cannot uniquely identify slice in manifest: ", which_slice)
  path <- switch(
    suffix,
    npmi = row$npmi_path[[1]],
    abc  = row$abc_path[[1]],
    sign = row$sign_path[[1]],
    biom = row$biom_path[[1]],
    stop("Unknown suffix: ", suffix)
  )
  if (!file.exists(path)) stop("Missing slice output: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

full_npmi <- read_slice_tbl("full", "npmi")
disc_npmi <- read_slice_tbl("discovery", "npmi")
hold_npmi <- read_slice_tbl("holdout", "npmi")

full_abc <- read_slice_tbl("full", "abc")
disc_abc <- read_slice_tbl("discovery", "abc")
hold_abc <- read_slice_tbl("holdout", "abc")

full_sign <- read_slice_tbl("full", "sign")
disc_sign <- read_slice_tbl("discovery", "sign")
hold_sign <- read_slice_tbl("holdout", "sign")

full_biom <- read_slice_tbl("full", "biom")
disc_biom <- read_slice_tbl("discovery", "biom")
hold_biom <- read_slice_tbl("holdout", "biom")

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
ensure_balance <- function(D) {
  if (!nrow(D)) return(D)
  if (!("balance" %in% names(D))) {
    if (all(c("pos_articles","neg_articles") %in% names(D))) {
      D <- D |> dplyr::mutate(balance = (pos_articles - neg_articles) / pmax(1, pos_articles + neg_articles))
    } else {
      D <- D |> dplyr::mutate(balance = NA_real_)
    }
  }
  D
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

wilson_prop_ci <- function(x, n){
  if (is.na(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  p <- x / n
  z <- 1.96
  denom <- 1 + z^2 / n
  center <- (p + z^2/(2*n)) / denom
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / denom
  c(low = center - half, high = center + half)
}

blank_plot <- function(title, subtitle = NULL){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No eligible data", size = 4) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 11)
}

topness_from_rank <- function(rank, n_total){
  ifelse(is.na(rank) | is.na(n_total) | n_total <= 0, NA_real_, 1 - ((rank - 1) / pmax(1, n_total - 1)))
}

sign_label <- function(x) ifelse(x > 0, "risk_up", ifelse(x < 0, "risk_down", "no_effect_or_mixed"))

# ---------------------------------------------------------
# 0) Full-corpus core check
# ---------------------------------------------------------
full_hot_core_available <- intersect(FULL_HOTSPOT_CORE_DESIRED, unique(full_npmi$A))
full_hot_core_check <- full_npmi |>
  dplyr::filter(A %in% full_hot_core_available) |>
  dplyr::group_by(A) |>
  dplyr::summarise(
    any_loose = any(npmi > 0 & lift > 1, na.rm = TRUE),
    best_npmi = max(npmi, na.rm = TRUE),
    best_lift = max(lift, na.rm = TRUE),
    .groups = "drop"
  )

full_biom_core_available <- intersect(FULL_BIOMARKER_CORE_DESIRED, unique(full_biom$A))
full_biom_core_check <- full_biom |>
  dplyr::filter(A %in% full_biom_core_available) |>
  dplyr::group_by(A) |>
  dplyr::summarise(
    any_present = dplyr::n() > 0,
    best_score = max(score, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(full_hot_core_check, file.path(TH_OUTPUT, "temporal_full_hotspot_core_check_v4_1.csv"))
readr::write_csv(full_biom_core_check, file.path(TH_OUTPUT, "temporal_full_biomarker_core_check_v4_1.csv"))

# ---------------------------------------------------------
# 1) Hotspot validation
# ---------------------------------------------------------
rank_hotspots <- function(D){
  if (!nrow(D)) return(tibble::tibble(C=character(), A=character(), npmi=double(), lift=double(), score_hot=double(), rank_hot=integer(), n_total=integer(), topness=double()))
  D |>
    dplyr::filter(C %in% OUTCOME_GROUPS) |>
    dplyr::mutate(log2lift = log2(pmax(lift, 1e-12))) |>
    dplyr::group_by(C) |>
    dplyr::mutate(
      rank_npmi = dplyr::dense_rank(dplyr::desc(npmi)),
      rank_lift = dplyr::dense_rank(dplyr::desc(log2lift)),
      score_hot = -(rank_npmi + rank_lift),
      rank_hot = dplyr::dense_rank(dplyr::desc(score_hot)),
      n_total = dplyr::n(),
      topness = topness_from_rank(rank_hot, n_total)
    ) |>
    dplyr::ungroup()
}

disc_hot_ranked <- rank_hotspots(disc_npmi)
hold_hot_ranked <- rank_hotspots(hold_npmi)

hotspot_primary <- disc_hot_ranked |>
  dplyr::group_by(C) |>
  dplyr::arrange(rank_hot, .by_group = TRUE) |>
  dplyr::slice_head(n = HOTSPOT_TOP_K) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    C, A,
    disc_rank = rank_hot,
    disc_topness = topness,
    disc_npmi = npmi,
    disc_lift = lift
  )

hotspot_validation <- hotspot_primary |>
  dplyr::left_join(
    hold_hot_ranked |>
      dplyr::select(C, A, hold_rank = rank_hot, hold_topness = topness, hold_npmi = npmi, hold_lift = lift),
    by = c("C","A")
  ) |>
  dplyr::mutate(
    strict_retained = !is.na(hold_npmi) & hold_npmi >= 0.05 & hold_lift >= 1.2,
    loose_retained  = !is.na(hold_npmi) & hold_npmi > 0 & hold_lift > 1,
    retention_class = dplyr::case_when(
      strict_retained ~ "Strict",
      loose_retained  ~ "Loose",
      TRUE            ~ "Not retained"
    )
  )

hotspot_summary <- hotspot_validation |>
  dplyr::group_by(C) |>
  dplyr::summarise(
    n_primary = dplyr::n(),
    strict_retained_n = sum(strict_retained, na.rm = TRUE),
    loose_retained_n  = sum(loose_retained, na.rm = TRUE),
    median_hold_topness = stats::median(hold_topness, na.rm = TRUE),
    spearman_rho = safe_spearman(disc_rank, hold_rank),
    .groups = "drop"
  )

# permutation baseline
perm_baseline_one <- function(cg){
  disc_terms <- hotspot_primary |> dplyr::filter(C == cg) |> dplyr::pull(A)
  hold_pool  <- hold_hot_ranked |> dplyr::filter(C == cg) |> dplyr::pull(A)

  if (!length(disc_terms) || !length(hold_pool)) {
    return(tibble::tibble(
      C = cg,
      observed_strict = NA_real_,
      expected_strict = NA_real_,
      p_perm = NA_real_
    ))
  }

  strict_hold_terms <- hold_hot_ranked |>
    dplyr::filter(C == cg, npmi >= 0.05, lift >= 1.2) |>
    dplyr::pull(A)

  observed <- sum(disc_terms %in% strict_hold_terms)
  k <- length(disc_terms)

  sims <- replicate(N_PERM, {
    sample_terms <- sample(hold_pool, size = min(k, length(hold_pool)), replace = FALSE)
    sum(sample_terms %in% strict_hold_terms)
  })

  tibble::tibble(
    C = cg,
    observed_strict = observed,
    expected_strict = mean(sims),
    p_perm = mean(sims >= observed)
  )
}

hotspot_perm <- dplyr::bind_rows(lapply(OUTCOME_GROUPS, perm_baseline_one))

readr::write_csv(hotspot_primary, file.path(TH_OUTPUT, "temporal_hotspot_primary_targets_v4_1.csv"))
readr::write_csv(hotspot_validation, file.path(TH_OUTPUT, "temporal_hotspot_validation_details_v4_1.csv"))
readr::write_csv(hotspot_summary, file.path(TH_OUTPUT, "temporal_hotspot_validation_summary_v4_1.csv"))
readr::write_csv(hotspot_perm, file.path(TH_OUTPUT, "temporal_hotspot_permutation_baseline_v4_1.csv"))

# ---------------------------------------------------------
# 2) AE-ILD bridge-intermediate validation
# ---------------------------------------------------------
aggregate_B <- function(abc_df){
  if (!nrow(abc_df)) return(tibble::tibble(B=character(), B_support=double(), B_rank=integer(), n_total=integer(), topness=double(), triads=integer()))
  abc_df |>
    dplyr::filter(C == "AE-ILD") |>
    dplyr::group_by(B) |>
    dplyr::summarise(
      B_support = sum(score_q, na.rm = TRUE),
      triads = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(B_support), dplyr::desc(triads)) |>
    dplyr::mutate(
      B_rank = dplyr::row_number(),
      n_total = dplyr::n(),
      topness = topness_from_rank(B_rank, n_total)
    )
}

disc_B <- aggregate_B(disc_abc)
hold_B <- aggregate_B(hold_abc)

bridge_primary <- disc_B |>
  dplyr::slice_head(n = BRIDGE_TOP_K) |>
  dplyr::transmute(
    B,
    disc_rank = B_rank,
    disc_topness = topness,
    disc_support = B_support,
    disc_triads = triads
  )

bridge_validation <- bridge_primary |>
  dplyr::left_join(
    hold_B |>
      dplyr::select(B, hold_rank = B_rank, hold_topness = topness, hold_support = B_support, hold_triads = triads),
    by = "B"
  ) |>
  dplyr::mutate(
    strict_retained = !is.na(hold_topness) & hold_topness >= 0.75,
    loose_retained  = !is.na(hold_topness) & hold_topness >= 0.50,
    retention_class = dplyr::case_when(
      strict_retained ~ "Top 25%",
      loose_retained  ~ "Top 50%",
      TRUE            ~ "Lower"
    )
  )

bridge_summary <- tibble::tibble(
  n_primary = nrow(bridge_validation),
  strict_retained_n = sum(bridge_validation$strict_retained, na.rm = TRUE),
  loose_retained_n  = sum(bridge_validation$loose_retained, na.rm = TRUE),
  median_hold_topness = stats::median(bridge_validation$hold_topness, na.rm = TRUE),
  spearman_rho = safe_spearman(bridge_validation$disc_rank, bridge_validation$hold_rank)
)

readr::write_csv(bridge_primary, file.path(TH_OUTPUT, "temporal_bridge_primary_targets_v4_1.csv"))
readr::write_csv(bridge_validation, file.path(TH_OUTPUT, "temporal_bridge_validation_details_v4_1.csv"))
readr::write_csv(bridge_summary, file.path(TH_OUTPUT, "temporal_bridge_validation_summary_v4_1.csv"))

# ---------------------------------------------------------
# 3) Signed-effect validation (stricter primary set)
# ---------------------------------------------------------
disc_sign <- ensure_balance(disc_sign)
hold_sign <- ensure_balance(hold_sign)

if (all(c("pos_articles","neg_articles") %in% names(disc_sign))) {
  disc_sign$directional_n <- as.integer(disc_sign$pos_articles + disc_sign$neg_articles)
} else {
  disc_sign$directional_n <- as.integer(disc_sign$articles)
}
disc_sign <- disc_sign |> dplyr::mutate(disc_label = sign_label(balance))

if (all(c("pos_articles","neg_articles") %in% names(hold_sign))) {
  hold_sign$directional_n <- as.integer(hold_sign$pos_articles + hold_sign$neg_articles)
} else {
  hold_sign$directional_n <- as.integer(hold_sign$articles)
}
hold_sign <- hold_sign |> dplyr::mutate(hold_label = sign_label(balance))

sign_primary_all <- disc_sign |>
  dplyr::filter(
    articles >= SIGN_MIN_ARTICLES,
    directional_n >= SIGN_MIN_DIRECTIONAL_N,
    abs(balance) >= SIGN_MIN_ABS_BALANCE
  ) |>
  dplyr::select(A, C, articles, directional_n, balance, disc_label) |>
  dplyr::distinct()

sign_primary_nonzero <- disc_sign |>
  dplyr::filter(
    articles >= SIGN_MIN_ARTICLES,
    directional_n >= SIGN_MIN_DIRECTIONAL_N,
    abs(balance) >= SIGN_MIN_ABS_BALANCE,
    disc_label != "no_effect_or_mixed"
  ) |>
  dplyr::select(A, C, articles, directional_n, balance, disc_label) |>
  dplyr::distinct()

eval_signed_set <- function(primary_df, set_name){
  if (!nrow(primary_df)) {
    return(tibble::tibble(
      set = set_name,
      n_primary = 0L,
      n_evaluable = 0L,
      evaluable_prop = NA_real_,
      concordance = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  cmp <- primary_df |>
    dplyr::left_join(
      hold_sign |>
        dplyr::select(A, C, hold_label, hold_balance = balance, hold_articles = articles),
      by = c("A","C")
    ) |>
    dplyr::mutate(
      evaluable = !is.na(hold_label),
      same_direction = evaluable & disc_label == hold_label
    )

  n_primary <- nrow(cmp)
  n_evaluable <- sum(cmp$evaluable, na.rm = TRUE)
  concordance <- if (n_evaluable > 0) mean(cmp$same_direction[cmp$evaluable]) else NA_real_
  ci <- wilson_prop_ci(sum(cmp$same_direction, na.rm = TRUE), n_evaluable)

  readr::write_csv(cmp, file.path(TH_OUTPUT, paste0("temporal_signed_validation_", gsub("[^A-Za-z0-9]+", "_", set_name), "_details_v4_1.csv")))

  tibble::tibble(
    set = set_name,
    n_primary = n_primary,
    n_evaluable = n_evaluable,
    evaluable_prop = n_evaluable / pmax(1, n_primary),
    concordance = concordance,
    ci_low = ci["low"],
    ci_high = ci["high"]
  )
}

signed_summary <- dplyr::bind_rows(
  eval_signed_set(sign_primary_all, "All primary pairs"),
  eval_signed_set(sign_primary_nonzero, "Non-zero primary pairs")
)
readr::write_csv(signed_summary, file.path(TH_OUTPUT, "temporal_signed_validation_summary_v4_1.csv"))

# ---------------------------------------------------------
# 4) Biomarker-core persistence
# ---------------------------------------------------------
rank_biom <- function(D){
  if (!nrow(D)) return(tibble::tibble(A=character(), C=character(), score=double(), rank_b=integer(), n_total=integer(), topness=double()))
  D |>
    dplyr::group_by(C) |>
    dplyr::arrange(dplyr::desc(score), .by_group = TRUE) |>
    dplyr::mutate(
      rank_b = dplyr::row_number(),
      n_total = dplyr::n(),
      topness = topness_from_rank(rank_b, n_total)
    ) |>
    dplyr::ungroup()
}

disc_biom_ranked <- rank_biom(disc_biom)
hold_biom_ranked <- rank_biom(hold_biom)

core_grid <- tidyr::expand_grid(
  A = CORE_BIOMARKERS,
  C = OUTCOME_GROUPS
)

biom_core_validation <- core_grid |>
  dplyr::left_join(
    disc_biom_ranked |>
      dplyr::select(A, C, disc_rank = rank_b, disc_topness = topness, disc_score = score),
    by = c("A","C")
  ) |>
  dplyr::left_join(
    hold_biom_ranked |>
      dplyr::select(A, C, hold_rank = rank_b, hold_topness = topness, hold_score = score),
    by = c("A","C")
  ) |>
  dplyr::mutate(
    discovery_present = !is.na(disc_rank),
    strict_retained = discovery_present & !is.na(hold_topness) & hold_topness >= 0.67,
    loose_retained  = discovery_present & !is.na(hold_topness) & hold_topness >= 0.50,
    retention_class = dplyr::case_when(
      !discovery_present ~ "Not prioritised in discovery",
      strict_retained    ~ "Strict",
      loose_retained     ~ "Loose",
      TRUE               ~ "Not retained"
    ),
    hold_topness_label = ifelse(is.na(hold_topness), "", sprintf("%.2f", hold_topness))
  )

biom_core_summary <- biom_core_validation |>
  dplyr::filter(discovery_present) |>
  dplyr::group_by(C) |>
  dplyr::summarise(
    n_discovery_core_pairs = dplyr::n(),
    strict_retained_n = sum(strict_retained, na.rm = TRUE),
    loose_retained_n  = sum(loose_retained, na.rm = TRUE),
    median_hold_topness = stats::median(hold_topness, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(biom_core_validation, file.path(TH_OUTPUT, "temporal_biomarker_core_validation_details_v4_1.csv"))
readr::write_csv(biom_core_summary, file.path(TH_OUTPUT, "temporal_biomarker_core_validation_summary_v4_1.csv"))

# ---------------------------------------------------------
# 5) Summary table S6 v4.1
# ---------------------------------------------------------
S6_full <- dplyr::bind_rows(
  tibble::tibble(
    section = "Full-corpus pipeline core check",
    item = "Hotspot core terms captured (loose threshold)",
    metric1 = sum(full_hot_core_check$any_loose, na.rm = TRUE),
    metric2 = nrow(full_hot_core_check),
    metric3 = NA_real_,
    metric4 = NA_real_,
    expected_strict = NA_real_,
    p_perm = NA_real_
  ),
  tibble::tibble(
    section = "Full-corpus pipeline core check",
    item = "Biomarker core terms captured",
    metric1 = sum(full_biom_core_check$any_present, na.rm = TRUE),
    metric2 = nrow(full_biom_core_check),
    metric3 = NA_real_,
    metric4 = NA_real_,
    expected_strict = NA_real_,
    p_perm = NA_real_
  )
)

S6_hot <- hotspot_summary |>
  dplyr::left_join(hotspot_perm, by = "C") |>
  dplyr::transmute(
    section = "Hotspot validation",
    item = C,
    metric1 = strict_retained_n,
    metric2 = loose_retained_n,
    metric3 = n_primary,
    metric4 = spearman_rho,
    expected_strict = expected_strict,
    p_perm = p_perm
  )

S6_bridge <- tibble::tibble(
  section = "AE-ILD bridge validation",
  item = "Bridge intermediates",
  metric1 = bridge_summary$strict_retained_n,
  metric2 = bridge_summary$loose_retained_n,
  metric3 = bridge_summary$n_primary,
  metric4 = bridge_summary$spearman_rho,
  expected_strict = NA_real_,
  p_perm = NA_real_
)

S6_signed <- signed_summary |>
  dplyr::transmute(
    section = "Signed-effect validation",
    item = set,
    metric1 = concordance,
    metric2 = evaluable_prop,
    metric3 = n_evaluable,
    metric4 = n_primary,
    expected_strict = NA_real_,
    p_perm = NA_real_
  )

S6_biom <- biom_core_summary |>
  dplyr::transmute(
    section = "Biomarker-core validation",
    item = C,
    metric1 = strict_retained_n,
    metric2 = loose_retained_n,
    metric3 = n_discovery_core_pairs,
    metric4 = median_hold_topness,
    expected_strict = NA_real_,
    p_perm = NA_real_
  )

S6 <- dplyr::bind_rows(S6_full, S6_hot, S6_bridge, S6_signed, S6_biom)
S6_path <- file.path(TH_OUTPUT, "Supplementary_Table_S6_temporal_holdout_summary_v4_1.csv")
readr::write_csv(S6, S6_path)

# ---------------------------------------------------------
# 6) Figure S5 v4.1
# ---------------------------------------------------------
# A: wider readable v2-like dot plot
pA_dat <- hotspot_validation |>
  dplyr::mutate(
    term_label = paste0(A, " [D", disc_rank, "]"),
    term_label = forcats::fct_reorder(term_label, disc_rank, .desc = TRUE),
    retention_class = factor(retention_class, levels = c("Strict","Loose","Not retained"))
  )

pA <- if (nrow(pA_dat)) {
  ggplot2::ggplot(pA_dat, ggplot2::aes(x = hold_topness, y = term_label, colour = retention_class)) +
    ggplot2::geom_vline(xintercept = c(0.50, 0.67), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(size = 2.6, alpha = 0.95) +
    ggplot2::facet_wrap(~ C, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    ggplot2::labs(
      title = "A. Discovery-defined hotspot retention",
      subtitle = "Point colour: retention class; x-position: holdout topness",
      x = "Holdout topness (1 = highest within outcome)",
      y = NULL,
      colour = "Retention"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.spacing.x = grid::unit(1.4, "lines"),
      plot.margin = ggplot2::margin(5.5, 18, 5.5, 5.5)
    )
} else {
  blank_plot("A. Discovery-defined hotspot retention")
}

# B: discovery importance retained
pB_dat <- bridge_validation |>
  dplyr::mutate(
    B = forcats::fct_reorder(B, disc_rank, .desc = TRUE),
    retention_class = factor(retention_class, levels = c("Top 25%","Top 50%","Lower"))
  )

pB <- if (nrow(pB_dat)) {
  ggplot2::ggplot(pB_dat, ggplot2::aes(x = hold_topness, y = B, colour = retention_class, size = disc_support)) +
    ggplot2::geom_vline(xintercept = c(0.50, 0.75), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_size_continuous(name = "Discovery\nimportance") +
    ggplot2::labs(
      title = "B. AE-ILD bridge-intermediate retention",
      x = "Holdout topness (1 = highest among B nodes)",
      y = NULL,
      colour = "Retention"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
} else {
  blank_plot("B. AE-ILD bridge-intermediate retention")
}

# C: stricter primary directional set
pC_dat <- signed_summary |>
  dplyr::mutate(
    set = factor(set, levels = c("All primary pairs","Non-zero primary pairs")),
    label = paste0(n_evaluable, "/", n_primary, " evaluable")
  )

pC <- if (nrow(pC_dat)) {
  ggplot2::ggplot(pC_dat, ggplot2::aes(x = set, y = concordance)) +
    ggplot2::geom_col(fill = "#59A14F") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.3, size = 3.4) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "C. Signed-effect concordance",
      x = NULL,
      y = "Concordance"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
} else {
  blank_plot("C. Signed-effect concordance")
}

# D: numeric labels for all cells
pD_dat <- biom_core_validation |>
  dplyr::mutate(
    A = factor(A, levels = rev(CORE_BIOMARKERS)),
    C = factor(C, levels = OUTCOME_GROUPS),
    retention_class = factor(retention_class, levels = c("Strict","Loose","Not retained","Not prioritised in discovery"))
  )

pD <- if (nrow(pD_dat)) {
  ggplot2::ggplot(pD_dat, ggplot2::aes(x = C, y = A, fill = retention_class)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = hold_topness_label), size = 3) +
    ggplot2::scale_fill_manual(values = c(
      "Strict" = "#F8766D",
      "Loose" = "#7CAE00",
      "Not retained" = "#00BFC4",
      "Not prioritised in discovery" = "#C77CFF"
    )) +
    ggplot2::labs(
      title = "D. Biomarker-core persistence",
      x = NULL,
      y = NULL,
      fill = "Retention"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
} else {
  blank_plot("D. Biomarker-core persistence")
}

figS5 <- ((pA | pB) + patchwork::plot_layout(widths = c(1.35, 1))) /
         ((pC | pD) + patchwork::plot_layout(widths = c(1, 1))) +
  patchwork::plot_annotation(
    title = "Supplementary Figure S5. Temporal holdout robustness of the RA-ILD literature atlas",
    subtitle = "Discovery: 1980-2020; Holdout: 2021-2025. Validation targets were fixed in discovery and evaluated unchanged in holdout."
  )

S5_stub <- file.path(TH_OUTPUT, "SupplementaryFigureS5_temporal_holdout_revised_v4_1")
ggplot2::ggsave(paste0(S5_stub, ".pdf"), figS5, width = 16, height = 10, device = grDevices::cairo_pdf)
ggplot2::ggsave(paste0(S5_stub, ".png"), figS5, width = 16, height = 10, dpi = 350)

readme_path <- file.path(TH_OUTPUT, "README_compare_v4_1.txt")
writeLines(c(
  "Temporal holdout comparison completed (v4.1).",
  paste0("Manifest: ", manifest_path),
  paste0("Summary table: ", S6_path),
  paste0("Figure S5 v4.1: ", paste0(S5_stub, ".pdf")),
  "",
  "Primary outputs:",
  "  temporal_hotspot_validation_details_v4_1.csv",
  "  temporal_bridge_validation_details_v4_1.csv",
  "  temporal_signed_validation_summary_v4_1.csv",
  "  temporal_biomarker_core_validation_details_v4_1.csv",
  "  temporal_hotspot_permutation_baseline_v4_1.csv"
), con = readme_path)

log_msg("WROTE S6 v4.1:", S6_path)
log_msg("WROTE S5 v4.1:", paste0(S5_stub, ".pdf"))
log_msg("=== DONE 18_TimeSlice_CompareAndSummarize_v4.1 ===")