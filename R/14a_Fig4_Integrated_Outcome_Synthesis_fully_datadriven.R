# =========================================================
# 14_Fig4_Integrated_Outcome_Synthesis_fully_datadriven_v5.R
# Figure 4 = Outcome-specific hot spot map in RA-ILD
#
# Revision:
#   - removes the manual 'Dominant theme' layer
#   - makes BOTH displayed layers data-supported
#       * Hot spots: selected from signed-effect outputs (prefer weighted summaries)
#       * Biomarker / exploratory layer: selected from biomarker support outputs
#   - keeps the left-side row labels outside the content boxes
#   - expands both displayed layers by ~5 additional terms per outcome
#   - widens the left label gutter and increases vertical spacing for readability
#   - writes provenance tables for both layers
# =========================================================

.find_setup <- function(fname = "00_setup_Final.R") {
  p1 <- file.path(getwd(), fname)
  if (file.exists(p1)) return(p1)

  cur <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:5) {
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    p <- file.path(parent, fname)
    if (file.exists(p)) return(p)
    cur <- parent
  }

  root_env <- Sys.getenv("RAILD_ROOT", unset = NA_character_)
  if (!is.na(root_env) && nzchar(root_env)) {
    p3 <- file.path(root_env, fname)
    if (file.exists(p3)) return(p3)
  }

  stop("Cannot find ", fname, ". Run from the project root or set RAILD_ROOT.")
}

source(.find_setup())

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("readr","dplyr","stringr","tibble","ggplot2","patchwork","grid","scales","purrr"))
library(grid)

OUTCOME_LEVELS <- c("AE-ILD", "Progression", "Mortality")
HOTSPOT_TERMS_PER_OUTCOME <- 18
LOWER_TERMS_PER_OUTCOME <- 18
HOTSPOT_TERMS_PER_LINE <- 3
LOWER_TERMS_PER_LINE <- 3
STRICT_DATA_DRIVEN <- TRUE

# Hot-spot layer should reflect non-biomarker, non-drug manuscript-facing terms.
HOTSPOT_EXCLUDE_CLASSES <- c(
  "biomarker", "molecular", "cell", "microbiome",
  "drug", "bio", "gene", "vaccine"
)
HOTSPOT_DROP_GENERIC <- c(
  "CYTOKINE_INFLAMMATION", "FIBROSIS_PATHWAY", "OXIDATIVE_STRESS_PATHWAY",
  "IMMUNE_CHECKPOINT", "DISEASE_ACTIVITY", "SEVERITY", "BASELINE_SEVERITY"
)

# Terms that are useful in broader atlas views but too generic/redundant for the compact lower layer.
LOWER_DROP_GENERIC <- c(
  "CYTOKINE_INFLAMMATION", "FIBROSIS_PATHWAY", "OXIDATIVE_STRESS_PATHWAY", "IMMUNE_CHECKPOINT"
)
LOWER_DROP_PREFER <- c(
  "ACPA_RF", "CRP_ESR", "NLR", "CAR", "MPO"
)

# Outcome card colours (visual only; no inferential meaning)
CARD_FILL <- tibble::tribble(
  ~outcome,       ~fill,
  "AE-ILD",      "#C74B50",
  "Progression", "#E08E2B",
  "Mortality",   "#4E79A7"
)

find_latest_file <- function(dir, regex){
  if (!dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

find_first_existing <- function(paths){
  ok <- paths[file.exists(paths)]
  if (length(ok)) ok[1] else NA_character_
}

first_existing_col <- function(df, choices){
  hit <- choices[choices %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}

num_from_choices <- function(df, choices){
  col <- first_existing_col(df, choices)
  if (is.na(col)) rep(NA_real_, nrow(df)) else suppressWarnings(as.numeric(df[[col]]))
}

chr_from_choices <- function(df, choices){
  col <- first_existing_col(df, choices)
  if (is.na(col)) rep(NA_character_, nrow(df)) else as.character(df[[col]])
}

canon_outcome <- function(x){
  xl <- stringr::str_to_lower(stringr::str_replace_all(as.character(x), "–", "-"))
  dplyr::case_when(
    stringr::str_detect(xl, "^ae-ild$|^aeild$|acute exacerb|exacerb") ~ "AE-ILD",
    stringr::str_detect(xl, "^progression$|progress|decline|fvc|dlco|worsen") ~ "Progression",
    stringr::str_detect(xl, "^mortality$|mortality|death|surviv") ~ "Mortality",
    TRUE ~ NA_character_
  )
}

norm_term <- function(x){
  x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_replace_all("–", "-") |>
    stringr::str_replace_all("\\s+", "_") |>
    toupper()
}

format_term_display <- function(x){
  x1 <- as.character(x)
  nx <- norm_term(x1)

  out <- dplyr::case_when(
    nx == "KL6" ~ "KL-6",
    nx %in% c("SP_D", "SPD") ~ "SP-D",
    nx == "MMP7" ~ "MMP-7",
    nx == "YKL40" ~ "YKL-40",
    nx == "ACPA_RF" ~ "ACPA/RF",
    nx == "IL17A" ~ "IL-17A",
    nx == "TNFA" ~ "TNF-α",
    nx == "TGFB" ~ "TGF-β",
    nx == "IFNG" ~ "IFN-γ",
    nx == "CXCL13" ~ "CXCL13",
    nx == "CCL18" ~ "CCL18",
    nx == "TH17" ~ "Th17",
    nx == "TREG" ~ "Treg",
    nx == "TH1" ~ "Th1",
    nx == "TH2" ~ "Th2",
    nx == "UIP" ~ "UIP",
    nx == "CPFE" ~ "CPFE",
    nx == "FVC_DECLINE" ~ "FVC decline",
    nx == "DLCO_DECLINE" ~ "DLCO decline",
    nx == "TRACTION_BRONCHIECTASIS" ~ "Traction bronchiectasis",
    nx == "FIBROSIS_EXTENT" ~ "Fibrosis extent",
    nx == "HONEYCOMBING" ~ "Honeycombing",
    nx == "MALE_SEX" ~ "Male sex",
    nx == "SMOKING" ~ "Smoking",
    nx == "CPI" ~ "CPI",
    nx == "GAP" ~ "GAP",
    TRUE ~ stringr::str_replace_all(x1, "_", " ")
  )

  out
}

collapse_terms <- function(x, per_line = 2){
  x <- unique(stats::na.omit(as.character(x)))
  if (!length(x)) return("NA")
  idx <- ceiling(seq_along(x) / per_line)
  lines <- split(x, idx)
  paste(vapply(lines, function(z) paste(z, collapse = " / "), character(1)), collapse = "\n")
}

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

build_dic_lookup <- function(){
  dic_path <- pick_dic()
  if (is.na(dic_path) || !file.exists(dic_path)) return(NULL)

  DIC <- readr::read_csv(dic_path, show_col_types = FALSE)
  if (!("term" %in% names(DIC))) return(NULL)

  out <- tibble::tibble(
    term_raw = as.character(DIC$term),
    term_norm = norm_term(DIC$term),
    class_low = if ("class" %in% names(DIC)) tolower(as.character(DIC$class)) else NA_character_,
    role_up = if ("role" %in% names(DIC)) toupper(as.character(DIC$role)) else NA_character_,
    dic_path = dic_path
  ) |>
    dplyr::filter(!is.na(term_norm), term_norm != "") |>
    dplyr::group_by(term_norm) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  out
}

prep_hotspot_candidates <- function(dic_lookup){
  f_cmp <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("signed_effects_weighted_vs_unweighted_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^signed_effects_weighted_vs_unweighted_.*\\.csv$")
  ))

  f_s3 <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("Supplementary_Table_S3_Weighted_vs_unweighted_signed_effect_summaries_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^Supplementary_Table_S3_Weighted_vs_unweighted_signed_effect_summaries_.*\\.csv$")
  ))

  f_weighted <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("signed_effects_summary_design_weighted_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^signed_effects_summary_design_weighted_.*\\.csv$")
  ))

  f_unweighted <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^signed_effects_summary_withCI_.*\\.csv$")
  ))

  out <- list()

  if (!is.na(f_cmp)) {
    X <- readr::read_csv(f_cmp, show_col_types = FALSE)
    a_col <- first_existing_col(X, c("A", "A term"))
    c_col <- first_existing_col(X, c("C", "Outcome"))
    if (!is.na(a_col) && !is.na(c_col)) {
      out[[length(out) + 1]] <- tibble::tibble(
        outcome = canon_outcome(X[[c_col]]),
        term_raw = as.character(X[[a_col]]),
        term_norm = norm_term(X[[a_col]]),
        display = format_term_display(X[[a_col]]),
        support = dplyr::coalesce(
          num_from_choices(X, c("articles_weighted", "Weighted articles")),
          num_from_choices(X, c("articles_unweighted", "Unweighted articles, n"))
        ),
        rank_value = dplyr::coalesce(
          num_from_choices(X, c("weighted_rank", "Weighted rank")),
          num_from_choices(X, c("unweighted_rank", "Unweighted rank"))
        ),
        aux_value = dplyr::coalesce(
          abs(num_from_choices(X, c("weighted_balance", "Weighted directional balance"))),
          abs(num_from_choices(X, c("unweighted_balance", "Unweighted directional balance")))
        ),
        source_type = "Signed-effects weighted vs unweighted comparison",
        source_file = f_cmp,
        source_priority = 1,
        metric_used = "weighted articles + |weighted balance|"
      )
    }
  }

  if (!is.na(f_s3) && is.na(f_cmp)) {
    X <- readr::read_csv(f_s3, show_col_types = FALSE)
    a_col <- first_existing_col(X, c("A", "A term"))
    c_col <- first_existing_col(X, c("C", "Outcome"))
    if (!is.na(a_col) && !is.na(c_col)) {
      out[[length(out) + 1]] <- tibble::tibble(
        outcome = canon_outcome(X[[c_col]]),
        term_raw = as.character(X[[a_col]]),
        term_norm = norm_term(X[[a_col]]),
        display = format_term_display(X[[a_col]]),
        support = dplyr::coalesce(
          num_from_choices(X, c("Weighted articles")),
          num_from_choices(X, c("Unweighted articles, n"))
        ),
        rank_value = dplyr::coalesce(
          num_from_choices(X, c("Weighted rank")),
          num_from_choices(X, c("Unweighted rank"))
        ),
        aux_value = dplyr::coalesce(
          abs(num_from_choices(X, c("Weighted directional balance"))),
          abs(num_from_choices(X, c("Unweighted directional balance")))
        ),
        source_type = "Supplementary Table S3",
        source_file = f_s3,
        source_priority = 2,
        metric_used = "Weighted articles + |weighted balance|"
      )
    }
  }

  if (!is.na(f_weighted)) {
    X <- readr::read_csv(f_weighted, show_col_types = FALSE)
    a_col <- first_existing_col(X, c("A", "A term"))
    c_col <- first_existing_col(X, c("C", "Outcome"))
    if (!is.na(a_col) && !is.na(c_col)) {
      out[[length(out) + 1]] <- tibble::tibble(
        outcome = canon_outcome(X[[c_col]]),
        term_raw = as.character(X[[a_col]]),
        term_norm = norm_term(X[[a_col]]),
        display = format_term_display(X[[a_col]]),
        support = num_from_choices(X, c("articles_weighted", "Weighted articles")),
        rank_value = NA_real_,
        aux_value = abs(num_from_choices(X, c("weighted_balance", "Weighted directional balance"))),
        source_type = "Internal weighted signed-effect summary",
        source_file = f_weighted,
        source_priority = 3,
        metric_used = "articles_weighted + |weighted balance|"
      )
    }
  }

  if (!is.na(f_unweighted)) {
    X <- readr::read_csv(f_unweighted, show_col_types = FALSE)
    a_col <- first_existing_col(X, c("A", "A term"))
    c_col <- first_existing_col(X, c("C", "Outcome"))
    if (!is.na(a_col) && !is.na(c_col)) {
      out[[length(out) + 1]] <- tibble::tibble(
        outcome = canon_outcome(X[[c_col]]),
        term_raw = as.character(X[[a_col]]),
        term_norm = norm_term(X[[a_col]]),
        display = format_term_display(X[[a_col]]),
        support = dplyr::coalesce(
          num_from_choices(X, c("articles", "Unweighted articles, n")),
          num_from_choices(X, c("directional_n"))
        ),
        rank_value = NA_real_,
        aux_value = abs(num_from_choices(X, c("balance", "Unweighted directional balance"))),
        source_type = "Unweighted signed-effect summary",
        source_file = f_unweighted,
        source_priority = 4,
        metric_used = "articles + |balance|"
      )
    }
  }

  if (!length(out)) return(NULL)

  X <- dplyr::bind_rows(out) |>
    dplyr::filter(!is.na(outcome), !is.na(term_norm), term_norm != "")

  if (!is.null(dic_lookup) && nrow(dic_lookup)) {
    X <- X |>
      dplyr::left_join(
        dic_lookup |>
          dplyr::select(term_norm, class_low, role_up, dic_path),
        by = "term_norm"
      )
  } else {
    X <- X |>
      dplyr::mutate(class_low = NA_character_, role_up = NA_character_, dic_path = NA_character_)
  }

  X |>
    dplyr::filter(!term_norm %in% HOTSPOT_DROP_GENERIC) |>
    dplyr::filter(is.na(class_low) | !class_low %in% HOTSPOT_EXCLUDE_CLASSES) |>
    dplyr::mutate(
      support = suppressWarnings(as.numeric(support)),
      rank_value = suppressWarnings(as.numeric(rank_value)),
      aux_value = suppressWarnings(as.numeric(aux_value)),
      preferred = dplyr::case_when(
        class_low %in% c("pattern", "radiology", "airway", "qct", "pft", "activity", "phenotype", "complication", "trend") ~ TRUE,
        class_low %in% c("host", "exposure", "event", "population", "system", "strategy", "nonpharm") ~ TRUE,
        TRUE ~ FALSE
      ),
      class_priority = dplyr::case_when(
        class_low %in% c("pattern", "radiology", "airway", "qct", "pft", "activity", "phenotype", "complication", "trend") ~ 1L,
        class_low %in% c("host", "exposure", "event", "population", "system", "strategy", "nonpharm") ~ 2L,
        is.na(class_low) ~ 3L,
        TRUE ~ 4L
      ),
      is_family_like = stringr::str_detect(term_norm, "_FAMILY$") | stringr::str_detect(term_norm, "PATHWAY|CHECKPOINT|INFLAMMATION|STRESS")
    ) |>
    dplyr::arrange(outcome, source_priority, class_priority, dplyr::desc(support), dplyr::desc(aux_value), rank_value, display) |>
    dplyr::group_by(outcome, term_norm) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
}

prep_biomarker_candidates <- function(){
  f_s4 <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("Supplementary_Table_S4_Design_tier_weighted_biomarker_support_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^Supplementary_Table_S4_Design_tier_weighted_biomarker_support_.*\\.csv$")
  ))

  f_biom_w <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("biomarker_metrics_design_weighted_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^biomarker_metrics_design_weighted_.*\\.csv$")
  ))

  f_atlas <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("biomarker_outcome_atlas_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^biomarker_outcome_atlas_.*\\.csv$")
  ))

  out <- list()

  if (!is.na(f_s4)) {
    S4 <- readr::read_csv(f_s4, show_col_types = FALSE)
    need <- c("Biomarker", "Outcome")
    if (all(need %in% names(S4))) {
      pmid_col <- if ("Weighted PMID support" %in% names(S4)) "Weighted PMID support" else NA_character_
      rank_col <- if ("Weighted rank" %in% names(S4)) "Weighted rank" else NA_character_
      auc_col  <- if ("Mean AUC (weighted)" %in% names(S4)) "Mean AUC (weighted)" else NA_character_

      out[[length(out) + 1]] <- S4 |>
        dplyr::transmute(
          outcome = canon_outcome(.data$Outcome),
          marker_raw = .data$Biomarker,
          marker_norm = norm_term(.data$Biomarker),
          display = format_term_display(.data$Biomarker),
          support = if (!is.na(pmid_col)) suppressWarnings(as.numeric(.data[[pmid_col]])) else NA_real_,
          rank_value = if (!is.na(rank_col)) suppressWarnings(as.numeric(.data[[rank_col]])) else NA_real_,
          aux_value = if (!is.na(auc_col)) suppressWarnings(as.numeric(.data[[auc_col]])) else NA_real_,
          source_type = "Supplementary Table S4",
          source_file = f_s4,
          source_priority = 1,
          metric_used = "Weighted PMID support"
        )
    }
  }

  if (!is.na(f_biom_w)) {
    BW <- readr::read_csv(f_biom_w, show_col_types = FALSE)
    need <- c("marker", "C")
    if (all(need %in% names(BW))) {
      out[[length(out) + 1]] <- BW |>
        dplyr::transmute(
          outcome = canon_outcome(.data$C),
          marker_raw = .data$marker,
          marker_norm = norm_term(.data$marker),
          display = format_term_display(.data$marker),
          support = if ("pmids_weighted" %in% names(BW)) suppressWarnings(as.numeric(.data$pmids_weighted)) else NA_real_,
          rank_value = if ("weighted_rank" %in% names(BW)) suppressWarnings(as.numeric(.data$weighted_rank)) else NA_real_,
          aux_value = if ("auc_weighted_mean" %in% names(BW)) suppressWarnings(as.numeric(.data$auc_weighted_mean)) else NA_real_,
          source_type = "Internal weighted biomarker support",
          source_file = f_biom_w,
          source_priority = 2,
          metric_used = "pmids_weighted"
        )
    }
  }

  if (!is.na(f_atlas)) {
    AT <- readr::read_csv(f_atlas, show_col_types = FALSE)
    need <- c("A", "C")
    if (all(need %in% names(AT))) {
      out[[length(out) + 1]] <- AT |>
        dplyr::transmute(
          outcome = canon_outcome(.data$C),
          marker_raw = .data$A,
          marker_norm = norm_term(.data$A),
          display = format_term_display(.data$A),
          support = dplyr::coalesce(
            if ("evidence" %in% names(AT)) suppressWarnings(as.numeric(.data$evidence)) else NA_real_,
            if ("articles" %in% names(AT)) suppressWarnings(as.numeric(.data$articles)) else NA_real_
          ),
          rank_value = NA_real_,
          aux_value = if ("balance" %in% names(AT)) suppressWarnings(as.numeric(.data$balance)) else NA_real_,
          source_type = "Biomarker outcome atlas",
          source_file = f_atlas,
          source_priority = 3,
          metric_used = if ("evidence" %in% names(AT)) "evidence" else "articles"
        )
    }
  }

  if (!length(out)) return(NULL)

  dplyr::bind_rows(out) |>
    dplyr::filter(!is.na(outcome), !is.na(marker_norm), marker_norm != "") |>
    dplyr::filter(!marker_norm %in% LOWER_DROP_GENERIC) |>
    dplyr::mutate(
      support = suppressWarnings(as.numeric(support)),
      rank_value = suppressWarnings(as.numeric(rank_value)),
      aux_value = suppressWarnings(as.numeric(aux_value)),
      preferred = !marker_norm %in% LOWER_DROP_PREFER,
      is_family_like = stringr::str_detect(marker_norm, "_FAMILY$")
    ) |>
    dplyr::arrange(outcome, source_priority, dplyr::desc(preferred), is_family_like, dplyr::desc(support), rank_value, dplyr::desc(aux_value), display) |>
    dplyr::group_by(outcome, marker_norm) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
}


prep_hotspot_candidates_from_cooc <- function(dic_lookup){
  f_cooc <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("cooc_npmi_lift_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^cooc_npmi_lift_.*\\.csv$")
  ))
  if (is.na(f_cooc)) return(NULL)

  X <- readr::read_csv(f_cooc, show_col_types = FALSE)

  term1_col <- first_existing_col(X, c("term1", "x", "lhs", "A", "node1", "term_x", "Var1"))
  term2_col <- first_existing_col(X, c("term2", "y", "rhs", "B", "C", "node2", "term_y", "Var2"))
  class1_col <- first_existing_col(X, c("class1", "class_x", "lhs_class", "A_class", "node1_class"))
  class2_col <- first_existing_col(X, c("class2", "class_y", "rhs_class", "B_class", "C_class", "node2_class"))

  if (is.na(term1_col) || is.na(term2_col)) return(NULL)

  t1 <- as.character(X[[term1_col]])
  t2 <- as.character(X[[term2_col]])
  c1 <- if (!is.na(class1_col)) tolower(as.character(X[[class1_col]])) else rep(NA_character_, nrow(X))
  c2 <- if (!is.na(class2_col)) tolower(as.character(X[[class2_col]])) else rep(NA_character_, nrow(X))

  support_vec <- dplyr::coalesce(
    num_from_choices(X, c("n11", "n_ac", "nAC", "count", "cooc", "co_mentions", "weighted_n11", "n11_weighted")),
    abs(num_from_choices(X, c("or", "OR")))
  )
  aux_vec <- dplyr::coalesce(
    num_from_choices(X, c("npmi", "NPMI")),
    num_from_choices(X, c("lift", "Lift")),
    abs(num_from_choices(X, c("or", "OR")))
  )
  rank_vec <- dplyr::coalesce(
    num_from_choices(X, c("q", "q_value", "fdr", "padj")),
    num_from_choices(X, c("p", "p_value"))
  )

  Y1 <- tibble::tibble(
    outcome = canon_outcome(t1),
    term_raw = t2,
    term_norm = norm_term(t2),
    display = format_term_display(t2),
    support = support_vec,
    rank_value = rank_vec,
    aux_value = aux_vec,
    source_type = "Pairwise co-mention summary",
    source_file = f_cooc,
    source_priority = 5,
    metric_used = "n11 + NPMI/lift",
    class_low = c2,
    role_up = NA_character_,
    dic_path = NA_character_
  )

  Y2 <- tibble::tibble(
    outcome = canon_outcome(t2),
    term_raw = t1,
    term_norm = norm_term(t1),
    display = format_term_display(t1),
    support = support_vec,
    rank_value = rank_vec,
    aux_value = aux_vec,
    source_type = "Pairwise co-mention summary",
    source_file = f_cooc,
    source_priority = 5,
    metric_used = "n11 + NPMI/lift",
    class_low = c1,
    role_up = NA_character_,
    dic_path = NA_character_
  )

  Y <- dplyr::bind_rows(Y1, Y2) |>
    dplyr::filter(!is.na(outcome), !is.na(term_norm), term_norm != "") |>
    dplyr::filter(is.na(canon_outcome(term_raw)))

  if (!nrow(Y)) return(NULL)

  if (!is.null(dic_lookup) && nrow(dic_lookup)) {
    Y <- Y |>
      dplyr::left_join(
        dic_lookup |>
          dplyr::select(term_norm, class_low_dic = class_low, role_up_dic = role_up, dic_path_dic = dic_path),
        by = "term_norm"
      ) |>
      dplyr::mutate(
        class_low = dplyr::coalesce(class_low, class_low_dic),
        role_up = dplyr::coalesce(role_up, role_up_dic),
        dic_path = dplyr::coalesce(dic_path, dic_path_dic)
      ) |>
      dplyr::select(-class_low_dic, -role_up_dic, -dic_path_dic)
  }

  Y |>
    dplyr::filter(!term_norm %in% HOTSPOT_DROP_GENERIC) |>
    dplyr::filter(is.na(class_low) | !class_low %in% HOTSPOT_EXCLUDE_CLASSES) |>
    dplyr::mutate(
      support = suppressWarnings(as.numeric(support)),
      rank_value = suppressWarnings(as.numeric(rank_value)),
      aux_value = suppressWarnings(as.numeric(aux_value)),
      preferred = dplyr::case_when(
        class_low %in% c("pattern", "radiology", "airway", "qct", "pft", "activity", "phenotype", "complication", "trend") ~ TRUE,
        class_low %in% c("host", "exposure", "event", "population", "system", "strategy", "nonpharm") ~ TRUE,
        TRUE ~ FALSE
      ),
      class_priority = dplyr::case_when(
        class_low %in% c("pattern", "radiology", "airway", "qct", "pft", "activity", "phenotype", "complication", "trend") ~ 1L,
        class_low %in% c("host", "exposure", "event", "population", "system", "strategy", "nonpharm") ~ 2L,
        is.na(class_low) ~ 3L,
        TRUE ~ 4L
      ),
      is_family_like = stringr::str_detect(term_norm, "_FAMILY$") | stringr::str_detect(term_norm, "PATHWAY|CHECKPOINT|INFLAMMATION|STRESS")
    ) |>
    dplyr::arrange(outcome, source_priority, class_priority, dplyr::desc(support), dplyr::desc(aux_value), rank_value, display) |>
    dplyr::group_by(outcome, term_norm) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
}

pick_hotspot_terms <- function(candidates, outcome_i, n_terms = 8){
  X <- candidates |>
    dplyr::filter(outcome == outcome_i)

  if (!nrow(X)) return(X)

  preferred <- X |>
    dplyr::filter(preferred, !is_family_like) |>
    dplyr::arrange(source_priority, class_priority, dplyr::desc(support), dplyr::desc(aux_value), rank_value, display)

  if (nrow(preferred) >= n_terms) {
    return(dplyr::slice_head(preferred, n = n_terms))
  }

  fallback <- X |>
    dplyr::arrange(source_priority, class_priority, dplyr::desc(preferred), is_family_like, dplyr::desc(support), dplyr::desc(aux_value), rank_value, display)

  dplyr::bind_rows(preferred, fallback) |>
    dplyr::distinct(term_norm, .keep_all = TRUE) |>
    dplyr::slice_head(n = n_terms)
}

pick_lower_terms <- function(candidates, outcome_i, n_terms = 8){
  X <- candidates |>
    dplyr::filter(outcome == outcome_i)

  if (!nrow(X)) return(X)

  preferred <- X |>
    dplyr::filter(preferred, !is_family_like) |>
    dplyr::arrange(source_priority, dplyr::desc(support), rank_value, dplyr::desc(aux_value), display)

  if (nrow(preferred) >= n_terms) {
    return(dplyr::slice_head(preferred, n = n_terms))
  }

  fallback <- X |>
    dplyr::arrange(source_priority, dplyr::desc(preferred), is_family_like, dplyr::desc(support), rank_value, dplyr::desc(aux_value), display)

  dplyr::bind_rows(preferred, fallback) |>
    dplyr::distinct(marker_norm, .keep_all = TRUE) |>
    dplyr::slice_head(n = n_terms)
}

# -------------------------
# Build data-supported layers
# -------------------------
dic_lookup <- build_dic_lookup()

hotspot_candidates <- prep_hotspot_candidates(dic_lookup)
cooc_hotspot_candidates <- prep_hotspot_candidates_from_cooc(dic_lookup)

if (is.null(hotspot_candidates) || !nrow(hotspot_candidates)) {
  hotspot_candidates <- cooc_hotspot_candidates
} else if (!is.null(cooc_hotspot_candidates) && nrow(cooc_hotspot_candidates)) {
  hotspot_candidates <- dplyr::bind_rows(hotspot_candidates, cooc_hotspot_candidates) |>
    dplyr::arrange(outcome, source_priority, class_priority, dplyr::desc(support), dplyr::desc(aux_value), rank_value, display) |>
    dplyr::group_by(outcome, term_norm) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
}

if (is.null(hotspot_candidates) || !nrow(hotspot_candidates)) {
  stop("No signed-effect or pairwise co-mention table was found for data-driven hot-spot generation in Figure 4.")
}

biomarker_candidates <- prep_biomarker_candidates()
if (is.null(biomarker_candidates) || !nrow(biomarker_candidates)) {
  stop("No biomarker support table was found for data-driven biomarker-layer generation in Figure 4.")
}

selected_hotspots <- dplyr::bind_rows(lapply(OUTCOME_LEVELS, function(oc) pick_hotspot_terms(hotspot_candidates, oc, n_terms = HOTSPOT_TERMS_PER_OUTCOME))) |>
  dplyr::group_by(outcome) |>
  dplyr::mutate(selection_order = dplyr::row_number()) |>
  dplyr::ungroup()

selected_lower <- dplyr::bind_rows(lapply(OUTCOME_LEVELS, function(oc) pick_lower_terms(biomarker_candidates, oc, n_terms = LOWER_TERMS_PER_OUTCOME))) |>
  dplyr::group_by(outcome) |>
  dplyr::mutate(selection_order = dplyr::row_number()) |>
  dplyr::ungroup()

hotspot_counts <- selected_hotspots |>
  dplyr::count(outcome, name = "selected_hotspots_n")
lower_counts <- selected_lower |>
  dplyr::count(outcome, name = "selected_lower_n")

missing_hotspot_outcomes <- setdiff(OUTCOME_LEVELS, unique(selected_hotspots$outcome))
missing_lower_outcomes <- setdiff(OUTCOME_LEVELS, unique(selected_lower$outcome))
if (length(missing_hotspot_outcomes)) {
  msg <- paste0("Data-driven hot-spot selection failed for outcome(s): ", paste(missing_hotspot_outcomes, collapse = ", "))
  if (STRICT_DATA_DRIVEN) stop(msg) else warning(msg)
}
if (length(missing_lower_outcomes)) {
  msg <- paste0("Data-driven biomarker-layer selection failed for outcome(s): ", paste(missing_lower_outcomes, collapse = ", "))
  if (STRICT_DATA_DRIVEN) stop(msg) else warning(msg)
}

hotspot_map <- selected_hotspots |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(
    hotspots = collapse_terms(display, per_line = HOTSPOT_TERMS_PER_LINE),
    hotspot_source = paste(unique(source_type), collapse = " | "),
    hotspot_n = dplyr::n(),
    .groups = "drop"
  )

lower_map <- selected_lower |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(
    lower = collapse_terms(display, per_line = LOWER_TERMS_PER_LINE),
    lower_source = paste(unique(source_type), collapse = " | "),
    lower_n = dplyr::n(),
    .groups = "drop"
  )

cards <- CARD_FILL |>
  dplyr::left_join(hotspot_map, by = "outcome") |>
  dplyr::left_join(lower_map, by = "outcome") |>
  dplyr::mutate(
    hotspots = ifelse(is.na(hotspots), "NA", hotspots),
    hotspot_source = ifelse(is.na(hotspot_source), "NA", hotspot_source),
    lower = ifelse(is.na(lower), "NA", lower),
    lower_source = ifelse(is.na(lower_source), "NA", lower_source)
  )

# -------------------------
# Plot
# -------------------------
make_label_panel <- function(){
  ggplot2::ggplot() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::annotate(
      "text",
      x = 0.03, y = 0.67, label = "Hot spots",
      hjust = 0, vjust = 0.5, size = 6.1
    ) +
    ggplot2::annotate(
      "text",
      x = 0.03, y = 0.16, label = "Biomarker / exploratory layer",
      hjust = 0, vjust = 0.5, size = 5.0
    )
}

make_card <- function(title, hotspots_txt, lower_txt, header_fill){
  ggplot2::ggplot() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::annotate(
      "label",
      x = 0.50, y = 0.935, label = title,
      size = 5.0, fontface = "bold",
      fill = header_fill, colour = "white",
      label.size = 0, label.padding = unit(0.38, "lines")
    ) +
    ggplot2::annotate(
      "label",
      x = 0.50, y = 0.65,
      label = hotspots_txt,
      size = 3.10, lineheight = 1.00,
      fill = "white", colour = "black",
      label.size = 0.35, label.padding = unit(0.44, "lines")
    ) +
    ggplot2::annotate(
      "label",
      x = 0.50, y = 0.15,
      label = lower_txt,
      size = 2.75, lineheight = 1.00,
      fill = scales::alpha(header_fill, 0.08), colour = "black",
      label.size = 0.25, label.padding = unit(0.32, "lines")
    )
}

label_panel <- make_label_panel()
plots <- lapply(seq_len(nrow(cards)), function(i){
  make_card(cards$outcome[i], cards$hotspots[i], cards$lower[i], cards$fill[i])
})

p <- label_panel + plots[[1]] + plots[[2]] + plots[[3]] +
  patchwork::plot_layout(widths = c(1.20, 1, 1, 1)) +
  patchwork::plot_annotation(
    title = "Outcome-specific hot spot map in RA-ILD",
    subtitle = "Hot spots were selected from signed-effect and pairwise co-mention outputs; biomarker/exploratory terms were selected from biomarker support outputs",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, margin = margin(b = 2)),
      plot.margin = margin(4, 10, 2, 12)
    )
  )

# -------------------------
# Save
# -------------------------
FIG_DIR <- if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG")) DIR_FIG else getwd()
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)

f_png <- file.path(FIG_DIR, "Figure4_outcome_specific_hotspot_map_fully_datadriven_v5.png")
f_pdf <- file.path(FIG_DIR, "Figure4_outcome_specific_hotspot_map_fully_datadriven_v5.pdf")
f_labels <- file.path(DIR_TABLE, sprintf("Figure4_outcome_specific_hotspot_map_labels_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_hotspot_prov <- file.path(DIR_TABLE, sprintf("Figure4_hotspot_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_biom_prov <- file.path(DIR_TABLE, sprintf("Figure4_biomarker_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG))

ggplot2::ggsave(f_png, p, width = 15.0, height = 8.8, dpi = 320, bg = "white")
ggplot2::ggsave(f_pdf, p, width = 15.0, height = 8.8, bg = "white")

readr::write_csv(cards, f_labels)
readr::write_csv(
  selected_hotspots |>
    dplyr::select(
      outcome,
      selection_order,
      display,
      term_raw,
      term_norm,
      class_low,
      role_up,
      support,
      rank_value,
      aux_value,
      metric_used,
      source_type,
      source_file,
      dic_path
    ),
  f_hotspot_prov
)
readr::write_csv(
  selected_lower |>
    dplyr::select(
      outcome,
      selection_order,
      display,
      marker_raw,
      marker_norm,
      support,
      rank_value,
      aux_value,
      metric_used,
      source_type,
      source_file
    ),
  f_biom_prov
)

message("Saved Figure 4 with fully data-supported layers (no Dominant theme):")
message(" - figure PNG: ", f_png)
message(" - figure PDF: ", f_pdf)
message(" - labels CSV: ", f_labels)
message(" - hotspot provenance CSV: ", f_hotspot_prov)
message(" - biomarker provenance CSV: ", f_biom_prov)
message("Selected hot-spot counts by outcome:")
print(hotspot_counts)
message("Selected biomarker-layer counts by outcome:")
print(lower_counts)
message("Selected hot spots by outcome:")
print(selected_hotspots |>
        dplyr::select(outcome, selection_order, display, class_low, source_type, support, aux_value))
message("Selected biomarker/exploratory terms by outcome:")
print(selected_lower |>
        dplyr::select(outcome, selection_order, display, source_type, support, rank_value))
