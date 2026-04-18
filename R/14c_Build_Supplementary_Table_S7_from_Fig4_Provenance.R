# =========================================================
# 14c_Build_Supplementary_Table_S7_from_Fig4_Provenance.R
# ---------------------------------------------------------
# Builds Supplementary Table S7 from current Figure 4 provenance outputs
#   - uses 14a provenance CSVs
#   - applies 14b-style sharedness logic
#   - writes manuscript-ready Panel A / Panel B / combined CSVs
# =========================================================

.find_setup <- function() {
  cand <- c(
    file.path(getwd(), "00_setup_Final.R"),
    file.path(getwd(), "00_setup.R"),
    file.path(Sys.getenv("RAILD_ROOT", unset = ""), "00_setup_Final.R"),
    file.path(Sys.getenv("RAILD_ROOT", unset = ""), "00_setup.R")
  )
  for (p in cand) {
    if (nzchar(p) && file.exists(p)) return(p)
  }
  stop("Cannot find 00_setup_Final.R or 00_setup.R.")
}
source(.find_setup())

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("readr","dplyr","stringr","tibble","purrr"))

find_first_existing <- function(paths){
  ok <- paths[!is.na(paths) & file.exists(paths)]
  if (length(ok)) ok[1] else NA_character_
}

find_latest_file <- function(dir, regex){
  if (is.na(dir) || !dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

canon_outcome <- function(x){
  xl <- stringr::str_to_lower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(xl, "^ae-ild$|acute") ~ "AE-ILD",
    stringr::str_detect(xl, "progress") ~ "Progression",
    stringr::str_detect(xl, "mort") ~ "Mortality",
    TRUE ~ as.character(x)
  )
}

clean_display <- function(x){
  x <- as.character(x)
  dplyr::case_when(
    x == "ACPA RF high" ~ "ACPA/RF high",
    x == "ACPA_RF high" ~ "ACPA/RF high",
    x == "SP A" ~ "SP-A",
    x == "SP_A" ~ "SP-A",
    x == "GM CSF" ~ "GM-CSF",
    x == "GM_CSF" ~ "GM-CSF",
    x == "IL8 CXCL8" ~ "IL-8/CXCL8",
    x == "IL8_CXCL8" ~ "IL-8/CXCL8",
    x == "WNT beta catenin" ~ "WNT/β-catenin",
    x == "WNT_beta_catenin" ~ "WNT/β-catenin",
    x == "drug induced ILD" ~ "drug-induced ILD",
    TRUE ~ x
  )
}

norm_key <- function(x){
  x |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_replace_all("–", "-") |>
    stringr::str_replace_all("[[:space:]/-]+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9_]+", "") |>
    toupper()
}

make_outcomes_present <- function(x, outcome_order){
  present <- outcome_order[outcome_order %in% unique(as.character(x))]
  paste(present, collapse = "; ")
}

source_code_hot <- function(source_type, source_file){
  s <- tolower(paste(source_type, basename(source_file)))
  dplyr::case_when(
    stringr::str_detect(s, "co-mention|cooc") ~ "Cooc",
    TRUE ~ "SE-wt"
  )
}

source_code_bio <- function(source_type, source_file){
  s <- tolower(paste(source_type, basename(source_file)))
  dplyr::case_when(
    stringr::str_detect(s, "supplementary table s4|weighted_biomarker_support|biomarker_metrics_design_weighted") ~ "S4",
    TRUE ~ "Atlas"
  )
}

round_or_na <- function(x, digits = 3){
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_real_, round(x, digits))
}

build_panel_table <- function(df, panel = c("A","B")){
  panel <- match.arg(panel)

  OUTCOME_ORDER <- c("AE-ILD", "Progression", "Mortality")

  if (panel == "A") {
    canonical_col <- if ("term_raw" %in% names(df)) "term_raw" else if ("term_norm" %in% names(df)) "term_norm" else "display"
    class_col <- if ("class_low" %in% names(df)) "class_low" else NA_character_
    source_fun <- source_code_hot
    preferred_shared <- c(
      "Male sex","Smoking","disease duration","drug-induced ILD",
      "Traction bronchiectasis","UIP","CPI index","Honeycombing"
    )
  } else {
    canonical_col <- if ("marker_raw" %in% names(df)) "marker_raw" else if ("marker_norm" %in% names(df)) "marker_norm" else "display"
    class_col <- NA_character_
    source_fun <- source_code_bio
    preferred_shared <- c("KL-6","SP-D","ACPA/RF","MMP-7","MPO")
  }

  X <- df |>
    dplyr::mutate(
      outcome = canon_outcome(outcome),
      display_clean = clean_display(display),
      canonical_term = as.character(.data[[canonical_col]]),
      canonical_term = ifelse(is.na(canonical_term) | canonical_term == "", display_clean, canonical_term),
      key = norm_key(canonical_term),
      source_code = source_fun(source_type, source_file),
      class_val = if (!is.na(class_col)) as.character(.data[[class_col]]) else "",
      support_num = round_or_na(support, 3),
      rank_num = round_or_na(rank_value, 3),
      aux_num = round_or_na(aux_value, 3)
    ) |>
    dplyr::filter(outcome %in% OUTCOME_ORDER, !is.na(key), key != "")

  counts <- X |>
    dplyr::group_by(key) |>
    dplyr::summarise(
      n_outcomes = dplyr::n_distinct(outcome),
      outcomes_present = make_outcomes_present(outcome, OUTCOME_ORDER),
      .groups = "drop"
    )

  X <- X |>
    dplyr::left_join(counts, by = "key") |>
    dplyr::mutate(
      Sharedness = dplyr::case_when(
        n_outcomes >= 3 ~ "All 3",
        n_outcomes == 2 ~ "2 outcomes",
        TRUE ~ "Outcome-specific"
      ),
      shared_order = match(display_clean, preferred_shared),
      shared_order = ifelse(is.na(shared_order), 999L, shared_order),
      sharedness_order = dplyr::case_when(
        Sharedness == "All 3" ~ 1L,
        Sharedness == "2 outcomes" ~ 2L,
        TRUE ~ 3L
      ),
      outcome_order = match(outcome, OUTCOME_ORDER)
    ) |>
    dplyr::arrange(
      sharedness_order,
      shared_order,
      outcome_order,
      selection_order,
      display_clean
    ) |>
    dplyr::transmute(
      Panel = ifelse(panel == "A",
                     "Panel A. Hot-spot layer provenance",
                     "Panel B. Biomarker / exploratory layer provenance"),
      Sharedness = Sharedness,
      `Outcomes present` = outcomes_present,
      Outcome = outcome,
      Order = selection_order,
      `Display term` = display_clean,
      `Canonical term` = canonical_term,
      Class = ifelse(panel == "A", class_val, ""),
      Support = support_num,
      Rank = rank_num,
      Aux = aux_num,
      Source = source_code
    )

  X
}

build_supplementary_table_S7_from_fig4_provenance <- function(verbose = TRUE){

  f_hot <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("Figure4_hotspot_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^Figure4_hotspot_layer_provenance_fully_datadriven_v5_.*\\.csv$")
  ))

  f_bio <- find_first_existing(c(
    file.path(DIR_TABLE, sprintf("Figure4_biomarker_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
    find_latest_file(DIR_TABLE, "^Figure4_biomarker_layer_provenance_fully_datadriven_v5_.*\\.csv$")
  ))

  if (is.na(f_hot) || is.na(f_bio)) {
    stop("Missing Figure 4 provenance CSV(s). Run 14a first.")
  }

  HOT <- readr::read_csv(f_hot, show_col_types = FALSE)
  BIO <- readr::read_csv(f_bio, show_col_types = FALSE)

  if (!all(c("outcome","selection_order","display","source_type","source_file") %in% names(HOT))) {
    stop("Hotspot provenance CSV is missing required columns.")
  }
  if (!all(c("outcome","selection_order","display","source_type","source_file") %in% names(BIO))) {
    stop("Biomarker provenance CSV is missing required columns.")
  }

  S7A <- build_panel_table(HOT, panel = "A")
  S7B <- build_panel_table(BIO, panel = "B")
  S7  <- dplyr::bind_rows(S7A, S7B)

  f_s7_a <- file.path(
    DIR_TABLE,
    sprintf("Supplementary_Table_S7_PanelA_Figure4_hotspot_provenance_%s__%s.csv", CORPUS_TAG, DIC_TAG)
  )
  f_s7_b <- file.path(
    DIR_TABLE,
    sprintf("Supplementary_Table_S7_PanelB_Figure4_biomarker_provenance_%s__%s.csv", CORPUS_TAG, DIC_TAG)
  )
  f_s7 <- file.path(
    DIR_TABLE,
    sprintf("Supplementary_Table_S7_Figure4_term_selection_provenance_%s__%s.csv", CORPUS_TAG, DIC_TAG)
  )
  f_legend <- file.path(
    DIR_TABLE,
    sprintf("Supplementary_Table_S7_legend_%s__%s.txt", CORPUS_TAG, DIC_TAG)
  )

  readr::write_csv(S7A, f_s7_a, na = "")
  readr::write_csv(S7B, f_s7_b, na = "")
  readr::write_csv(S7,  f_s7,   na = "")

  legend_lines <- c(
    "Supplementary Table S7. Figure 4 term-selection provenance.",
    "Term-level provenance for the fully data-driven Figure 4 pipeline, regrouped by cross-outcome sharedness.",
    "",
    "Panel A. Hot-spot layer provenance.",
    "Outcome-specific hot-spot terms selected from signed-effect and pairwise co-mention outputs.",
    "",
    "Panel B. Biomarker / exploratory layer provenance.",
    "Outcome-specific biomarker / exploratory terms selected from biomarker support outputs.",
    "",
    "Abbreviations and source codes.",
    "SE-wt = signed-effects weighted vs unweighted comparison or weighted signed-effect-derived source.",
    "Cooc = pairwise co-mention summary.",
    "S4 = Supplementary Table S4 weighted biomarker support.",
    "Atlas = biomarker outcome atlas.",
    "Support, rank, and aux are the numeric fields exported by the Figure 4 fully data-driven pipeline from the corresponding upstream source tables.",
    "Sharedness was added downstream from the selected Figure 4 term lists."
  )
  writeLines(legend_lines, con = f_legend)

  if (isTRUE(verbose)) {
    message("Using hotspot provenance : ", f_hot)
    message("Using biomarker provenance: ", f_bio)
    message("WROTE: ", f_s7_a)
    message("WROTE: ", f_s7_b)
    message("WROTE: ", f_s7)
    message("WROTE: ", f_legend)
  }

  invisible(list(
    hotspot_provenance = f_hot,
    biomarker_provenance = f_bio,
    s7_panel_a = f_s7_a,
    s7_panel_b = f_s7_b,
    s7_combined = f_s7,
    legend = f_legend
  ))
}

# 実行
res_s7 <- build_supplementary_table_S7_from_fig4_provenance(verbose = TRUE)