# =========================================================
# 08_Sensitivity_Nonreview_Final.R  (public-ready, no-mixing design)
# - Sensitivity analysis excluding Review / Case Reports
# - Outputs drug x outcome co-occurrence table (counts + ratio)
# =========================================================

# ---- bootstrap: locate project setup (portable paths) ----
setup_candidates <- c(
  file.path(getwd(), "00_setup_Final.R"),
  file.path(Sys.getenv("RAILD_ROOT"), "00_setup_Final.R")
)
setup_path <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_path) || setup_path == "") {
  stop("Cannot find 00_setup_Final.R. Run from project root or set RAILD_ROOT to the project folder.")
}
source(setup_path)



quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","tibble","purrr"))

# -------------------------
# Parameters
# -------------------------
OUTCOMES_FOCUS <- c("AE-ILD","progression","mortality","FVC_decline","DLCO_decline")
MIN_A_MENTIONS <- 10   # min mentions of A in subset to report ratio
MIN_AC         <- 3    # min co-occurrence to report ratio

# -------------------------
# Load inputs (STRICT, tagged)
# -------------------------
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_hits_tag)) {
  stop("Tagged hits_matrix not found: ", f_hits_tag,
       "\nRun 02_build_hits_matrix.R with tagged saving enabled.")
}
HM <- readr::read_csv(f_hits_tag, show_col_types = FALSE)
log_msg("08 loaded hits_matrix:", f_hits_tag, " n=", nrow(HM))

dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
DIC <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","class") %in% names(DIC))) stop("Dictionary must have term/class columns.")
log_msg("08 using dictionary:", basename(dic_path), " n=", nrow(DIC))

# -------------------------
# Preconditions
# -------------------------
need_cols <- c("pmid","pubtype")
miss <- setdiff(need_cols, names(HM))
if (length(miss)) stop("hits_matrix missing required columns: ", paste(miss, collapse=", "))

binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

# -------------------------
# 1) Subset: exclude Review / Case Reports
# -------------------------
NR <- HM %>%
  dplyr::filter(
    !grepl("(?i)Review", pubtype),
    !grepl("(?i)Case\\s*Reports?", pubtype)
  )

log_msg("08 non-review/case subset n=", nrow(NR), " (from ", nrow(HM), ")")

if (nrow(NR) < 50) {
  log_msg("08: subset too small (<50). Still writing empty outputs for reproducibility.")
}

# -------------------------
# 2) Define drugs and outcomes
# -------------------------
hit_terms <- sub("^hit__", "", names(NR)[startsWith(names(NR), "hit__")])

drugs <- intersect(DIC$term[DIC$class == "drug"], hit_terms)
outcomes_all <- intersect(DIC$term[DIC$class == "outcome"], hit_terms)

# Use focus outcomes if present; otherwise fallback to all outcomes
outcomes <- intersect(OUTCOMES_FOCUS, outcomes_all)
if (!length(outcomes)) outcomes <- outcomes_all

log_msg("08 drugs n=", length(drugs), " outcomes n=", length(outcomes))

if (!length(drugs) || !length(outcomes)) {
  # Write empty result (public-ready behavior)
  out_empty <- tibble::tibble(A=character(), C=character(), nA=integer(), nC=integer(), nAC=integer(), ratio=double())
  f_out <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  readr::write_csv(out_empty, f_out)
  log_msg("WROTE (empty):", f_out)
  log_msg("=== DONE 08_sensitivity_nonreview ===")
  quit(save="no", status=0)
}

# -------------------------
# 3) Compute co-occurrence + ratio  (FIXED: always creates nA/nC/nAC/ratio)
# -------------------------
get_col <- function(term){
  nm <- paste0("hit__", term)
  if (nm %in% names(NR)) binv(NR[[nm]]) else integer(nrow(NR))
}

out <- purrr::map_dfr(drugs, function(a){
  vA <- get_col(a)
  nA <- as.integer(sum(vA, na.rm = TRUE))

  purrr::map_dfr(outcomes, function(c){
    vC  <- get_col(c)
    nC  <- as.integer(sum(vC, na.rm = TRUE))
    nAC <- as.integer(sum(vA & vC, na.rm = TRUE))
    ratio <- if (nA > 0) (nAC / nA) else NA_real_

    tibble::tibble(
      A = a, C = c,
      nA = nA, nC = nC, nAC = nAC,
      ratio = ratio
    )
  })
})

# Safety check: stop here if expected columns are missing, to make the cause explicit
need_cols <- c("A","C","nA","nC","nAC","ratio")
miss <- setdiff(need_cols, names(out))
if (length(miss)) stop("Internal error: missing columns in out: ", paste(miss, collapse=", "))

# filter for readability
out_main <- out %>%
  dplyr::filter(nA >= MIN_A_MENTIONS, nAC >= MIN_AC) %>%
  dplyr::arrange(dplyr::desc(nAC), dplyr::desc(ratio), dplyr::desc(nA))


# -------------------------
# 4) Output (tagged)
# -------------------------
f_full <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_full_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_main <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))

readr::write_csv(out, f_full)
readr::write_csv(out_main, f_main)

log_msg("WROTE:", f_full, " n=", nrow(out))
log_msg("WROTE:", f_main, " n=", nrow(out_main))
log_msg("=== DONE 08_sensitivity_nonreview ===")