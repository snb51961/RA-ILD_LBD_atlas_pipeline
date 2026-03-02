# =========================================================
# 05_ac_npmi.R  (public-ready, no-mixing design)
# - Compute A<->C co-occurrence probability, lift, NPMI
# - Outputs a single tagged CSV for downstream figures
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("dplyr","readr","stringr","tibble","purrr"))

# -------------------------
# 1) Load inputs (STRICT, tagged)
# -------------------------
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_hits_tag)) {
  stop("Tagged hits_matrix not found: ", f_hits_tag,
       "\nRun 02_build_hits_matrix.R with tagged saving enabled.")
}
df <- readr::read_csv(f_hits_tag, show_col_types = FALSE)
log_msg("05 loaded hits_matrix:", f_hits_tag, " n=", nrow(df))

dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
dic <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","regex","class") %in% names(dic))) stop("Dictionary must have term/regex/class columns.")
log_msg("05 using dictionary:", basename(dic_path), " n=", nrow(dic))

# _var mapping: exclude child terms (same policy)
var_childs <- unique(dic$term[grepl("_var$", dic$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) %>%
  dplyr::distinct(child, .keep_all = TRUE)

# Helpers
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))
hit_cols <- names(df)[startsWith(names(df), "hit__")]
terms_all <- sub("^hit__", "", hit_cols)

get_hit <- function(t) {
  nm <- paste0("hit__", t)
  if (nm %in% names(df)) as.integer(ifelse(is.na(df[[nm]]), 0L, df[[nm]])) else integer(nrow(df))
}

# -------------------------
# 2) Define A and C sets
# -------------------------
C_terms <- dic$term[dic$class == "outcome"]
C_terms <- intersect(C_terms, terms_all)
stopifnot(length(C_terms) > 0)

# "A" here is "all non-outcome terms" (same as your original 7) block)
A_terms_for_npmi <- setdiff(terms_all, C_terms)
A_terms_for_npmi <- setdiff(A_terms_for_npmi, term_groups$child)  # exclude _var children

log_msg(sprintf("05 sets |A|=%d |C|=%d", length(A_terms_for_npmi), length(C_terms)))

# -------------------------
# 3) Probabilities and NPMI/Lift
# -------------------------
p_of <- function(v) mean(binv(v))

# Precompute marginals
pA <- sapply(A_terms_for_npmi, function(t) p_of(get_hit(t)))
pC <- sapply(C_terms,          function(t) p_of(get_hit(t)))

# Compute all pairs (A x C)
co_npmi <- purrr::map_dfr(A_terms_for_npmi, function(a){
  vA <- get_hit(a)
  purrr::map_dfr(C_terms, function(c){
    vC <- get_hit(c)
    pAc <- mean(binv(vA) & binv(vC))

    pA_ <- pA[[a]]
    pC_ <- pC[[c]]

    # small epsilon for numerical stability
    eps <- 1e-12
    pmi  <- log((pAc + eps) / (pA_ * pC_ + eps))
    npmi <- if (pAc > 0) pmi / (-log(pAc + eps)) else 0

    tibble::tibble(
      A = a, C = c,
      pA = pA_, pC = pC_, pAC = pAc,
      lift = (pAc + eps) / (pA_ * pC_ + eps),
      npmi = npmi
    )
  })
})

# -------------------------
# 4) Output (tagged, fixed name)
# -------------------------
out_npmi <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(co_npmi, out_npmi)
log_msg("WROTE:", out_npmi, " n=", nrow(co_npmi))

log_msg("=== DONE 05_ac_npmi ===")
