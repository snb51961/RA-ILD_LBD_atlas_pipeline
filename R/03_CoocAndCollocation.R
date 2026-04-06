# =========================================================
# 03_CoocAndCollocation_Final.R
# RA-ILD pipeline (public-ready)
# - Co-occurrence: each class vs outcomes, and biomarker vs outcomes
# - Collocations: bigram/trigram from corpus text
# - Strict I/O tagging to avoid mixing corpora/dictionaries
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("dplyr","readr","stringr","tibble","purrr","quanteda","quanteda.textstats"))

# -------------------------
# 0) Helpers
# -------------------------
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

# Prefer tagged inputs first (prevents corpus/dictionary mixing).
# Ideally 02 saves tagged outputs; otherwise fall back to legacy hits_matrix_YYYYMMDD.
find_hits_matrix <- function(){
  # tagged path (recommended)
  f_tagged <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  if (file.exists(f_tagged)) return(f_tagged)

  # fallback: today's legacy
  f_legacy_today <- file.path(DIR_PROC, sprintf("hits_matrix_%s.csv", format(Sys.Date(), "%Y%m%d")))
  if (file.exists(f_legacy_today)) return(f_legacy_today)

  # fallback: latest legacy (only within data_proc; still safer than mixing output dirs)
  cand <- Sys.glob(file.path(DIR_PROC, "hits_matrix_*.csv"))
  if (!length(cand)) stop("hits_matrix file not found. Run 02_build_hits_matrix.R first.")
  cand[which.max(file.info(cand)$mtime)]
}

# Dictionary loader (one true source)
load_dic <- function(){
  dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
  if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
  dic <- readr::read_csv(dic_path, show_col_types = FALSE)
  # minimal checks
  if (!all(c("term","regex","class") %in% names(dic))) {
    stop("Dictionary must contain columns: term, regex, class")
  }
  if (any(duplicated(dic$term))) stop("Duplicated term in dictionary.")
  log_msg("03 using dictionary:", basename(dic_path), " n=", nrow(dic))
  log_msg("DIC head:", paste(head(dic$term, 6), collapse=" / "))
  dic
}

# "_var" child -> parent mapping (same as your logic)
make_term_groups <- function(dic){
  var_childs <- unique(dic$term[grepl("_var$", dic$term)])
  tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) %>%
    dplyr::distinct(child, .keep_all = TRUE)
}

# safe select matrix for a term set
safe_hit_matrix <- function(df, terms){
  cols <- paste0("hit__", terms)
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(NULL)
  as.matrix(df[, cols, drop = FALSE])
}

# -------------------------
# 1) Load inputs
# -------------------------
f_hits <- find_hits_matrix()
df <- readr::read_csv(f_hits, show_col_types = FALSE)
log_msg("03 loaded hits_matrix:", f_hits, " n=", nrow(df))

# Ensure required columns exist
need_cols <- c("year","pmid")
miss <- setdiff(need_cols, names(df))
if (length(miss)) stop("hits_matrix missing required columns: ", paste(miss, collapse=", "))

# text field is not in hits_matrix by default; collocations need articles text.
# We'll look for articles_{today}.csv (legacy) or a tagged articles file if you later add it.
find_articles <- function(){
  # tagged optional
  f_tagged <- file.path(DIR_TABLE, sprintf("articles_%s__%s.csv", CORPUS_TAG, format(Sys.Date(), "%Y%m%d")))
  if (file.exists(f_tagged)) return(f_tagged)

  f_legacy_today <- file.path(DIR_PROC, sprintf("articles_%s.csv", format(Sys.Date(), "%Y%m%d")))
  if (file.exists(f_legacy_today)) return(f_legacy_today)

  cand <- Sys.glob(file.path(DIR_PROC, "articles_*.csv"))
  if (!length(cand)) return(NA_character_)
  cand[which.max(file.info(cand)$mtime)]
}

f_articles <- find_articles()
articles <- if (!is.na(f_articles) && file.exists(f_articles)) {
  readr::read_csv(f_articles, show_col_types = FALSE)
} else {
  tibble::tibble()
}
if (nrow(articles)) {
  log_msg("03 loaded articles:", f_articles, " n=", nrow(articles))
} else {
  log_msg("03 articles not found -> collocations will be skipped.")
}

dic <- load_dic()
term_groups <- make_term_groups(dic)

get_terms <- function(cls) dic |> dplyr::filter(class == cls) |> dplyr::pull(term)

# Outcomes must exist
ae_terms <- dic |> dplyr::filter(class == "outcome") |> dplyr::pull(term)
if (!length(ae_terms)) stop("No outcome terms found in dictionary.")
ae_main <- if ("AE-ILD" %in% ae_terms) "AE-ILD" else ae_terms[1]
log_msg("03 outcomes:", paste(ae_terms, collapse=", "), " | AE_main=", ae_main)

# -------------------------
# 2) Co-occurrence: class vs outcomes
# -------------------------
build_for_class <- function(cls){
  terms <- get_terms(cls)
  if (!length(terms)) return(invisible(NULL))

  # exclude _var children (same as your policy)
  terms <- setdiff(terms, term_groups$child)

  X_g  <- safe_hit_matrix(df, terms)
  X_ae <- safe_hit_matrix(df, ae_terms)
  if (is.null(X_g) || is.null(X_ae)) return(invisible(NULL))

  co_mat <- t(X_g) %*% X_ae
  co_df <- as.data.frame(as.table(co_mat)) |>
    dplyr::rename(term = Var1, outcome = Var2, count = Freq) |>
    dplyr::mutate(
      term = stringr::str_remove(term, "^hit__"),
      outcome = stringr::str_remove(outcome, "^hit__")
    ) |>
    dplyr::arrange(dplyr::desc(count))

  # Tagged output (no timestamp guessing)
  out_name <- sprintf("cooc_%s_outcomes_%s__%s.csv", cls, CORPUS_TAG, DIC_TAG)
  out_path <- file.path(DIR_TABLE, out_name)
  readr::write_csv(co_df, out_path)
  log_msg("WROTE:", out_path)
  invisible(co_df)
}

# Run for all classes you used (keep identical order)
CLASSES <- c("bio","drug","gene","exposure","host","event","nonpharm","vaccine",
             "pattern","airway","qct","biomarker","pft","activity","phenotype",
             "population","system","radiology","cell","molecular","microbiome","complication","trend")

invisible(purrr::map(CLASSES, build_for_class))

# Biomarker-like special (your old cooc_bio_AE)
bio_terms <- get_terms("bio")
if (length(bio_terms)) {
  bio_terms <- setdiff(bio_terms, term_groups$child)
  X_bio <- safe_hit_matrix(df, bio_terms)
  X_ae  <- safe_hit_matrix(df, ae_terms)
  if (!is.null(X_bio) && !is.null(X_ae)) {
    co_mat <- t(X_bio) %*% X_ae
    co_df <- as.data.frame(as.table(co_mat)) |>
      dplyr::rename(bio = Var1, ae = Var2, count = Freq) |>
      dplyr::arrange(dplyr::desc(count))
    out_path <- file.path(DIR_TABLE, sprintf("cooc_bio_AE_%s__%s.csv", CORPUS_TAG, DIC_TAG))
    readr::write_csv(co_df, out_path)
    log_msg("WROTE:", out_path)
  }
}

# -------------------------
# 3) Collocations (token-based) - optional
# -------------------------
if (nrow(articles) && "text" %in% names(articles)) {
  # Domain stopwords: same philosophy as your original
  sw_general <- quanteda::stopwords("en")
  sw_domain  <- c("rheumatoid","arthritis","interstitial","lung","disease","ild","ip","pulmonary","fibrosis",
                  "patients","patient","study","studies","review","case","cases","report","reports",
                  "introduction","conclusion","methods","background","objective","aim","result","results","purpose",
                  "mg","ml","day","days","week","weeks","year","years")

  corp <- quanteda::corpus(articles, text_field = "text")
  toks <- corp |>
    quanteda::tokens(remove_punct=TRUE, remove_numbers=TRUE, remove_symbols=TRUE) |>
    quanteda::tokens_tolower() |>
    quanteda::tokens_remove(c(sw_general, sw_domain))

  coll2 <- quanteda.textstats::textstat_collocations(toks, size=2, min_count=10, smoothing=0.5)
  coll3 <- quanteda.textstats::textstat_collocations(toks, size=3, min_count=8,  smoothing=0.5)

  clean_noise <- function(D){
    if (is.null(D) || nrow(D) == 0) return(D)
    D |>
      dplyr::mutate(phrase = stringr::str_replace_all(collocation, "_", " ")) |>
      dplyr::filter(!grepl("\\d", phrase), nchar(phrase) >= 5, !grepl("(^[a-z]$|^[-_]+$)", phrase))
  }

  coll2c <- clean_noise(coll2)
  coll3c <- clean_noise(coll3)

  f2 <- file.path(DIR_TABLE, sprintf("collocation_bigram_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  f3 <- file.path(DIR_TABLE, sprintf("collocation_trigram_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  readr::write_csv(coll2c |> dplyr::arrange(dplyr::desc(lambda)), f2)
  readr::write_csv(coll3c |> dplyr::arrange(dplyr::desc(lambda)), f3)
  log_msg("WROTE:", f2)
  log_msg("WROTE:", f3)
} else {
  log_msg("03 collocations skipped (articles/text not available).")
}

log_msg("=== DONE 03_cooc_and_collocation ===")