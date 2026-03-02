# =========================================================
# 02_BuildHitsMatrix_Final.R  (based on your hardened code)
# - Loads articles_YYYYMMDD.csv
# - Loads dictionary: ra_ild_dictionary_analysis_v1_genetics.csv
# - Builds hit__ matrix with _var parent merge + mutual exclusivity fixes
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("dplyr","readr","stringr","tibble"))

STAMP_DATE <- format(Sys.Date(), "%Y%m%d")

# Load articles saved by 01 (use today's file if present, otherwise pick the latest).
f_today <- file.path(DIR_PROC, sprintf("articles_%s.csv", STAMP_DATE))

if (file.exists(f_today)) {
  f_articles_csv <- f_today
} else {
  cand <- Sys.glob(file.path(DIR_PROC, "articles_*.csv"))
  if (length(cand) == 0) stop("No articles_*.csv found in ", DIR_PROC, " . Run 01_fetch_pubmed.R first.")
  # Assumption: lexicographic sort makes the largest YYYYMMDD the latest.
  f_articles_csv <- sort(cand, decreasing = TRUE)[1]
}

articles <- readr::read_csv(f_articles_csv, show_col_types = FALSE)
log_msg("02 loaded articles:", f_articles_csv, " n=", nrow(articles))


# ------------------ Load dictionary (fixed filename). ------------------
dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary file not found: ", dic_path)

dic <- readr::read_csv(dic_path, show_col_types = FALSE)
log_msg("Using dictionary:", basename(dic_path), "(n=", nrow(dic), ")")
# Log dictionary size here each run (expected n~159).
log_msg("DIC head:", paste(head(dic$term, 6), collapse=" / "))

# Prepare _var alias -> parent-term aggregation.
var_childs <- unique(dic$term[grepl("_var$", dic$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) %>%
  dplyr::distinct(child, .keep_all = TRUE)

# ------------------ Hits matrix ------------------------
articles <- articles |> dplyr::mutate(text = ifelse(is.na(text), "", text))

match_one <- function(pattern, texts) {
  ok <- !is.na(texts); out <- integer(length(texts))
  if (any(ok)) out[ok] <- as.integer(grepl(pattern, texts[ok], ignore.case = FALSE, perl = TRUE))
  out
}

for (i in seq_len(nrow(dic))) {
  nm <- paste0("hit__", dic$term[i])
  if (!nm %in% names(articles)) {
    articles[[nm]] <- match_one(dic$regex[i], articles$text)
  }
}

# Parent aggregation for _var terms.
if (nrow(term_groups)) {
  for (k in seq_len(nrow(term_groups))){
    ch <- term_groups$child[k]; pa <- term_groups$parent[k]
    ch_col <- paste0("hit__", ch)
    pa_col <- paste0("hit__", pa)
    if (ch_col %in% names(articles)) {
      if (!pa_col %in% names(articles)) {
        articles[[pa_col]] <- articles[[ch_col]]
      } else {
        articles[[pa_col]] <- as.integer((articles[[pa_col]]==1L) | (articles[[ch_col]]==1L))
      }
    }
  }
}

# Mutual exclusivity / consistency fixes: ACPA/CRP/ESR.
safe_col <- function(nm) {
  if (nm %in% names(articles)) as.integer(articles[[nm]])
  else integer(nrow(articles))
}
z <- function(nm) safe_col(nm)

if (all(c("hit__ACPA_neg","hit__ACPA_pos") %in% names(articles))) {
  idx_neg <- which(z("hit__ACPA_neg")==1L)
  if (length(idx_neg)) {
    articles$hit__ACPA_pos[idx_neg]  <- 0L
    if ("hit__ACPA_RF_high" %in% names(articles)) articles$hit__ACPA_RF_high[idx_neg] <- 0L
  }
}
if (all(c("hit__CRP_high","hit__CRP_low") %in% names(articles))) {
  both <- which(z("hit__CRP_high")==1L & z("hit__CRP_low")==1L)
  if (length(both)) articles$hit__CRP_low[both] <- 0L
}
if (all(c("hit__ESR_high","hit__ESR_low") %in% names(articles))) {
  both <- which(z("hit__ESR_high")==1L & z("hit__ESR_low")==1L)
  if (length(both)) articles$hit__ESR_low[both] <- 0L
}

# Save (backward compatible filenames).
f_hits <- file.path(DIR_PROC, sprintf("hits_matrix_%s.csv", STAMP_DATE))
hits_df <- articles |> dplyr::select(pmid, year, journal, pubtype, dplyr::starts_with("hit__"))
readr::write_csv(hits_df, f_hits)
log_msg("Saved hits:", basename(f_hits))

# Also save tag-stamped outputs (prevents mixed dictionary/corpus).
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(hits_df, f_hits_tag)
log_msg("Saved hits (tagged):", basename(f_hits_tag))

log_msg("=== DONE 02_build_hits_matrix ===")
