
# =========================================================
# 17_TimeSlice_BuildAndAnalyze.R
# Temporal holdout build script for the RA-ILD atlas project
#
# Fixed inputs:
#   temporal_holdout/input/articles_main_original_frozen.csv
#   temporal_holdout/input/ra_ild_dictionary_analysis_v1_genetics.csv
#
# Fixed outputs:
#   temporal_holdout/output/full/
#   temporal_holdout/output/discovery/
#   temporal_holdout/output/holdout/
#   temporal_holdout/output/temporal_holdout_manifest.csv
#
# Notes
#   - Self-contained: does NOT source 00_setup or any of 02/04/05/06/13.
#   - Main manuscript outputs remain untouched.
#   - Full slice is produced for internal pipeline-reproducibility checks.
# =========================================================

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","tibble","purrr","tidyr"))

# ---------------------------------------------------------
# Path resolution anchored to THIS script
# ---------------------------------------------------------
.this_file <- function(default_name) {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(getwd(), default_name), winslash = "/", mustWork = FALSE)
}
SCRIPT_FILE <- .this_file("17_TimeSlice_BuildAndAnalyze.R")
SCRIPT_DIR  <- dirname(SCRIPT_FILE)
TH_ROOT     <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = FALSE)
TH_INPUT    <- file.path(TH_ROOT, "input")
TH_OUTPUT   <- file.path(TH_ROOT, "output")
TH_LOG      <- file.path(TH_OUTPUT, "log")
dir.create(TH_OUTPUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TH_LOG, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(TH_LOG, paste0("17_build_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  message(msg)
  write(msg, file = log_file, append = TRUE)
}

# ---------------------------------------------------------
# Fixed inputs
# ---------------------------------------------------------
ARTICLES_FILE <- file.path(TH_INPUT, "articles_main_original_frozen.csv")
DIC_FILE      <- file.path(TH_INPUT, "ra_ild_dictionary_analysis_v1_genetics.csv")

if (!file.exists(ARTICLES_FILE)) stop("Missing input corpus CSV: ", ARTICLES_FILE)
if (!file.exists(DIC_FILE)) stop("Missing dictionary CSV: ", DIC_FILE)

# ---------------------------------------------------------
# Slice definitions
# ---------------------------------------------------------
FULL_START      <- 1980L
FULL_END        <- 2025L
DISCOVERY_START <- 1980L
DISCOVERY_END   <- 2020L
HOLDOUT_START   <- 2021L
HOLDOUT_END     <- 2025L

OUTCOME_GROUPS <- c("AE-ILD", "progression", "mortality")

HOTSPOT_TOP_K  <- 12L
BRIDGE_TOP_K   <- 10L
SIGNED_MIN_ARTICLES <- 3L
BIOMARKER_CORE <- c("KL6","SP_D","ACPA_RF","MMP7","APRIL","Fibrosis_pathway","Cytokine_inflammation")


# ---------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------
# ABC ranking
LIFT_MIN_AB <- 1.20
NPMI_MIN_AB <- 0.05
LIFT_MIN_BC <- 1.20
NPMI_MIN_BC <- 0.05
NOVELTY_LAMBDA <- 1.0
A_CLASSES <- c("drug","bio","gene","exposure","host","event","nonpharm","vaccine","population","system","strategy")
B_CLASSES <- c("pattern","radiology","airway","qct","pft","activity","biomarker","molecular","microbiome","cell","complication","phenotype","trend")
TERM_PRIOR <- tibble::tibble(
  term  = c("Oxidative_stress","oxidative_stress","Oxidative_stress_var"),
  prior = c(0.2, 0.2, 0.2)
)

# Signed effects
CHAR_WINDOW  <- 180
MAX_CONTEXTS <- 3
WINDOW <- 1
A_CLASSES_SIGNED <- unique(c(A_CLASSES, "biomarker", "molecular", "microbiome"))
AE_EXPAND_RE <- "(?i)((acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation)|acute[-\\s]*on[-\\s]*chronic).{0,80}(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b)|(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b).{0,80}acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation))"
NEG_RE <- "(?i)\\b(no|not|without|neither|never|absence\\s+of|lack\\s+of)\\b"
HEDGE_RE <- "(?i)\\b(may|might|tend\\s+to|trend(s)?\\s+toward|borderline|possibly|suggests?)\\b"
SERO_NEG_RE <- "(?i)\\b(seronegative|ACPA\\s*negative|anti-?CCP\\s*negative|RF\\s*negative)\\b"
rx <- list(
  inc = "(?i)(associated\\s+with\\s+(an?\\s+)?increase(d)?|increase(d|s)?\\s+(in|of)|higher\\s+risk|elevat(ed|es)|worsen(ed|ing)|exacerbat(ed|es)|risk\\s+factor|predict(or|ive)\\s+of|trigger(s|ed)?)",
  dec = "(?i)(associated\\s+with\\s+(a\\s+)?decrease(d)?|decrease(d|s)?\\s+(in|of)|lower\\s+risk|reduc(ed|es)|protective|improv(ed|es|ement)|ameliorat(ed|es)|attenuat(ed|es)|inhibit(ed|s))",
  null= "(?i)(not\\s+associated|no\\s+(significant\\s+)?association|no\\s+significant\\s+(difference|effect)|nonsignificant|non[-\\s]?significant|did\\s+not\\s+(find|show)|failed\\s+to\\s+(show|demonstrate)|similar\\s+(rates?|risk)|comparable)",
  meas= "(?i)\\b((?:a\\s*OR|aOR|OR|HR|aHR|RR|IRR|SHR|adj(?:usted)?\\s*(?:odds|hazard|risk)\\s*ratio))\\b\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)",
  ci  = "(?i)95%\\s*CI\\s*[:=]?\\s*\\(?\\s*([0-9]+\\.?[0-9]*)\\s*[-–,]\\s*([0-9]+\\.?[0-9]*)\\s*\\)?",
  p   = "(?i)\\bp\\s*[<=>]\\s*0\\.?0*([0-9]\\d?)"
)

# Biomarker atlas
B_CLASSES_ATLAS <- c("biomarker","molecular","cell","microbiome")
META_TERMS <- c("Cytokine_inflammation","Fibrosis_pathway","Oxidative_stress_pathway","Immune_checkpoint")
MUST_KEEP_BIOMARKERS <- c("KL6","SP_D","ACPA_RF","CRP_ESR","MMP7","MPO","NLR","CAR","CXCL9","IL17A","APRIL")

# ---------------------------------------------------------
# Load inputs
# ---------------------------------------------------------
ART <- readr::read_csv(ARTICLES_FILE, show_col_types = FALSE)
DIC <- readr::read_csv(DIC_FILE, show_col_types = FALSE)

need_art <- c("pmid","title","abstract","year","pubtype")
need_dic <- c("term","regex","class")
miss_art <- setdiff(need_art, names(ART))
miss_dic <- setdiff(need_dic, names(DIC))
if (length(miss_art)) stop("Corpus CSV missing columns: ", paste(miss_art, collapse = ", "))
if (length(miss_dic)) stop("Dictionary CSV missing columns: ", paste(miss_dic, collapse = ", "))

ART <- ART |>
  dplyr::mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    title = as.character(title),
    abstract = as.character(abstract),
    pubtype = as.character(pubtype)
  )
if (!("text" %in% names(ART))) {
  ART$text <- paste(ART$title, ART$abstract, sep = " ")
} else {
  ART$text <- ifelse(is.na(ART$text), paste(ART$title, ART$abstract, sep = " "), ART$text)
}

DIC <- DIC |>
  dplyr::mutate(
    term  = as.character(term),
    regex = as.character(regex),
    class = tolower(trimws(as.character(class)))
  )

if (anyDuplicated(DIC$term)) stop("Duplicated dictionary terms detected. Please fix the dictionary input.")

var_childs <- unique(DIC$term[grepl("_var$", DIC$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) |>
  dplyr::distinct(child, .keep_all = TRUE)

canon_outcome <- function(x){
  xl <- stringr::str_to_lower(stringr::str_replace_all(x, "–", "-"))
  dplyr::case_when(
    stringr::str_detect(xl, "ae-ild|aeild|acute exacerb|exacerb") ~ "AE-ILD",
    stringr::str_detect(xl, "progress|decline|fvc|dlco|worsen")   ~ "progression",
    stringr::str_detect(xl, "mortality|death|surviv")             ~ "mortality",
    TRUE ~ NA_character_
  )
}

OUTCOME_DIC <- DIC |>
  dplyr::filter(class == "outcome") |>
  dplyr::mutate(group = canon_outcome(term)) |>
  dplyr::filter(!is.na(group), group %in% OUTCOME_GROUPS)

if (nrow(OUTCOME_DIC) == 0) stop("No outcome terms could be mapped to canonical outcome groups.")

OUTCOME_TERM_MAP <- split(OUTCOME_DIC$term, OUTCOME_DIC$group)
OUTCOME_REGEX_MAP <- lapply(split(OUTCOME_DIC$regex, OUTCOME_DIC$group), function(xs){
  xs <- unique(xs[!is.na(xs) & nzchar(xs)])
  if (!length(xs)) return(NA_character_)
  paste0("(", paste(xs, collapse = "|"), ")")
})

log_msg("ARTICLES_FILE:", ARTICLES_FILE)
log_msg("DIC_FILE     :", DIC_FILE)
log_msg("n_articles(all) =", nrow(ART))
log_msg("year range =", paste(range(ART$year, na.rm = TRUE), collapse = " - "))
log_msg("Outcome groups available:", paste(names(OUTCOME_TERM_MAP), collapse = ", "))

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))
p_of <- function(v) mean(binv(v))
safelog <- function(z){ x <- suppressWarnings(log(z)); x[!is.finite(x)] <- 0; x }
prior_of <- function(t){ v <- TERM_PRIOR$prior[match(t, TERM_PRIOR$term)]; ifelse(is.na(v), 1.0, v) }
safe_parent <- function(x) ifelse(grepl("_var$", x), sub("_var$", "", x), x)

pair_stats <- function(x, y){
  x <- binv(x); y <- binv(y)
  n11 <- sum(x==1 & y==1); n10 <- sum(x==1 & y==0); n01 <- sum(x==0 & y==1); n00 <- sum(x==0 & y==0)
  N <- n11 + n10 + n01 + n00
  if (N == 0) return(list(or = NA_real_, p = 1, n11 = 0L, lift = 0, npmi = 0))
  OR <- ((n11+0.5)*(n00+0.5))/((n10+0.5)*(n01+0.5))
  p  <- tryCatch(suppressWarnings(fisher.test(matrix(c(n11,n10,n01,n00),2,byrow = TRUE))$p.value), error = function(e) 1)
  pA <- (n11+n10)/N; pY <- (n11+n01)/N; pAc <- n11/N
  pmi <- if (pAc > 0 && pA > 0 && pY > 0) log(pAc/(pA*pY)) else 0
  npmi <- if (pAc > 0) pmi/(-log(pAc)) else 0
  lift <- if (pA > 0 && pY > 0) pAc/(pA*pY) else 0
  list(or = as.numeric(OR), p = p, n11 = as.integer(n11), lift = lift, npmi = npmi)
}

split_sentences <- function(txt){
  if (is.na(txt) || !nzchar(txt)) return(character(0))
  t <- gsub("[\\r\\n]+", " ", txt)
  t <- gsub("([A-Za-z])\\.(\\s*[A-Za-z])", "\\1. \\2", t)
  xs <- unlist(strsplit(t, "(?<=[\\.!?;:])\\s+|\\s+(?=Results?:|Conclusions?:|Background:)", perl = TRUE))
  xs <- trimws(xs)
  xs[nzchar(xs)]
}

all_spans <- function(pattern, text){
  if (is.na(text) || !nzchar(text)) return(list(starts = integer(0), ends = integer(0)))
  m <- gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  if (length(m) == 0 || m[[1]][1] == -1) return(list(starts = integer(0), ends = integer(0)))
  starts <- as.integer(m[[1]])
  lens <- attr(m[[1]], "match.length")
  list(starts = starts, ends = starts + lens - 1L)
}

span_contexts <- function(text, reA, reC, char_window = CHAR_WINDOW, max_ctx = MAX_CONTEXTS){
  if (is.na(text) || !nzchar(text)) return(character(0))
  sa <- all_spans(reA, text); sc <- all_spans(reC, text)
  if (!length(sa$starts) || !length(sc$starts)) return(character(0))
  comb <- expand.grid(i = seq_along(sa$starts), j = seq_along(sc$starts))
  comb$dist <- mapply(function(i,j){
    a_mid <- (sa$starts[i] + sa$ends[i]) / 2
    c_mid <- (sc$starts[j] + sc$ends[j]) / 2
    abs(a_mid - c_mid)
  }, comb$i, comb$j)
  comb <- comb[order(comb$dist), , drop = FALSE]
  ctxs <- character(0)
  used <- matrix(FALSE, nrow = length(sa$starts), ncol = length(sc$starts))
  for (k in seq_len(nrow(comb))){
    i <- comb$i[k]; j <- comb$j[k]
    if (used[i,j]) next
    used[i,] <- TRUE; used[,j] <- TRUE
    mid1 <- min(sa$starts[i], sc$starts[j])
    mid2 <- max(sa$ends[i],   sc$ends[j])
    bgn  <- max(1L, mid1 - char_window)
    end  <- min(nchar(text), mid2 + char_window)
    ctxs <- c(ctxs, substr(text, bgn, end))
    if (length(ctxs) >= max_ctx) break
  }
  unique(ctxs)
}

des_w <- function(pubtype, year){
  year <- suppressWarnings(as.integer(year))
  w <- 1
  if (!is.na(pubtype) && grepl("Randomized|Controlled Clinical Trial", pubtype, ignore.case = TRUE)) w <- w + 0.8
  if (!is.na(pubtype) && grepl("Cohort|Case-Control", pubtype, ignore.case = TRUE)) w <- w + 0.3
  if (!is.na(pubtype) && grepl("Case Reports?", pubtype, ignore.case = TRUE)) w <- w - 0.3
  if (!is.na(year)) w <- w + 0.1 * pmax(0, year - 2015) / 10
  pmax(w, 0.2)
}

infer_one <- function(s){
  null_flag <- grepl(rx$null, s, perl = TRUE)
  inc_flag  <- grepl(rx$inc,  s, perl = TRUE)
  dec_flag  <- grepl(rx$dec,  s, perl = TRUE)
  neg_flag  <- grepl(NEG_RE,  s, perl = TRUE)

  if (neg_flag && (inc_flag || dec_flag)) {
    inc_flag <- FALSE
    dec_flag <- FALSE
    null_flag <- TRUE
  }

  hedge_flag <- grepl(HEDGE_RE, s, perl = TRUE)
  m <- stringr::str_match(s, rx$meas)
  val <- suppressWarnings(as.numeric(if (is.null(m)) NA else m[,3]))
  ci  <- stringr::str_match(s, rx$ci)
  lo  <- suppressWarnings(as.numeric(ci[,2]))
  hi  <- suppressWarnings(as.numeric(ci[,3]))
  p_sig <- grepl(rx$p, s, perl = TRUE)

  sig_by_ci <- if (!is.na(lo) && !is.na(hi)) if (lo > 1 || hi < 1) TRUE else FALSE else NA
  sig <- if (!is.na(sig_by_ci)) sig_by_ci else if (p_sig) TRUE else NA

  sgn <- NA_real_
  weight <- 1
  if (!is.na(val)) sgn <- if (val > 1) +1 else if (val < 1) -1 else 0
  if (is.na(sgn)) sgn <- if (null_flag) 0 else if (inc_flag && !dec_flag) +1 else if (dec_flag && !inc_flag) -1 else NA

  if (!is.na(sig) && sig) weight <- weight + 1
  if (grepl("(?i)(adjusted|multivaria(te|ted))", s, perl = TRUE)) weight <- weight + 0.5
  if (hedge_flag) weight <- weight * 0.6
  if (is.na(sgn) && grepl(SERO_NEG_RE, s, perl = TRUE)) { sgn <- -1; weight <- weight * 0.8 }

  list(sign = sgn, weight = weight)
}

wilson_balance_ci <- function(pos, neg){
  n <- pos + neg
  if (is.na(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  p <- pos / n
  z <- 1.96
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  half   <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / denom
  lo <- center - half
  hi <- center + half
  c(low = 2*lo - 1, high = 2*hi - 1)
}

grp_col_name <- function(g) paste0("grp__", gsub("[^A-Za-z0-9]+", "_", g))
get_group_hit <- function(df, group){
  nm <- grp_col_name(group)
  if (nm %in% names(df)) binv(df[[nm]]) else integer(nrow(df))
}

# ---------------------------------------------------------
# Core builders
# ---------------------------------------------------------
build_hits_matrix <- function(articles_df, dic, outcome_term_map){
  articles_df <- articles_df |> dplyr::mutate(text = ifelse(is.na(text), "", text))

  match_one <- function(pattern, texts) {
    ok <- !is.na(texts)
    out <- integer(length(texts))
    if (any(ok)) out[ok] <- as.integer(grepl(pattern, texts[ok], ignore.case = FALSE, perl = TRUE))
    out
  }

  out <- articles_df
  for (i in seq_len(nrow(dic))) {
    nm <- paste0("hit__", dic$term[i])
    if (!nm %in% names(out)) out[[nm]] <- match_one(dic$regex[i], out$text)
  }

  # Parent aggregation for _var terms
  if (nrow(term_groups)) {
    for (k in seq_len(nrow(term_groups))) {
      ch <- term_groups$child[k]
      pa <- term_groups$parent[k]
      ch_col <- paste0("hit__", ch)
      pa_col <- paste0("hit__", pa)
      if (ch_col %in% names(out)) {
        if (!pa_col %in% names(out)) out[[pa_col]] <- out[[ch_col]]
        else out[[pa_col]] <- as.integer((out[[pa_col]] == 1L) | (out[[ch_col]] == 1L))
      }
    }
  }

  # Canonical grouped outcomes
  for (g in names(outcome_term_map)) {
    hit_cols <- paste0("hit__", outcome_term_map[[g]])
    hit_cols <- hit_cols[hit_cols %in% names(out)]
    nm <- grp_col_name(g)
    if (length(hit_cols)) out[[nm]] <- as.integer(rowSums(as.data.frame(out[hit_cols])) > 0)
    else out[[nm]] <- integer(nrow(out))
  }

  out
}

compute_ac_npmi_grouped <- function(hm, dic){
  hit_cols <- names(hm)[startsWith(names(hm), "hit__")]
  terms_all <- sub("^hit__", "", hit_cols)
  raw_outcomes <- intersect(dic$term[dic$class == "outcome"], terms_all)
  A_terms <- setdiff(terms_all, raw_outcomes)
  A_terms <- setdiff(A_terms, term_groups$child)

  get_hit <- function(t) {
    nm <- paste0("hit__", t)
    if (nm %in% names(hm)) binv(hm[[nm]]) else integer(nrow(hm))
  }

  pA <- sapply(A_terms, function(t) p_of(get_hit(t)))

  purrr::map_dfr(A_terms, function(a){
    vA <- get_hit(a)
    purrr::map_dfr(OUTCOME_GROUPS, function(cg){
      vC <- get_group_hit(hm, cg)
      pAc <- mean(binv(vA) & binv(vC))
      pA_ <- pA[[a]]
      pC_ <- p_of(vC)
      eps <- 1e-12
      pmi  <- log((pAc + eps) / (pA_ * pC_ + eps))
      npmi <- if (pAc > 0) pmi / (-log(pAc + eps)) else 0
      tibble::tibble(
        A = a,
        C = cg,
        pA = pA_,
        pC = pC_,
        pAC = pAc,
        lift = (pAc + eps) / (pA_ * pC_ + eps),
        npmi = npmi
      )
    })
  })
}

compute_abc_rankings_ae <- function(hm, dic){
  hit_cols <- names(hm)[startsWith(names(hm), "hit__")]
  terms_all <- sub("^hit__", "", hit_cols)

  A_set_full <- intersect(terms_all, dic$term[dic$class %in% A_CLASSES])
  B_set_full <- intersect(terms_all, dic$term[dic$class %in% B_CLASSES])
  A_set <- setdiff(A_set_full, term_groups$child)
  B_set <- setdiff(B_set_full, term_groups$child)

  get_hit <- function(t) {
    nm <- paste0("hit__", t)
    if (nm %in% names(hm)) binv(hm[[nm]]) else integer(nrow(hm))
  }

  c_vec <- get_group_hit(hm, "AE-ILD")
  N_docs <- nrow(hm)
  df_B <- sapply(B_set, function(b) sum(get_hit(b)))
  idf_B <- log((N_docs + 1) / (df_B + 1)) + 1
  idfB_of <- function(b) idf_B[[b]]

  res <- list()
  for (a in A_set) {
    a_vec <- get_hit(a)
    for (b in B_set) {
      b_vec <- get_hit(b)
      AB <- pair_stats(a_vec, b_vec)
      BC <- pair_stats(b_vec, c_vec)
      AC <- pair_stats(a_vec, c_vec)
      passAB <- (AB$lift >= LIFT_MIN_AB) && (AB$npmi >= NPMI_MIN_AB)
      passBC <- (BC$lift >= LIFT_MIN_BC) && (BC$npmi >= NPMI_MIN_BC)
      if (!(passAB && passBC)) next
      res[[length(res)+1]] <- data.frame(
        A = a, B = b, C = "AE-ILD",
        AB_or = AB$or, AB_p = AB$p, AB_n11 = AB$n11, AB_lift = AB$lift, AB_npmi = AB$npmi,
        BC_or = BC$or, BC_p = BC$p, BC_n11 = BC$n11, BC_lift = BC$lift, BC_npmi = BC$npmi,
        AC_or = AC$or, AC_p = AC$p, AC_n11 = AC$n11,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(res)) return(tibble::tibble())

  dplyr::bind_rows(res) |>
    dplyr::mutate(
      AB_log = safelog(AB_or),
      BC_log = safelog(BC_or),
      AC_log = safelog(AC_or),
      AB_q = p.adjust(AB_p, method = "BH"),
      BC_q = p.adjust(BC_p, method = "BH"),
      AC_q = p.adjust(AC_p, method = "BH"),
      AB_qsig = -log10(pmax(AB_q, 1e-300)),
      BC_qsig = -log10(pmax(BC_q, 1e-300)),
      AC_pen = pmax(AC_log, 0) * as.integer(-log10(pmax(AC_p, 1e-300)) >= 3),
      priorA = vapply(A, prior_of, numeric(1)),
      priorB = vapply(B, prior_of, numeric(1)),
      idfB = vapply(B, idfB_of, numeric(1)),
      score_q = (AB_log * AB_qsig * idfB * priorA * priorB + BC_log * BC_qsig * idfB * priorB) - (AC_pen * NOVELTY_LAMBDA)
    ) |>
    dplyr::arrange(dplyr::desc(score_q))
}

compute_signed_summary_grouped <- function(hm, dic, outcome_regex_map){
  terms_all <- sub("^hit__", "", names(hm)[startsWith(names(hm), "hit__")])
  A_terms <- intersect(terms_all, dic$term[dic$class %in% A_CLASSES_SIGNED])

  get_hit <- function(t) {
    nm <- paste0("hit__", t)
    if (nm %in% names(hm)) binv(hm[[nm]]) else integer(nrow(hm))
  }

  analyze_pair <- function(a_term, outcome_group){
    mask <- which(get_hit(a_term) == 1 & get_group_hit(hm, outcome_group) == 1)
    c_re <- outcome_regex_map[[outcome_group]]
    if (is.null(c_re) || is.na(c_re)) return(NULL)
    if (identical(outcome_group, "AE-ILD")) {
      c_extra <- grepl(AE_EXPAND_RE, hm$text, ignore.case = TRUE, perl = TRUE)
      mask <- which(get_hit(a_term) == 1 & (get_group_hit(hm, outcome_group) == 1 | c_extra))
      c_re <- paste0("(", c_re, "|", AE_EXPAND_RE, ")")
    }
    if (!length(mask)) return(NULL)

    a_re <- dic$regex[dic$term == a_term][1]
    if (is.na(a_re) || !nzchar(a_re)) return(NULL)

    out_art <- list()
    for (idx in mask) {
      abs_txt <- hm$text[idx]
      if (is.na(abs_txt) || !nzchar(abs_txt)) next

      ss <- split_sentences(abs_txt)
      contexts <- character(0)
      src_weight <- 1

      Sa <- which(grepl(a_re, ss, ignore.case = TRUE, perl = TRUE))
      Sc <- which(grepl(c_re, ss, ignore.case = TRUE, perl = TRUE))

      if (length(Sa) && length(Sc)) {
        idxs <- sort(unique(c(Sa, Sc)))
        for (k in idxs) {
          i1 <- max(1, k - WINDOW)
          i2 <- min(length(ss), k + WINDOW)
          seg <- paste(ss[i1:i2], collapse = " ")
          if (any(grepl(a_re, seg, ignore.case = TRUE, perl = TRUE)) &&
              any(grepl(c_re, seg, ignore.case = TRUE, perl = TRUE))) {
            contexts <- c(contexts, seg)
          }
        }
      }

      if (!length(contexts)) {
        ctx_span <- span_contexts(abs_txt, a_re, c_re, CHAR_WINDOW, MAX_CONTEXTS)
        if (length(ctx_span)) {
          contexts <- c(contexts, ctx_span)
          src_weight <- 0.8
        }
      }

      if (!length(contexts) &&
          grepl(a_re, abs_txt, ignore.case = TRUE, perl = TRUE) &&
          grepl(c_re, abs_txt, ignore.case = TRUE, perl = TRUE)) {
        contexts <- c(contexts, abs_txt)
        src_weight <- 0.4
      }

      if (!length(contexts)) next

      art_w <- des_w(hm$pubtype[idx], hm$year[idx])
      scores <- vapply(unique(contexts), function(s) {
        r <- infer_one(s)
        ifelse(is.na(r$sign), 0, r$sign * r$weight * src_weight * art_w)
      }, numeric(1))

      sc <- sum(scores, na.rm = TRUE)
      lbl <- if (sc > 0) "risk_up" else if (sc < 0) "risk_down" else "no_effect_or_mixed"
      out_art[[length(out_art)+1]] <- data.frame(
        pmid = hm$pmid[idx],
        A = a_term,
        C = outcome_group,
        score = sc,
        label = lbl,
        stringsAsFactors = FALSE
      )
    }

    if (length(out_art)) dplyr::bind_rows(out_art) else NULL
  }

  art_res <- list()
  for (a in unique(A_terms)) {
    for (cg in OUTCOME_GROUPS) {
      z <- analyze_pair(a, cg)
      if (!is.null(z)) art_res[[length(art_res)+1]] <- z
    }
  }

  signed_art <- if (length(art_res)) dplyr::bind_rows(art_res) else tibble::tibble()
  if (!nrow(signed_art)) {
    return(tibble::tibble(
      A = character(), C = character(),
      articles = integer(), pos_articles = integer(), neg_articles = integer(),
      null_or_mix = integer(), net_score = double(), pos_ratio = double(),
      balance = double(), balance_low = double(), balance_high = double()
    ))
  }

  sum_tab <- signed_art |>
    dplyr::group_by(A, C) |>
    dplyr::summarise(
      articles = dplyr::n(),
      pos_articles = sum(label == "risk_up", na.rm = TRUE),
      neg_articles = sum(label == "risk_down", na.rm = TRUE),
      null_or_mix  = sum(label == "no_effect_or_mixed", na.rm = TRUE),
      net_score = sum(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pos_ratio = pos_articles / pmax(1, pos_articles + neg_articles),
      balance = (pos_articles - neg_articles) / pmax(1, pos_articles + neg_articles)
    )

  ci_mat <- t(mapply(wilson_balance_ci, sum_tab$pos_articles, sum_tab$neg_articles))
  sum_tab$balance_low  <- ci_mat[, "low"]
  sum_tab$balance_high <- ci_mat[, "high"]

  sum_tab |>
    dplyr::arrange(C, dplyr::desc(articles), dplyr::desc(abs(balance)))
}

compute_biomarker_atlas_grouped <- function(signed_summary, dic){
  if (!nrow(signed_summary)) {
    return(tibble::tibble(
      A = character(), C = character(),
      articles = double(), balance = double(),
      evidence = double(), score = double()
    ))
  }

  BIO_TERMS <- unique(c(
    as.character(dic$term[dic$class %in% B_CLASSES_ATLAS]),
    META_TERMS,
    MUST_KEEP_BIOMARKERS
  ))

  signed_summary |>
    dplyr::mutate(
      A_parent = safe_parent(A),
      C = canon_outcome(C),
      articles = as.numeric(articles),
      pos_articles = as.numeric(pos_articles),
      neg_articles = as.numeric(neg_articles),
      balance = as.numeric(balance),
      evidence = log10(pmax(pos_articles + neg_articles, 1)),
      score = (abs(balance) + 0.20) * (evidence + 0.50)
    ) |>
    dplyr::filter(!is.na(C), C %in% OUTCOME_GROUPS) |>
    dplyr::filter(grepl("_family$", A_parent) | (A_parent %in% BIO_TERMS)) |>
    dplyr::transmute(A = A_parent, C, articles, balance, evidence, score) |>
    dplyr::distinct()
}

# ---------------------------------------------------------
# Slice runner
# ---------------------------------------------------------
slices <- tibble::tribble(
  ~slice, ~year_min, ~year_max,
  "full", FULL_START, FULL_END,
  "discovery", DISCOVERY_START, DISCOVERY_END,
  "holdout", HOLDOUT_START, HOLDOUT_END
)

run_one_slice <- function(slice_name, year_min, year_max){
  log_msg("=== Slice:", slice_name, "(", year_min, "-", year_max, ") ===")

  art_slice <- ART |>
    dplyr::filter(!is.na(year), year >= year_min, year <= year_max)

  n_articles <- nrow(art_slice)
  log_msg("n_articles =", n_articles)
  if (n_articles == 0) stop("No articles in slice: ", slice_name)

  hm   <- build_hits_matrix(art_slice, DIC, OUTCOME_TERM_MAP)
  npmi <- compute_ac_npmi_grouped(hm, DIC)
  abc  <- compute_abc_rankings_ae(hm, DIC)
  sign <- compute_signed_summary_grouped(hm, DIC, OUTCOME_REGEX_MAP)
  biom <- compute_biomarker_atlas_grouped(sign, DIC)

  slice_dir <- file.path(TH_OUTPUT, slice_name)
  dir.create(slice_dir, recursive = TRUE, showWarnings = FALSE)

  f_hm   <- file.path(slice_dir, sprintf("hits_matrix_%s_%d_%d.csv", slice_name, year_min, year_max))
  f_npmi <- file.path(slice_dir, sprintf("coherence_grouped_%s_%d_%d.csv", slice_name, year_min, year_max))
  f_abc  <- file.path(slice_dir, sprintf("abc_rankings_ae_%s_%d_%d.csv", slice_name, year_min, year_max))
  f_sign <- file.path(slice_dir, sprintf("signed_effects_summary_grouped_%s_%d_%d.csv", slice_name, year_min, year_max))
  f_biom <- file.path(slice_dir, sprintf("biomarker_atlas_grouped_%s_%d_%d.csv", slice_name, year_min, year_max))

  readr::write_csv(hm |> dplyr::select(-dplyr::any_of(c("title","abstract","text"))), f_hm)
  readr::write_csv(npmi, f_npmi)
  readr::write_csv(abc, f_abc)
  readr::write_csv(sign, f_sign)
  readr::write_csv(biom, f_biom)

  tibble::tibble(
    slice = slice_name,
    year_min = year_min,
    year_max = year_max,
    n_articles = n_articles,
    hits_path = f_hm,
    npmi_path = f_npmi,
    abc_path  = f_abc,
    sign_path = f_sign,
    biom_path = f_biom
  )
}

manifest <- purrr::pmap_dfr(slices, run_one_slice)
manifest_path <- file.path(TH_OUTPUT, "temporal_holdout_manifest.csv")
readr::write_csv(manifest, manifest_path)

counts_path <- file.path(TH_OUTPUT, "temporal_holdout_slice_counts.csv")
readr::write_csv(manifest |> dplyr::select(slice, year_min, year_max, n_articles), counts_path)


# ---------------------------------------------------------
# Discovery target tables for downstream validation plotting
# ---------------------------------------------------------
disc_row <- manifest |> dplyr::filter(slice == "discovery")
if (nrow(disc_row) == 1) {
  disc_npmi <- readr::read_csv(disc_row$npmi_path[[1]], show_col_types = FALSE)
  disc_abc  <- readr::read_csv(disc_row$abc_path[[1]],  show_col_types = FALSE)
  disc_sign <- readr::read_csv(disc_row$sign_path[[1]], show_col_types = FALSE)
  disc_biom <- readr::read_csv(disc_row$biom_path[[1]], show_col_types = FALSE)

  disc_hotspot_targets <- disc_npmi |>
    dplyr::mutate(
      log2lift = log2(pmax(lift, 1e-12)),
      rank_npmi = dplyr::dense_rank(dplyr::desc(npmi)),
      rank_lift = dplyr::dense_rank(dplyr::desc(log2lift)),
      hotspot_score = -(rank_npmi + rank_lift)
    ) |>
    dplyr::group_by(C) |>
    dplyr::arrange(dplyr::desc(hotspot_score), dplyr::desc(npmi), .by_group = TRUE) |>
    dplyr::mutate(discovery_rank = dplyr::row_number()) |>
    dplyr::slice_head(n = HOTSPOT_TOP_K) |>
    dplyr::ungroup() |>
    dplyr::select(C, A, npmi, lift, hotspot_score, discovery_rank)

  disc_bridge_targets <- disc_abc |>
    dplyr::filter(C == "AE-ILD") |>
    dplyr::arrange(dplyr::desc(score_q)) |>
    dplyr::slice_head(n = 100) |>
    dplyr::group_by(B) |>
    dplyr::summarise(
      B_support = sum(score_q, na.rm = TRUE),
      triads = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(B_support), dplyr::desc(triads)) |>
    dplyr::mutate(discovery_rank = dplyr::row_number()) |>
    dplyr::slice_head(n = BRIDGE_TOP_K)

  disc_signed_primary <- disc_sign |>
    dplyr::filter(articles >= SIGNED_MIN_ARTICLES, balance != 0) |>
    dplyr::arrange(C, dplyr::desc(articles), dplyr::desc(abs(balance)), A) |>
    dplyr::mutate(discovery_rank = dplyr::row_number())

  disc_biomarker_core <- tidyr::expand_grid(
    A = BIOMARKER_CORE,
    C = OUTCOME_GROUPS
  ) |>
    dplyr::left_join(
      disc_biom |>
        dplyr::group_by(C) |>
        dplyr::arrange(dplyr::desc(score), dplyr::desc(evidence), .by_group = TRUE) |>
        dplyr::mutate(discovery_rank = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::select(A, C, discovery_rank, score, evidence),
      by = c("A","C")
    )

  readr::write_csv(disc_hotspot_targets, file.path(TH_OUTPUT, "discovery_targets_hotspots.csv"))
  readr::write_csv(disc_bridge_targets,  file.path(TH_OUTPUT, "discovery_targets_bridge_B.csv"))
  readr::write_csv(disc_signed_primary,  file.path(TH_OUTPUT, "discovery_targets_signed_primary.csv"))
  readr::write_csv(disc_biomarker_core,  file.path(TH_OUTPUT, "discovery_targets_biomarker_core.csv"))
}

readme_path <- file.path(TH_OUTPUT, "README_build.txt")
writeLines(c(
  "Temporal holdout build completed.",
  paste0("ARTICLES_FILE: ", ARTICLES_FILE),
  paste0("DIC_FILE     : ", DIC_FILE),
  paste0("Manifest     : ", manifest_path),
  paste0("Counts       : ", counts_path),
  "Discovery target tables written to output/ (hotspots, bridge_B, signed_primary, biomarker_core)"
), con = readme_path)

log_msg("WROTE manifest:", manifest_path)
log_msg("WROTE counts  :", counts_path)
log_msg("=== DONE 17_TimeSlice_BuildAndAnalyze ===")
