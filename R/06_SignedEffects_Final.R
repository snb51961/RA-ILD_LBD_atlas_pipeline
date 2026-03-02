# =========================================================
# 06_signed_effects.R  (public-ready, no-mixing design)
# - Signed effects (A -> C) from abstracts/sentences with simple cues
# - Outputs: sentence-level, article-level, summary(+Wilson CI), numeric, biomarker metrics
# - Strict: uses tagged hits_matrix + dictionary fixed
# - Safe: articles text is restricted to PMIDs present in hits_matrix (prevents mixing)
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","tibble","purrr"))

# -------------------------
# 0) Parameters (tunable)
# -------------------------
OUTCOMES_FOR_SIGNED <- c("AE-ILD","progression","mortality")
PSEUDO_FUNC_OUTCOME <- TRUE

CHAR_WINDOW  <- 180
MAX_CONTEXTS <- 3
SAME_SENTENCE_ONLY <- FALSE
WINDOW <- 1

# If computation is heavy, you can cap A terms here (keep NULL to disable)
MAX_A_TERMS <- NULL  # e.g., 120

# -------------------------
# 1) Load inputs (STRICT, tagged hits)
# -------------------------
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_hits_tag)) {
  stop("Tagged hits_matrix not found: ", f_hits_tag,
       "\nRun 02_build_hits_matrix.R with tagged saving enabled.")
}
HM <- readr::read_csv(f_hits_tag, show_col_types = FALSE)
log_msg("06 loaded hits_matrix:", f_hits_tag, " n=", nrow(HM))

dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
dic <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","regex","class") %in% names(dic))) stop("Dictionary must have term/regex/class columns.")

# ---- Normalize dictionary class labels (public robustness) ----
dic <- dic %>%
  mutate(
    class_raw = class,
    class = tolower(trimws(class)),
    class = gsub("\\s+", "_", class),
    class = gsub("-", "_", class)
  ) %>%
  mutate(
    class = dplyr::case_when(
      class %in% c("outcomes","endpoint","endpoints") ~ "outcome",
      stringr::str_detect(class, "outcome|endpoint|exacerb|ae|mort|death|surviv|progress|decline|fvc|dlco") ~ "outcome",
      class %in% c("drug","drugs","treatment","therapy","medication","medicine") ~ "drug",
      stringr::str_detect(class, "drug|treat|therap|dmard|biologic|jak|inhibitor") ~ "drug",
      class %in% c("biomarker","marker","protein","molecular") ~ "bio",
      stringr::str_detect(class, "bio|biomark|marker|protein|molecular") ~ "bio",
      class %in% c("gene","genetics","variant","snp") ~ "gene",
      stringr::str_detect(class, "gene|genetic|variant|snp|gwas") ~ "gene",
      TRUE ~ class
    )
  )
log_msg("06 using dictionary:", basename(dic_path), " n=", nrow(dic))

# -------------------------
# 2) Load articles text safely (prevent mixing)
# -------------------------
# We deliberately restrict to PMIDs that exist in hits_matrix.
# If a wrong articles file is loaded, match-rate will be low and we stop.

find_articles_csv <- function(){
  # Prefer a same-day legacy file if present
  f_today <- file.path(DIR_PROC, sprintf("articles_%s.csv", format(Sys.Date(), "%Y%m%d")))
  if (file.exists(f_today)) return(f_today)

  # Otherwise pick latest in data_proc
  cand <- Sys.glob(file.path(DIR_PROC, "articles_*.csv"))
  if (!length(cand)) return(NA_character_)
  cand[which.max(file.info(cand)$mtime)]
}

f_art <- find_articles_csv()
if (is.na(f_art) || !file.exists(f_art)) {
  stop("articles_*.csv not found in data_proc. Run 01_fetch_pubmed.R first.")
}
ART <- readr::read_csv(f_art, show_col_types = FALSE)

# minimal columns needed
needA <- c("pmid","title","abstract")
missA <- setdiff(needA, names(ART))
if (length(missA)) stop("articles file missing required columns: ", paste(missA, collapse=", "))

ART <- ART %>%
  mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    text = paste(title, abstract, sep=" ")
  )

HM <- HM %>% mutate(pmid = as.character(pmid))

# restrict to hits_matrix PMIDs (anti-mixing guard)
ART2 <- ART %>% semi_join(HM %>% select(pmid), by="pmid")
match_rate <- nrow(ART2) / max(1, nrow(HM))
log_msg(sprintf("06 articles match-rate vs hits_matrix: %.3f", match_rate))

if (match_rate < 0.90) {
  stop(
    "Low PMID match-rate between hits_matrix and articles (", round(match_rate,3), "). ",
    "This suggests you may be mixing runs. Re-run 01 then 02, and run 06 on the same run outputs."
  )
}

# merge pubtype/journal/year/title/text into one DF used for text scanning
DF <- HM %>%
  left_join(ART2 %>% select(pmid, title, abstract, text, year), by="pmid") %>%
  mutate(text = ifelse(is.na(text), "", text))

# Robust year handling (covers year.x/year.y after joins)
if ("year.x" %in% names(DF) || "year.y" %in% names(DF)) {
  DF <- DF %>%
    mutate(year = dplyr::coalesce(.data$year.x, .data$year.y)) %>%
    select(-dplyr::any_of(c("year.x","year.y")))
}
if (!("year" %in% names(DF))) DF$year <- NA_integer_
DF <- DF %>% mutate(year = suppressWarnings(as.integer(.data$year)))
DF <- DF %>% mutate(text = ifelse(is.na(.data$text), "", .data$text))

# Ensure year exists
if (!("year" %in% names(DF))) DF$year <- NA_integer_


# ---- HARD ASSERT: DF must have text (fail fast) ----
stopifnot("text" %in% names(DF))
log_msg("06 DF text nonempty=", sum(nzchar(DF$text)), " empty=", sum(!nzchar(DF$text)))

# ---- Export to GlobalEnv (RStudio Source uses local=TRUE by default) ----
assign("DF", DF, envir = .GlobalEnv)
assign("df", DF, envir = .GlobalEnv)  # backward compatibility / convenience
assign("dic", dic, envir = .GlobalEnv)

# -------------------------
# 3) Term sets
# -------------------------
hit_cols <- names(DF)[startsWith(names(DF), "hit__")]
terms_all <- sub("^hit__", "", hit_cols)

A_CLASSES <- c("drug","bio","gene","exposure","host","event","nonpharm","vaccine","population","system","strategy")
A_CLASSES_SIGNED <- unique(c(A_CLASSES, "biomarker","molecular","microbiome"))

A_TERMS <- intersect(terms_all, dic$term[dic$class %in% A_CLASSES_SIGNED])

if (!length(A_TERMS)) {
  log_msg("06: A_TERMS empty after class filter; falling back to all terms in hits_matrix.")
  A_TERMS <- terms_all
}

C_TERMS <- intersect(terms_all, dic$term[dic$class %in% c("outcome")])

if (!length(C_TERMS)) {
  log_msg("06: C_TERMS empty after class filter; trying outcomes by name match.")
  # fallback: use OUTCOMES_FOR_SIGNED if present in terms_all
  C_TERMS <- intersect(terms_all, OUTCOMES_FOR_SIGNED)
}
C_TERMS <- intersect(C_TERMS, OUTCOMES_FOR_SIGNED)

if (!length(C_TERMS)) {
  log_msg("06: No C_TERMS matched OUTCOMES_FOR_SIGNED; will fall back to any outcome.")
  C_TERMS <- intersect(terms_all, dic$term[dic$class %in% c("outcome")])

if (!length(C_TERMS)) {
  log_msg("06: C_TERMS empty after class filter; trying outcomes by name match.")
  # fallback: use OUTCOMES_FOR_SIGNED if present in terms_all
  C_TERMS <- intersect(terms_all, OUTCOMES_FOR_SIGNED)
}
}

if (isTRUE(PSEUDO_FUNC_OUTCOME)) {
  col01 <- function(nm) { if (nm %in% names(DF)) as.integer(DF[[nm]]==1L) else integer(nrow(DF)) }
  make_pseudo <- function(name, base1, base2, desc){
    b1 <- paste0("hit__", base1); b2 <- paste0("hit__", base2)
    v  <- as.integer((col01(b1)==1L) & (col01(b2)==1L))
    nm <- paste0("hit__", name)
    if (!nm %in% names(DF)) DF[[nm]] <<- v

    # Create a *real* text regex for the pseudo outcome so sentence extraction can find contexts.
    # Use dictionary regex for base terms if available; otherwise fall back to the raw term.
    r1 <- dic$regex[match(base1, dic$term)]
    r2 <- dic$regex[match(base2, dic$term)]
    if (is.na(r1) || !nzchar(r1)) r1 <- base1
    if (is.na(r2) || !nzchar(r2)) r2 <- base2
    pseudo_regex <- paste0("(", r1, ".{0,40}", r2, "|", r2, ".{0,40}", r1, ")")

    if (!name %in% dic$term) {
      dic <<- dplyr::bind_rows(dic, tibble::tibble(
        term=name,
        regex=pseudo_regex,
        class="outcome",
        notes=desc
      ))
    }
    name
  }
  oc1 <- make_pseudo("FVC_decline",  "FVC",  "decline", "functional outcome")
  oc2 <- make_pseudo("DLCO_decline", "DLCO", "decline", "functional outcome")
  C_TERMS <- unique(c(C_TERMS, oc1, oc2))
}

# ---- Diagnostics for empty outputs ----
log_msg("06 preview C_TERMS:", paste(head(C_TERMS, 10), collapse=", "))
log_msg("06 preview A_TERMS:", paste(head(A_TERMS, 10), collapse=", "))
if (all(C_TERMS %in% c("FVC_decline","DLCO_decline"))) {
  log_msg("06 NOTE: C_TERMS contains only pseudo functional outcomes; regex-based sentence matching may be sparse.")
}


if (!length(A_TERMS) || !length(C_TERMS)) stop("A_TERMS or C_TERMS is empty; cannot run signed effects.")

# Optional cap for speed
if (!is.null(MAX_A_TERMS) && length(A_TERMS) > MAX_A_TERMS) {
  # keep top by frequency in hits_matrix
  freqA <- sapply(A_TERMS, function(a){
    nm <- paste0("hit__", a); if (nm %in% names(DF)) sum(DF[[nm]], na.rm=TRUE) else 0L
  })
  A_TERMS <- names(sort(freqA, decreasing=TRUE))[seq_len(MAX_A_TERMS)]
  log_msg("06 capped A_TERMS to top by frequency:", length(A_TERMS))
}

log_msg(sprintf("06 sets |A|=%d |C|=%d", length(A_TERMS), length(C_TERMS)))
assign("A_TERMS", A_TERMS, envir = .GlobalEnv)
assign("C_TERMS", C_TERMS, envir = .GlobalEnv)


# -------------------------
# 4) Core helpers (from your hardened block)
# -------------------------
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v>0L))

get_hit <- function(t) {
  nm <- paste0("hit__", t)
  if (nm %in% names(DF)) as.integer(ifelse(is.na(DF[[nm]]), 0L, DF[[nm]])) else integer(nrow(DF))
}

# AE expansion regex
AE_EXPAND_RE <- "(?i)((acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation)|acute[-\\s]*on[-\\s]*chronic).{0,80}(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b)|(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b).{0,80}acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation))"
NEG_RE <- "(?i)\\b(no|not|without|neither|never|absence\\s+of|lack\\s+of)\\b"
HEDGE_RE <- "(?i)\\b(may|might|tend\\s+to|trend(s)?\\s+toward|borderline|possibly|suggests?)\\b"

split_sentences <- function(txt){
  if (is.na(txt) || !nzchar(txt)) return(character(0))
  t <- gsub("[\\r\\n]+"," ", txt)
  t <- gsub("([A-Za-z])\\.(\\s*[A-Za-z])", "\\1. \\2", t)
  xs <- unlist(strsplit(t, "(?<=[\\.!?;:])\\s+|\\s+(?=Results?:|Conclusions?:|Background:)", perl=TRUE))
  xs <- trimws(xs)
  xs[nzchar(xs)]
}

all_spans <- function(pattern, text){
  if (is.na(text) || !nzchar(text)) return(list(starts=integer(0), ends=integer(0)))
  m <- gregexpr(pattern, text, perl=TRUE, ignore.case=TRUE)
  if (length(m) == 0 || m[[1]][1] == -1) return(list(starts=integer(0), ends=integer(0)))
  starts <- as.integer(m[[1]]); lens <- attr(m[[1]], "match.length")
  list(starts=starts, ends=starts + lens - 1L)
}

span_contexts <- function(text, reA, reC, char_window=CHAR_WINDOW, max_ctx=MAX_CONTEXTS){
  if (is.na(text) || !nzchar(text)) return(character(0))
  sa <- all_spans(reA, text); sc <- all_spans(reC, text)
  if (!length(sa$starts) || !length(sc$starts)) return(character(0))
  comb <- expand.grid(i=seq_along(sa$starts), j=seq_along(sc$starts))
  comb$dist <- mapply(function(i,j){
    a_mid <- (sa$starts[i] + sa$ends[i]) / 2
    c_mid <- (sc$starts[j] + sc$ends[j]) / 2
    abs(a_mid - c_mid)
  }, comb$i, comb$j)
  comb <- comb[order(comb$dist), , drop=FALSE]

  ctxs <- character(0)
  used <- matrix(FALSE, nrow=length(sa$starts), ncol=length(sc$starts))
  for (k in seq_len(nrow(comb))){
    i <- comb$i[k]; j <- comb$j[k]
    if (used[i,j]) next
    used[i,] <- TRUE; used[,j] <- TRUE
    mid1 <- min(sa$starts[i], sc$starts[j])
    # FIXED BUG: use sc$ends[j] (not sc$ends[i])
    mid2 <- max(sa$ends[i],   sc$ends[j])
    bgn  <- max(1L, mid1 - char_window)
    end  <- min(nchar(text), mid2 + char_window)
    seg  <- substr(text, bgn, end)
    ctxs <- c(ctxs, seg)
    if (length(ctxs) >= max_ctx) break
  }
  unique(ctxs)
}

rx <- list(
  inc = "(?i)(associated\\s+with\\s+(an?\\s+)?increase(d)?|increase(d|s)?\\s+(in|of)|higher\\s+risk|elevat(ed|es)|worsen(ed|ing)|exacerbat(ed|es)|risk\\s+factor|predict(or|ive)\\s+of|trigger(s|ed)?)",
  dec = "(?i)(associated\\s+with\\s+(a\\s+)?decrease(d)?|decrease(d|s)?\\s+(in|of)|lower\\s+risk|reduc(ed|es)|protective|improv(ed|es|ement)|ameliorat(ed|es)|attenuat(ed|es)|inhibit(ed|s))",
  null= "(?i)(not\\s+associated|no\\s+(significant\\s+)?association|no\\s+significant\\s+(difference|effect)|nonsignificant|non[-\\s]?significant|did\\s+not\\s+(find|show)|failed\\s+to\\s+(show|demonstrate)|similar\\s+(rates?|risk)|comparable)",
  meas= "(?i)\\b((?:a\\s*OR|aOR|OR|HR|aHR|RR|IRR|SHR|adj(?:usted)?\\s*(?:odds|hazard|risk)\\s*ratio))\\b\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)",
  ci  = "(?i)95%\\s*CI\\s*[:=]?\\s*\\(?\\s*([0-9]+\\.?[0-9]*)\\s*[-–,]\\s*([0-9]+\\.?[0-9]*)\\s*\\)?",
  p   = "(?i)\\bp\\s*[<=>]\\s*0\\.?0*([0-9]\\d?)"
)
rx$auc   <- "(?i)\\bAUC\\b\\s*[:=]?\\s*([0-9]\\.?[0-9]*)"
rx$cut   <- "(?i)(cut[-\\s]*off|threshold)\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)\\s*(U\\s*/\\s*mL|U/mL|mg/L|ng/mL|pg/mL)?"
rx$sens  <- "(?i)(sensitivity)\\s*[:=]?\\s*([0-9]{1,3}\\.?[0-9]*)\\s*%?"
rx$spec  <- "(?i)(specificity)\\s*[:=]?\\s*([0-9]{1,3}\\.?[0-9]*)\\s*%?"
SERO_NEG_RE <- "(?i)\\b(seronegative|ACPA\\s*negative|anti-?CCP\\s*negative|RF\\s*negative)\\b"

des_w <- function(pubtype, year){
  # Defensive: accept length-0 / NULL / non-numeric
  if (length(year) == 0 || is.null(year)) year <- NA_integer_
  year <- suppressWarnings(as.integer(year))

  w <- 1
  if (!is.na(pubtype) && grepl("Randomized|Controlled Clinical Trial", pubtype, ignore.case=TRUE)) w <- w + 0.8
  if (!is.na(pubtype) && grepl("Cohort|Case-Control", pubtype, ignore.case=TRUE)) w <- w + 0.3
  if (!is.na(pubtype) && grepl("Case Reports?", pubtype, ignore.case=TRUE)) w <- w - 0.3
  if (!is.na(year)) w <- w + 0.1 * pmax(0, year - 2015) / 10
  pmax(w, 0.2)
}


infer_one <- function(s){
  null_flag <- grepl(rx$null, s, perl=TRUE)
  inc_flag  <- grepl(rx$inc,  s, perl=TRUE)
  dec_flag  <- grepl(rx$dec,  s, perl=TRUE)
  neg_flag  <- grepl(NEG_RE, s, perl=TRUE)
  if (neg_flag && (inc_flag || dec_flag)) { inc_flag <- FALSE; dec_flag <- FALSE; null_flag <- TRUE }
  hedge_flag <- grepl(HEDGE_RE, s, perl=TRUE)

  m  <- stringr::str_match(s, rx$meas)
  val<- suppressWarnings(as.numeric(if (is.null(m)) NA else m[,3]))
  meas_name <- if (!is.null(m) && !is.na(m[,2])) m[,2] else NA

  ci <- stringr::str_match(s, rx$ci)
  lo <- suppressWarnings(as.numeric(ci[,2])); hi <- suppressWarnings(as.numeric(ci[,3]))

  p_sig <- grepl(rx$p, s, perl=TRUE)
  m_auc <- stringr::str_match(s, rx$auc); AUC  <- suppressWarnings(as.numeric(m_auc[,2]))
  m_cut <- stringr::str_match(s, rx$cut); CUT  <- suppressWarnings(as.numeric(m_cut[,2])); UNIT <- ifelse(is.na(m_cut[,3]), NA_character_, m_cut[,3])
  m_sen <- stringr::str_match(s, rx$sens); SENS <- suppressWarnings(as.numeric(m_sen[,3]))
  m_spe <- stringr::str_match(s, rx$spec); SPEC <- suppressWarnings(as.numeric(m_spe[,3]))

  sig_by_ci <- if (!is.na(lo) && !is.na(hi)) if (lo>1 || hi<1) TRUE else FALSE else NA
  sig <- if (!is.na(sig_by_ci)) sig_by_ci else if (p_sig) TRUE else NA

  sgn <- NA_real_; weight <- 1
  if (!is.na(val)) sgn <- if (val > 1) +1 else if (val < 1) -1 else 0
  if (is.na(sgn)) {
    sign_cue <- if (null_flag) 0 else if (inc_flag && !dec_flag) +1 else if (dec_flag && !inc_flag) -1 else NA
    sgn <- sign_cue
  }
  if (!is.na(sig) && sig) weight <- weight + 1
  if (grepl("(?i)(adjusted|multivaria(te|ted))", s, perl=TRUE)) weight <- weight + 0.5
  if (hedge_flag) weight <- weight * 0.6

  seroneg_flag <- grepl(SERO_NEG_RE, s, perl=TRUE)
  if (is.na(sgn) && seroneg_flag) { sgn <- -1; weight <- weight * 0.8 }

  list(
    sign=sgn, weight=weight,
    inc=as.integer(inc_flag), dec=as.integer(dec_flag), null=as.integer(null_flag),
    measure=if (!is.na(meas_name)) meas_name else NA_character_,
    value=val, ci_low=lo, ci_high=hi, p_sig=as.integer(p_sig),
    auc=AUC, cutoff=CUT, unit=UNIT, sens=SENS, spec=SPEC
  )
}

analyze_pair <- function(a_term, c_term, window=WINDOW, same_sentence_only=SAME_SENTENCE_ONLY){
  mask <- which(get_hit(a_term)==1 & get_hit(c_term)==1)

  # Outcome expand: AE-ILD
  if (identical(c_term, "AE-ILD")) {
    c_extra <- grepl(AE_EXPAND_RE, DF$text, ignore.case=TRUE, perl=TRUE)
    mask <- which(get_hit(a_term)==1 & (get_hit(c_term)==1 | c_extra))
  }

  if (!length(mask)) return(NULL)

  # regex for A and C
  a_re <- dic$regex[dic$term==a_term][1]
  c_re <- if (identical(c_term, "AE-ILD")) AE_EXPAND_RE else dic$regex[dic$term==c_term][1]
  if (is.na(a_re) || is.na(c_re)) return(NULL)

  out_sent <- list(); out_art <- list()

  for (idx in mask){
    pmid <- DF$pmid[idx]; yr <- DF$year[idx]
    abs_txt <- DF$text[idx]
    if (is.na(abs_txt) || !nzchar(abs_txt)) next

    ss <- split_sentences(abs_txt)
    contexts <- character(0); src_weight <- 1

    Sa <- which(grepl(a_re, ss, ignore.case=TRUE, perl=TRUE))
    Sc <- which(grepl(c_re, ss, ignore.case=TRUE, perl=TRUE))

    if (length(Sa) && length(Sc)) {
      if (same_sentence_only) {
        S_both <- intersect(Sa, Sc)
        if (length(S_both)) contexts <- c(contexts, ss[S_both])
      } else {
        idxs <- sort(unique(c(Sa, Sc)))
        for (k in idxs){
          i1 <- max(1, k - window); i2 <- min(length(ss), k + window)
          seg <- paste(ss[i1:i2], collapse=" ")
          if ( any(grepl(a_re, seg, ignore.case=TRUE, perl=TRUE)) &&
               any(grepl(c_re, seg, ignore.case=TRUE, perl=TRUE)) ) {
            contexts <- c(contexts, seg)
          }
        }
      }
    }

    if (!length(contexts)) {
      ctx_span <- span_contexts(abs_txt, a_re, c_re, CHAR_WINDOW, MAX_CONTEXTS)
      if (length(ctx_span)) { contexts <- c(contexts, ctx_span); src_weight <- 0.8 }
    }

    if (!length(contexts)) {
      if (grepl(a_re, abs_txt, ignore.case=TRUE, perl=TRUE) &&
          grepl(c_re, abs_txt, ignore.case=TRUE, perl=TRUE)) {
        contexts <- c(contexts, abs_txt); src_weight <- 0.4
      }
    }
    if (!length(contexts)) next

    art_w <- des_w(DF$pubtype[idx], DF$year[idx])

    S <- dplyr::bind_rows(lapply(unique(contexts), function(s){
      r <- infer_one(s)
      r$weight <- r$weight * src_weight * art_w
      data.frame(
        pmid=pmid, year=yr, A=a_term, C=c_term, sentence=s,
        sign=r$sign, weight=r$weight, inc=r$inc, dec=r$dec, null=r$null,
        measure=r$measure, value=r$value, ci_low=r$ci_low, ci_high=r$ci_high,
        p_sig=r$p_sig, auc=r$auc, cutoff=r$cutoff, unit=r$unit, sens=r$sens, spec=r$spec,
        stringsAsFactors=FALSE
      )
    }))
    if (!nrow(S)) next
    out_sent[[length(out_sent)+1]] <- S

    sc <- sum(ifelse(is.na(S$sign), 0, S$sign * S$weight))
    lbl <- if (sc>0) "risk_up" else if (sc<0) "risk_down" else "no_effect_or_mixed"
    out_art[[length(out_art)+1]] <- data.frame(
      pmid=pmid, year=yr, A=a_term, C=c_term, n_sent=nrow(S),
      pos=sum(S$sign==1, na.rm=TRUE),
      neg=sum(S$sign==-1, na.rm=TRUE),
      null=sum(S$sign==0, na.rm=TRUE),
      score=sc, label=lbl, stringsAsFactors=FALSE
    )
  }

  list(
    sent = if (length(out_sent)) dplyr::bind_rows(out_sent) else NULL,
    art  = if (length(out_art))  dplyr::bind_rows(out_art)  else NULL
  )
}

# -------------------------
# 5) Run scanning
# -------------------------
sent_res <- list(); art_res <- list()
for (a in unique(A_TERMS)){
  for (c in unique(C_TERMS)){
    z <- analyze_pair(a, c, window=WINDOW, same_sentence_only=SAME_SENTENCE_ONLY)
    if (!is.null(z$sent)) sent_res[[length(sent_res)+1]] <- z$sent
    if (!is.null(z$art))  art_res[[length(art_res)+1]]  <- z$art
  }
}

signed_sent <- if (length(sent_res)) dplyr::bind_rows(sent_res) else
  tibble::tibble(pmid=character(), year=integer(), A=character(), C=character(), sentence=character(),
                 sign=double(), weight=double(), inc=integer(), dec=integer(), null=integer(),
                 measure=character(), value=double(), ci_low=double(), ci_high=double(), p_sig=integer(),
                 auc=double(), cutoff=double(), unit=character(), sens=double(), spec=double())

signed_art  <- if (length(art_res))  dplyr::bind_rows(art_res) else
  tibble::tibble(pmid=character(), year=integer(), A=character(), C=character(),
                 n_sent=integer(), pos=integer(), neg=integer(), null=integer(),
                 score=double(), label=character())

# -------------------------
# 6) Save outputs (tagged, fixed names)
# -------------------------
f_signed_sent <- file.path(DIR_TABLE, sprintf("signed_effects_sentence_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_signed_art  <- file.path(DIR_TABLE, sprintf("signed_effects_article_%s__%s.csv",  CORPUS_TAG, DIC_TAG))
readr::write_csv(signed_sent, f_signed_sent)
readr::write_csv(signed_art,  f_signed_art)
log_msg("WROTE:", f_signed_sent, " n=", nrow(signed_sent))
log_msg("WROTE:", f_signed_art,  " n=", nrow(signed_art))


# -------------------------
# =========================
# 07 PATCH: summary + Wilson CI (schema-fixed, legacy-compatible)
# =========================

# 1) Fallback when label is missing (generate from score)
if (nrow(signed_art) > 0) {
  if (!("label" %in% names(signed_art))) {
    if (!("score" %in% names(signed_art))) {
      stop("signed_art lacks both 'label' and 'score'. Cannot summarize signed effects.")
    }
    signed_art <- signed_art |>
      dplyr::mutate(
        label = dplyr::case_when(
          score > 0 ~ "risk_up",
          score < 0 ~ "risk_down",
          TRUE      ~ "no_effect_or_mixed"
        )
      )
  }
}

# 2) Aggregate into legacy schema (do not include pos_neg in outputs)
if (nrow(signed_art) > 0) {
  sum_tab <- signed_art |>
    dplyr::group_by(A, C) |>
    dplyr::summarise(
      articles     = dplyr::n(),
      pos_articles = sum(label == "risk_up", na.rm = TRUE),
      neg_articles = sum(label == "risk_down", na.rm = TRUE),
      null_or_mix  = sum(label == "no_effect_or_mixed", na.rm = TRUE),
      net_score    = sum(score, na.rm = TRUE),
      .groups      = "drop"
    ) |>
    dplyr::mutate(
      pos_ratio = pos_articles / pmax(1, pos_articles + neg_articles),
      balance   = (pos_articles - neg_articles) / pmax(1, pos_articles + neg_articles)
    ) |>
    dplyr::arrange(dplyr::desc(articles))
} else {
  # Ensure columns exist even when there are 0 rows (backward compatible)
  sum_tab <- tibble::tibble(
    A = character(), C = character(),
    articles = integer(),
    pos_articles = integer(), neg_articles = integer(), null_or_mix = integer(),
    net_score = double(),
    pos_ratio = double(), balance = double()
  )
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

if (nrow(sum_tab) > 0) {
  ci_mat <- t(mapply(wilson_balance_ci, sum_tab$pos_articles, sum_tab$neg_articles))
  sum_tab$balance_low  <- ci_mat[, "low"]
  sum_tab$balance_high <- ci_mat[, "high"]
} else {
  sum_tab$balance_low  <- double()
  sum_tab$balance_high <- double()
}

# 3) Fix column types and order (most important for reproducibility)
sum_tab <- sum_tab |>
  dplyr::mutate(
    articles     = as.integer(articles),
    pos_articles = as.integer(pos_articles),
    neg_articles = as.integer(neg_articles),
    null_or_mix  = as.integer(null_or_mix),
    net_score    = as.numeric(net_score),
    pos_ratio    = as.numeric(pos_ratio),
    balance      = as.numeric(balance),
    balance_low  = as.numeric(balance_low),
    balance_high = as.numeric(balance_high)
  ) |>
  dplyr::select(
    A, C,
    articles, pos_articles, neg_articles, null_or_mix,
    net_score, pos_ratio, balance,
    balance_low, balance_high
  )

# 4) Spec checks (fail fast here if needed)
need_cols <- c("A","C","articles","pos_articles","neg_articles","null_or_mix",
               "net_score","pos_ratio","balance","balance_low","balance_high")
stopifnot(all(need_cols %in% names(sum_tab)))


# 5) Save summary_withCI (tagged, fixed name)
f_signed_sum <- file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(sum_tab, f_signed_sum)
log_msg("WROTE:", f_signed_sum, " n=", nrow(sum_tab))

# -------------------------
# 8) Numeric extraction (sentence-level)
# -------------------------
num_cols <- c("pmid","year","A","C","measure","value","ci_low","ci_high","p_sig","auc","cutoff","unit","sens","spec","sentence")
signed_num <- signed_sent %>% dplyr::select(dplyr::any_of(num_cols))
f_signed_num <- file.path(DIR_TABLE, sprintf("signed_effects_numeric_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(signed_num, f_signed_num)
log_msg("WROTE:", f_signed_num, " n=", nrow(signed_num))

# -------------------------
# 9) Biomarker metrics (same as your set)
# -------------------------
biomarkers <- c("KL6","KL6_var","SP_D","SPD_var","NLR","NLR_var","CAR","CAR_var")
BM <- signed_sent %>%
  dplyr::filter(A %in% biomarkers) %>%
  dplyr::mutate(marker = A) %>%
  dplyr::select(marker, pmid, year, C, auc, cutoff, unit, sens, spec, measure, value, ci_low, ci_high, sentence) %>%
  dplyr::arrange(dplyr::desc(auc), dplyr::desc(sens), dplyr::desc(spec))
f_biom <- file.path(DIR_TABLE, sprintf("biomarker_metrics_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(BM, f_biom)
log_msg("WROTE:", f_biom, " n=", nrow(BM))

log_msg("=== DONE 06_signed_effects ===")

