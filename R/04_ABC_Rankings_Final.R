# =========================================================
# 04_ABC_Rankings_Final.R  (public-ready, no-mixing design)
# - ABC rankings (A->B->C) with lift/NPMI thresholds
# - FDR(BH), IDF(B), TERM_PRIOR
# - Evidence PMIDs for top pairs (AB / BC)
# - Network edges for downstream plotting (no figure here)
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("dplyr","readr","stringr","tibble","purrr"))

# -------------------------
# 0) Parameters (keep consistent with your hardened defaults)
# -------------------------
LIFT_MIN_AB <- 1.20     # A→B
NPMI_MIN_AB <- 0.05
LIFT_MIN_BC <- 1.20     # B→C
NPMI_MIN_BC <- 0.05
NOVELTY_LAMBDA <- 1.0   # AC penalty strength
MAX_PMIDS_PER_PAIR <- 15
TOP_EVIDENCE_PAIRS <- 30

# Classes (same as your code)
A_CLASSES <- c("drug","bio","gene","exposure","host","event","nonpharm","vaccine","population","system","strategy")
B_CLASSES <- c("pattern","radiology","airway","qct","pft","activity","biomarker","molecular","microbiome","cell",
               "complication","phenotype","trend")

# TERM_PRIOR (same spirit; can be extended later)
TERM_PRIOR <- tibble::tibble(
  term  = c("Oxidative_stress","oxidative_stress","Oxidative_stress_var"),
  prior = c(0.2, 0.2, 0.2)
)

# -------------------------
# 1) Load inputs (STRICT, tagged)
# -------------------------
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_hits_tag)) {
  stop("Tagged hits_matrix not found: ", f_hits_tag,
       "\nRun 02_build_hits_matrix.R with tagged saving enabled.")
}
df <- readr::read_csv(f_hits_tag, show_col_types = FALSE)
log_msg("04 loaded hits_matrix:", f_hits_tag, " n=", nrow(df))

dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
dic <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","regex","class") %in% names(dic))) stop("Dictionary must have term/regex/class columns.")
log_msg("04 using dictionary:", basename(dic_path), " n=", nrow(dic))
log_msg("DIC head:", paste(head(dic$term, 6), collapse=" / "))

# _var mapping (exclude child terms from ABC sets, as you do)
var_childs <- unique(dic$term[grepl("_var$", dic$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) %>%
  dplyr::distinct(child, .keep_all = TRUE)

# -------------------------
# 2) Term sets A/B/C
# -------------------------
hit_cols <- names(df)[startsWith(names(df), "hit__")]
term_of  <- function(x) sub("^hit__", "", x)

A_set_full <- intersect(term_of(hit_cols), dic$term[dic$class %in% A_CLASSES])
B_set_full <- intersect(term_of(hit_cols), dic$term[dic$class %in% B_CLASSES])
A_set <- setdiff(A_set_full, term_groups$child)
B_set <- setdiff(B_set_full, term_groups$child)

C_set <- intersect(term_of(hit_cols), dic$term[dic$class %in% c("outcome")])
stopifnot(length(C_set) > 0)
C_main <- if ("AE-ILD" %in% C_set) "AE-ILD" else C_set[1]
log_msg(sprintf("ABC targets |A|=%d |B|=%d C_main=%s", length(A_set), length(B_set), C_main))

# Safe hit getter (length always matches df rows)
get_hit <- function(t) {
  nm <- paste0("hit__", t)
  if (nm %in% names(df)) as.integer(ifelse(is.na(df[[nm]]), 0L, df[[nm]])) else integer(nrow(df))
}
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

# Prior and IDF
prior_of <- function(t){
  v <- TERM_PRIOR$prior[match(t, TERM_PRIOR$term)]
  ifelse(is.na(v), 1.0, v)
}

N_docs <- nrow(df)
df_B <- sapply(B_set, function(b) sum(binv(get_hit(b))))
idf_B <- log((N_docs + 1) / (df_B + 1)) + 1
idfB_of <- function(b) idf_B[[b]]

# Pair statistics
pair_stats <- function(x, y){
  x <- binv(x); y <- binv(y)
  n11 <- sum(x==1 & y==1); n10 <- sum(x==1 & y==0); n01 <- sum(x==0 & y==1); n00 <- sum(x==0 & y==0)
  N   <- n11+n10+n01+n00
  OR  <- ((n11+0.5)*(n00+0.5))/((n10+0.5)*(n01+0.5))
  p   <- tryCatch(suppressWarnings(fisher.test(matrix(c(n11,n10,n01,n00),2,byrow=TRUE))$p.value),
                  error=function(e) 1)
  pA  <- (n11+n10)/N; pY <- (n11+n01)/N; pAc <- n11/N
  pmi <- if (pAc>0 && pA>0 && pY>0) log((pAc)/(pA*pY)) else 0
  npmi<- if (pAc>0) pmi/(-log(pAc)) else 0
  lift<- if (pA>0 && pY>0) (pAc)/(pA*pY) else 0
  list(or=as.numeric(OR), p=p, n11=n11, lift=lift, npmi=npmi)
}
safelog <- function(z){ x <- suppressWarnings(log(z)); x[!is.finite(x)] <- 0; x }

# -------------------------
# 3) ABC ranking core
# -------------------------
res <- list()
c_vec <- get_hit(C_main)

for (a in A_set){
  a_vec <- get_hit(a)
  for (b in B_set){
    b_vec <- get_hit(b)
    AB <- pair_stats(a_vec, b_vec)
    BC <- pair_stats(b_vec, c_vec)
    AC <- pair_stats(a_vec, c_vec)

    passAB <- (AB$lift >= LIFT_MIN_AB) && (AB$npmi >= NPMI_MIN_AB)
    passBC <- (BC$lift >= LIFT_MIN_BC) && (BC$npmi >= NPMI_MIN_BC)
    if (!(passAB && passBC)) next

    res[[length(res)+1]] <- data.frame(
      A=a, B=b, C=C_main,
      AB_or=AB$or, AB_p=AB$p, AB_n11=AB$n11, AB_lift=AB$lift, AB_npmi=AB$npmi,
      BC_or=BC$or, BC_p=BC$p, BC_n11=BC$n11, BC_lift=BC$lift, BC_npmi=BC$npmi,
      AC_or=AC$or, AC_p=AC$p, AC_n11=AC$n11,
      stringsAsFactors = FALSE
    )
  }
}

tab <- if (!length(res)) {
  tibble::tibble()
} else {
  dplyr::bind_rows(res) |>
    dplyr::mutate(
      AB_log  = safelog(AB_or),
      BC_log  = safelog(BC_or),
      AC_log  = safelog(AC_or),
      AB_sig  = -log10(pmax(AB_p, 1e-300)),
      BC_sig  = -log10(pmax(BC_p, 1e-300)),
      AC_pen  = pmax(AC_log, 0) * as.integer(-log10(pmax(AC_p, 1e-300)) >= 3),
      AB_q    = p.adjust(AB_p, method="BH"),
      BC_q    = p.adjust(BC_p, method="BH"),
      AC_q    = p.adjust(AC_p, method="BH"),
      AB_qsig = -log10(pmax(AB_q, 1e-300)),
      BC_qsig = -log10(pmax(BC_q, 1e-300)),
      priorA  = vapply(A, prior_of, numeric(1)),
      priorB  = vapply(B, prior_of, numeric(1)),
      idfB    = vapply(B, idfB_of, numeric(1)),
      score_q = (AB_log*AB_qsig*idfB*priorA*priorB + BC_log*BC_qsig*idfB*priorB) - (AC_pen*NOVELTY_LAMBDA)
    ) |>
    dplyr::arrange(dplyr::desc(score_q)) |>
    dplyr::mutate(dplyr::across(c(AB_or,BC_or,AC_or,AB_lift,BC_lift,AB_npmi,BC_npmi,score_q), ~ round(., 4)))
}

# Output (tagged, fixed name)
out_rank <- file.path(DIR_TABLE, sprintf("abc_rankings_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(tab, out_rank)
log_msg("WROTE:", out_rank, " n=", nrow(tab))

# -------------------------
# 4) Evidence PMIDs (AB / BC) for top pairs
# -------------------------
pmids_for <- function(t) df$pmid[get_hit(t) == 1L]

topN <- min(TOP_EVIDENCE_PAIRS, nrow(tab))
top_pairs <- if (topN > 0) dplyr::slice_head(tab, n=topN) else tab

ev_AB <- if (nrow(top_pairs) > 0) {
  purrr::map_dfr(seq_len(nrow(top_pairs)), function(i){
    a <- top_pairs$A[i]; b <- top_pairs$B[i]
    tibble::tibble(A=a, B=b, pmid=head(intersect(pmids_for(a), pmids_for(b)), MAX_PMIDS_PER_PAIR))
  })
} else tibble::tibble(A=character(), B=character(), pmid=character())

ev_BC <- if (nrow(top_pairs) > 0) {
  purrr::map_dfr(seq_len(nrow(top_pairs)), function(i){
    b <- top_pairs$B[i]; c <- C_main
    tibble::tibble(B=b, C=c, pmid=head(intersect(pmids_for(b), pmids_for(c)), MAX_PMIDS_PER_PAIR))
  })
} else tibble::tibble(B=character(), C=character(), pmid=character())

join_titles <- function(ev) {
  if (!nrow(ev)) return(ev)
  have <- intersect(c("pmid","year","journal","title"), names(df))
  if (!all(c("pmid","year","journal") %in% have)) return(ev)  # safe fallback
  dplyr::left_join(ev, df |> dplyr::select(dplyr::any_of(c("pmid","year","journal","title"))), by="pmid") |>
    dplyr::arrange(dplyr::desc(year))
}

evAB_path <- file.path(DIR_TABLE, sprintf("abc_evidence_AB_%s__%s.csv", CORPUS_TAG, DIC_TAG))
evBC_path <- file.path(DIR_TABLE, sprintf("abc_evidence_BC_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(join_titles(ev_AB), evAB_path)
readr::write_csv(join_titles(ev_BC), evBC_path)
log_msg("WROTE:", evAB_path, " n=", nrow(ev_AB))
log_msg("WROTE:", evBC_path, " n=", nrow(ev_BC))

# -------------------------
# 5) Network edges (safe even if empty)
# -------------------------
empty_edges <- tibble::tibble(from=character(), to=character(), w=integer(), kind=character())

edge_AB <- if (nrow(ev_AB) > 0) {
  ev_AB |> dplyr::count(A, B, name="w") |> dplyr::transmute(from=A, to=B, w=w, kind="AB")
} else empty_edges

edge_BC <- if (nrow(ev_BC) > 0) {
  ev_BC |> dplyr::count(B, name="w") |> dplyr::transmute(from=B, to=C_main, w=w, kind="BC")
} else empty_edges

edges <- dplyr::bind_rows(edge_AB, edge_BC)

out_edges <- file.path(DIR_TABLE, sprintf("abc_edges_top%d_%s__%s.csv", topN, CORPUS_TAG, DIC_TAG))
readr::write_csv(edges, out_edges)
log_msg("WROTE:", out_edges, " n=", nrow(edges))

log_msg("=== DONE 04_abc_rankings ===")
