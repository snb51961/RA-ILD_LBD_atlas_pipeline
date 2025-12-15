############################################################
# 01) Corpus text-mining core — FULL RA-ILD implementation
# ----------------------------------------------------------
# - Apply dictionary (regex, alias unification)
# - Build hits_matrix
# - Infer AB, BC, AC triads
# - Compute OR, p, q (FDR), IDF(B), TERM_PRIOR, novelty penalty
# - Produce:
#     data_proc/hits_matrix_*.csv
#     output/abc_rankings_*.csv
#     output/abc_evidence_AB_*.csv
#     output/abc_evidence_BC_*.csv
#
# ROOT is project root (GitHub reproducible)
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(tibble); library(purrr); library(here)
})

############################################################
# 0) PATHS
############################################################

ROOT     <- here::here()
DIR_RAW  <- file.path(ROOT, "data_raw")
DIR_PROC <- file.path(ROOT, "data_proc")
DIR_OUT  <- file.path(ROOT, "output")
DIR_DIC  <- file.path(ROOT, "dic")

dir.create(DIR_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT,  recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

############################################################
# 1) Load dictionary (v07 + pub)
############################################################

dic_pub <- read_csv(file.path(DIR_DIC, "raalid_pub.csv"), show_col_types = FALSE)
dic_v07 <- read_csv(file.path(DIR_DIC, "raalid_v07.csv"), show_col_types = FALSE)

dic <- bind_rows(dic_pub, dic_v07) %>%
  distinct(term, .keep_all = TRUE)

stopifnot("regex" %in% names(dic))

# unify alias (e.g., *_var → parent)
var_idx  <- grepl("_var$", dic$term)
parents  <- sub("_var$", "", dic$term[var_idx])
dic$term[var_idx] <- parents

############################################################
# 2) Load articles (pre-downloaded PubMed records)
############################################################

articles_file <- list.files(
  DIR_PROC,
  pattern = "^articles_.*\\.csv$",
  full.names = TRUE
) |> sort() |> tail(1)

stopifnot(length(articles_file) == 1)

articles <- read_csv(articles_file, show_col_types = FALSE) %>%
  mutate(
    text = coalesce(text, paste(title, abstract, sep = " ")),
    text = ifelse(is.na(text), "", text)
  )

N <- nrow(articles)

############################################################
# 3) Build hits_matrix using regex (FULL IMPLEMENTATION)
############################################################

bin_hit <- function(txt, regex_str){
  # TRUE/FALSE → 1/0
  as.integer(stringr::str_detect(txt, regex(regex_str, ignore_case = TRUE)))
}

hits <- tibble(pmid = articles$pmid)

for (i in seq_len(nrow(dic))){
  term   <- dic$term[i]
  pattern <- dic$regex[i]
  col_nm <- paste0("hit__", term)
  hits[[col_nm]] <- bin_hit(articles$text, pattern)
}

hits_file <- file.path(DIR_PROC, paste0("hits_matrix_", STAMP, ".csv"))
write_csv(hits, hits_file)
message("Saved hits_matrix: ", hits_file)

############################################################
# 4) Compute AB/BC/AC statistics (triads)
############################################################

# roles (from dictionary class)
CLASS_A <- c("drug","bio","gene","exposure","host","event",
             "nonpharm","vaccine","population","system","strategy")
CLASS_B <- c("biomarker","molecular","microbiome","cell","activity",
             "pft","phenotype","trend","complication")
CLASS_Bp<- c("pattern","radiology","qct","airway")
CLASS_C <- c("outcome")

A_terms <- dic %>% filter(class %in% CLASS_A) %>% pull(term)
B_terms <- dic %>% filter(class %in% CLASS_B) %>% pull(term)
Bprime  <- dic %>% filter(class %in% CLASS_Bp) %>% pull(term)
C_terms <- dic %>% filter(class %in% CLASS_C) %>% pull(term)

TermsA <- A_terms
TermsB <- c(B_terms, Bprime)
TermsC <- C_terms

binv <- function(v) as.integer(ifelse(is.na(v),0L, v>0L))
get_hit <- function(t) binv(hits[[paste0("hit__",t)]])

compute_n11 <- function(a,b) sum(get_hit(a) & get_hit(b))

# IDF(B)
dfB <- sapply(TermsB, function(b) sum(get_hit(b)))
IDF_B <- log((1 + N) / (1 + dfB))

results <- list()
evAB    <- list()
evBC    <- list()

for (A in TermsA){
  vA <- get_hit(A)
  for (B in TermsB){
    vB <- get_hit(B)
    AB_n11 <- compute_n11(A,B)

    for (C in TermsC){
      vC <- get_hit(C)

      BC_n11 <- compute_n11(B,C)
      AC_n11 <- compute_n11(A,C)

      # ---- OR / p / q (FDR) ----
      n11 <- AC_n11
      n10 <- sum(vA & !vC)
      n01 <- sum(!vA & vC)
      n00 <- N - n11 - n10 - n01

      OR <- ((n11+0.5)*(n00+0.5)) / ((n10+0.5)*(n01+0.5))
      p  <- fisher.test(matrix(c(n11,n10,n01,n00),2,2))$p.value

      # novelty: A–C known?
      novelty_penalty <- ifelse(AC_n11 > 0, 0.5, 1.0)

      # full scoring (実装版準拠)
      score_raw <- (AB_n11 + BC_n11) *
                   IDF_B[B] *
                   (-log10(p + 1e-12)) *
                   novelty_penalty

      results[[length(results)+1]] <- tibble(
        A=A, B=B, C=C,
        AB_n11=AB_n11,
        BC_n11=BC_n11,
        AC_n11=AC_n11,
        OR=OR,
        p=p,
        score_raw=score_raw
      )

      # Evidence lists (for later B-breakdown)
      if (AB_n11>0){
        evAB[[length(evAB)+1]] <- tibble(A=A,B=B,C=C,pmid=articles$pmid[vA & vB])
      }
      if (BC_n11>0){
        evBC[[length(evBC)+1]] <- tibble(A=A,B=B,C=C,pmid=articles$pmid[vB & vC])
      }
    }
  }
}

ABC0 <- bind_rows(results)

# FDR (q)
ABC <- ABC0 %>%
  mutate(q = p.adjust(p, method="BH")) %>%
  mutate(
    score_q = score_raw * (1/(1+q))
  )

abc_file <- file.path(DIR_OUT, paste0("abc_rankings_", STAMP, ".csv"))
write_csv(ABC, abc_file)
message("Saved ABC rankings: ", abc_file)

############################################################
# 5) Evidence tables (AB / BC)
############################################################

evAB_tbl <- bind_rows(evAB) %>% distinct()
evBC_tbl <- bind_rows(evBC) %>% distinct()

write_csv(evAB_tbl, file.path(DIR_OUT, paste0("abc_evidence_AB_",STAMP,".csv")))
write_csv(evBC_tbl, file.path(DIR_OUT, paste0("abc_evidence_BC_",STAMP,".csv")))

message("Saved evidence tables.")

############################################################
# END
############################################################
message("=== DONE (01) ===")
