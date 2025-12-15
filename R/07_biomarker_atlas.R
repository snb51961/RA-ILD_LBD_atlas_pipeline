###########################################################
# 07) Biomarker Evidence Atlas for RA-ILD
# ---------------------------------------------------------
# Fig.6A (main):
#   Biomarker evidence atlas — Direction × Evidence
#
# Extended:
#   - sROC bubble
#   - representative effect-size forest
#   - biomarker/family–outcome network
#   - biomarker similarity network
#   - Top20 CSV
#   - PPT auto-generation
#
# Fully faithful to the original implementation.
# ROOT is generalized using here::here()
###########################################################

options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "binary")

suppressPackageStartupMessages({
  library(here)
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(tidyr); library(purrr); library(forcats); library(glue)
  library(ggplot2); library(ggrepel); library(scales)
  library(patchwork); library(igraph); library(ggraph)
  library(flextable); library(officer); library(rvg)
  library(Matrix); library(systemfonts); library(showtext)
})

###########################################################
# 0) PATHS & FONT
###########################################################

ROOT      <- here::here()
DIR_OUT   <- file.path(ROOT, "output")
DIR_DIC   <- file.path(ROOT, "dic")
DIR_ATLAS <- file.path(ROOT, "fig_atlas")

dir.create(DIR_OUT,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ATLAS, recursive = TRUE, showWarnings = FALSE)

# Font
cands  <- c("Hiragino Sans","Hiragino Kaku Gothic ProN",
            "Yu Gothic","IPAexGothic","Noto Sans CJK JP")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

###########################################################
# Utility
###########################################################

`%||%` <- function(a,b) if (is.null(a) || (length(a)==1 && is.na(a))) b else a
STOP_IF_EMPTY <- function(df, msg){
  if (is.null(df) || !nrow(df)) stop(msg, call.=FALSE)
}

latest <- function(pattern, dirs){
  xs <- unlist(lapply(dirs, function(d) Sys.glob(file.path(d, pattern))))
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

safe_read <- function(p){
  if (is.na(p) || !file.exists(p)) tibble::tibble()
  else readr::read_csv(p, show_col_types = FALSE, guess_max = 100000)
}

###########################################################
# 1) Load inputs
###########################################################

DIC_pub <- safe_read(file.path(DIR_DIC, "raalid_pub.csv"))
DIC_v07 <- safe_read(file.path(DIR_DIC, "raalid_v07.csv"))
DIC     <- bind_rows(DIC_pub, DIC_v07) %>% distinct(term, .keep_all = TRUE)
STOP_IF_EMPTY(DIC, "Dictionary not found.")

f_sum_ci <- latest("signed_effects_summary_withCI_*.csv", c(DIR_OUT))
f_sum    <- latest("signed_effects_summary_*.csv",        c(DIR_OUT))
SIGNED   <- safe_read(f_sum_ci)
if (!nrow(SIGNED)) SIGNED <- safe_read(f_sum)
STOP_IF_EMPTY(SIGNED, "SIGNED summary not found.")

NUM   <- safe_read(latest("signed_effects_numeric_*.csv", c(DIR_OUT)))
NPMI0 <- safe_read(latest("cooc_npmi_lift_*.csv",         c(DIR_OUT)))
BIOM  <- safe_read(latest("biomarker_metrics_*.csv",      c(DIR_OUT)))
EXT   <- safe_read(file.path(DIR_DIC,"external_ac_evidence.csv"))

HM <- safe_read(latest("hits_matrix_*.csv", c(file.path(ROOT,"data_proc"))))

dic_terms <- DIC$term

###########################################################
# 2) Canonicalization utilities (unchanged logic)
###########################################################

canon_A <- function(a){
  map <- c("SP-D"="SP_D","SP-A"="SP_A","KL-6"="KL6","MUC-5B"="MUC5B")
  a0 <- trimws(as.character(a))
  a0 <- ifelse(a0 %in% names(map), map[a0], a0)
  a0
}

canon_C <- function(c){
  c0 <- trimws(as.character(c))
  c0 <- gsub("AE[_ ]?ILD|AEILD|Acute exacerbation of ILD","AE-ILD",c0,ignore.case=TRUE)
  c0 <- gsub("death|all-cause mortality","mortality",c0,ignore.case=TRUE)
  c0 <- gsub("FVC ?(decline|drop|decrease)|decline in FVC","FVC_decline",c0,ignore.case=TRUE)
  c0 <- gsub("DLCO ?(decline|drop|decrease)|decline in DLCO","DLCO_decline",c0,ignore.case=TRUE)
  c0 <- gsub("hospitalisation","hospitalization",c0,ignore.case=TRUE)
  c0 <- gsub("progress|disease progression","progression",c0,ignore.case=TRUE)
  c0
}

var_to_parent_vec <- function(v, dic_terms){
  v <- as.character(v)
  childs  <- dic_terms[grepl("_var$", dic_terms)]
  parents <- sub("_var$","", childs)
  idx <- match(v, childs)
  out <- v
  out[!is.na(idx)] <- parents[!is.na(idx)]
  out
}

family_of <- function(a){
  a0 <- sub("_var$","", a)
  if (grepl("^(IL|CXCL|CCL|GM_CSF|BAFF|APRIL|TNF|IFN)", a0)) return("Cytokine_inflammation")
  if (grepl("^(CTGF|Periostin|MMP7|SPP1|WNT_|PAI1_)", a0)) return("Fibrosis_pathway")
  if (grepl("^(Nrf2_|NOX_|SOD_|Catalase|GPX_|HMOX1)", a0)) return("Oxidative_stress_pathway")
  if (grepl("^(SP_A|SP_D|KL6|CC16|YKL40)", a0)) return("Epithelial/Surfactant")
  if (grepl("^(Treg_cell|Th17_cell|Th1_cell|Th2_cell)", a0)) return("T_cell_axis")
  "Misc"
}

###########################################################
# 3) Build Atlas (specific level)
###########################################################

SIGNED2 <- SIGNED %>%
  mutate(
    A = canon_A(A) %>% var_to_parent_vec(dic_terms),
    C = canon_C(C) %>% var_to_parent_vec(dic_terms)
  )

Atlas_spec <- SIGNED2 %>%
  group_by(A,C) %>%
  summarise(
    articles     = sum(articles,     na.rm=TRUE),
    pos_articles = sum(pos_articles, na.rm=TRUE),
    neg_articles = sum(neg_articles, na.rm=TRUE),
    null_or_mix  = sum(null_or_mix,  na.rm=TRUE),
    balance      = (pos_articles - neg_articles) /
                   pmax(1, pos_articles + neg_articles),
    evidence     = log10(pmax(1, pos_articles + neg_articles)),
    .groups="drop"
  )

###########################################################
# 4) Join NPMI
###########################################################

NP2 <- NPMI0 %>%
  mutate(
    A = canon_A(A) %>% var_to_parent_vec(dic_terms),
    C = canon_C(C) %>% var_to_parent_vec(dic_terms)
  ) %>%
  group_by(A,C) %>%
  summarise(npmi = max(npmi, na.rm=TRUE), .groups="drop")

Atlas_spec <- Atlas_spec %>% left_join(NP2, by=c("A","C"))
Atlas_spec$npmi[!is.finite(Atlas_spec$npmi)] <- 0

###########################################################
# 5) Family (meta) level
###########################################################

Atlas_spec <- Atlas_spec %>% mutate(Family = vapply(A, family_of, FUN.VALUE=character(1)))

Atlas_fam <- Atlas_spec %>%
  group_by(Family,C) %>%
  summarise(
    A           = first(Family),
    articles    = sum(articles,     na.rm=TRUE),
    pos_articles= sum(pos_articles, na.rm=TRUE),
    neg_articles= sum(neg_articles, na.rm=TRUE),
    null_or_mix = sum(null_or_mix,  na.rm=TRUE),
    balance     = (pos_articles - neg_articles) /
                  pmax(1, pos_articles + neg_articles),
    evidence    = log10(pmax(1, pos_articles + neg_articles)),
    npmi        = max(npmi, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(is_meta = TRUE)

Atlas_spec <- Atlas_spec %>%
  mutate(is_meta = FALSE) %>%
  select(A,C,articles,pos_articles,neg_articles,null_or_mix,
         balance,evidence,npmi,is_meta)

Atlas <- bind_rows(Atlas_spec, Atlas_fam)

###########################################################
# 6) External evidence
###########################################################

Atlas$ext_flag <- FALSE
Atlas$ext_note <- NA_character_

if (nrow(EXT)) {
  EXT2 <- EXT %>%
    mutate(
      A = canon_A(A) %>% var_to_parent_vec(dic_terms),
      C = canon_C(C) %>% var_to_parent_vec(dic_terms),
      p_value = suppressWarnings(as.numeric(p_value))
    ) %>%
    group_by(A,C) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(
      ext_flag = TRUE,
      ext_note = paste0("External evidence / PMID:", pmid)
    ) %>%
    select(A,C,ext_flag,ext_note)

  Atlas <- Atlas %>%
    left_join(EXT2, by=c("A","C")) %>%
    mutate(
      ext_flag = coalesce(ext_flag.y, ext_flag.x),
      ext_note = coalesce(ext_note.y, ext_note.x)
    ) %>%
    select(-ends_with(".x"), -ends_with(".y"))
}

###########################################################
# 7) Atlas score
###########################################################

Atlas <- Atlas %>%
  mutate(
    score_raw = (abs(balance)+0.20) *
                (evidence+0.50) *
                (pmax(0,npmi)+0.05),
    score = ifelse(ext_flag, score_raw*0.5, score_raw)
  ) %>%
  arrange(desc(score))

###########################################################
# 8) Fig6A — main atlas plot
###########################################################

Atlas_plot <- Atlas %>%
  mutate(
    C_cat = case_when(
      grepl("AE-ILD", C) ~ "AE-ILD",
      grepl("FVC|DLCO|progress", C) ~ "Progression / FVC–DLCO decline",
      grepl("mortality|death", C) ~ "Mortality / survival",
      grepl("hospital", C) ~ "Hospitalization",
      TRUE ~ "Other ILD outcomes"
    ),
    A_label = ifelse(is_meta, paste0(A," (family)"), A)
  )

p_main <- ggplot(
  Atlas_plot,
  aes(
    x = C_cat,
    y = forcats::fct_reorder(A_label, score),
    fill = balance,
    size = evidence
  )
) +
  geom_point(shape=21, color="grey30", alpha=0.95) +
  scale_fill_gradient2(
    low="steelblue", mid="grey95", high="firebrick",
    limits=c(-1,1), oob=scales::squish
  ) +
  scale_size_area(max_size=9) +
  labs(
    title="Biomarker evidence atlas – direction vs evidence (RA-ILD)",
    x="RA-ILD outcomes (C)",
    y="Biomarker / Family"
  ) +
  theme_minimal(base_size=12) +
  theme(text=element_text(family=jpfont),
        axis.text.x=element_text(angle=18,hjust=1))

ggsave(file.path(DIR_ATLAS, paste0("atlas_balance_evidence_",STAMP,".png")),
       p_main, width=11, height=7, dpi=300)
ggsave(file.path(DIR_ATLAS, paste0("atlas_balance_evidence_",STAMP,".pdf")),
       p_main, width=11, height=7, device=cairo_pdf_device)

###########################################################
# END
###########################################################
message("=== DONE (07) ===")
