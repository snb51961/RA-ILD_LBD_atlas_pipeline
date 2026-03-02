# =========================================================
# 07_ae_ratio_trendbreak.R  (public-ready, no-mixing design)
# - Detect change-points in AE-ILD ratio trend for top A terms
# - ratio(year) = (A & AE-ILD) / A
# - A terms exclude outcomes (dictionary class-based) to avoid tautology
# - Filters small denominators per year to avoid 0/100% spikes
# - QC plot: Top3 with red dashed lines, point size = den
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","tibble","purrr","changepoint","ggplot2","scales","patchwork"))

# -------------------------
# Parameters
# -------------------------
TOPK_FOR_RATIO_CPT <- 15
MIN_YEARS_FOR_CPT  <- 6
MIN_DEN_PER_YEAR   <- 5    # <- IMPORTANT: removes unstable spike years

# -------------------------
# Load inputs (STRICT, tagged)
# -------------------------
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_hits_tag)) {
  stop("Tagged hits_matrix not found: ", f_hits_tag,
       "\nRun 02_build_hits_matrix.R with tagged saving enabled.")
}
df <- readr::read_csv(f_hits_tag, show_col_types = FALSE)
log_msg("07 loaded hits_matrix:", f_hits_tag, " n=", nrow(df))

dic_path <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")
if (!file.exists(dic_path)) stop("Dictionary not found: ", dic_path)
dic <- readr::read_csv(dic_path, show_col_types = FALSE)
if (!all(c("term","regex","class") %in% names(dic))) stop("Dictionary must have term/regex/class columns.")
log_msg("07 using dictionary:", basename(dic_path), " n=", nrow(dic))

# -------------------------
# Helpers
# -------------------------
hit_cols   <- names(df)[startsWith(names(df), "hit__")]
terms_all  <- sub("^hit__", "", hit_cols)

binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

get_hit <- function(t){
  nm <- paste0("hit__", t)
  if (nm %in% names(df)) binv(df[[nm]]) else integer(nrow(df))
}

make_yearly <- function(v) {
  df |>
    dplyr::transmute(year = as.integer(year), hit = as.integer(v)) |>
    dplyr::filter(!is.na(year)) |>
    dplyr::group_by(year) |>
    dplyr::summarise(n = sum(hit, na.rm=TRUE), .groups="drop") |>
    dplyr::arrange(year)
}

# -------------------------
# Preconditions
# -------------------------
if (!("AE-ILD" %in% terms_all)) {
  log_msg("07: AE-ILD not found in hits_matrix -> skipped.")
  quit(save="no", status=0)
}

# -------------------------
# Define A candidates (IMPORTANT FIX)
# - exclude outcomes (class == outcome)
# - exclude _var children
# -------------------------
C_terms_all <- dic$term[dic$class == "outcome"]

var_childs <- unique(dic$term[grepl("_var$", dic$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) %>%
  dplyr::distinct(child, .keep_all = TRUE)

A_all <- setdiff(terms_all, C_terms_all)
A_all <- setdiff(A_all, term_groups$child)

# Safety: AE-ILD itself should never be in A_all
A_all <- setdiff(A_all, "AE-ILD")

log_msg("07 A_all size (outcomes excluded):", length(A_all))

# -------------------------
# Pick top A by total co-occurrence with AE-ILD
# -------------------------
tot_AE <- sapply(A_all, function(a){
  sum(get_hit(a) & get_hit("AE-ILD"))
})

topA <- names(sort(tot_AE, decreasing=TRUE))[seq_len(min(TOPK_FOR_RATIO_CPT, length(tot_AE)))]
log_msg("07 topA:", paste(topA, collapse=", "))

# -------------------------
# Build TB with changepoints
# -------------------------
tb_list <- list()

for (a in topA){
  yr_co <- make_yearly(get_hit(a) & get_hit("AE-ILD"))  # numerator
  yr_a  <- make_yearly(get_hit(a))                      # denominator

  D <- yr_co |>
    dplyr::left_join(yr_a |> dplyr::rename(den = n), by="year") |>
    dplyr::mutate(ratio = dplyr::if_else(den > 0, n/den, NA_real_)) |>
    dplyr::filter(!is.na(ratio)) |>
    dplyr::filter(den >= MIN_DEN_PER_YEAR)               # IMPORTANT: stabilize

  cps <- tryCatch({
    if (nrow(D) >= MIN_YEARS_FOR_CPT) {
      cpt <- changepoint::cpt.mean(D$ratio, method="PELT", penalty="BIC")
      idx <- changepoint::cpts(cpt)
      if (length(idx)) D$year[idx] else integer(0)
    } else integer(0)
  }, error=function(e) integer(0))

  tb_list[[length(tb_list)+1]] <- tibble::tibble(
    term  = a,
    year  = D$year,
    num   = D$n,
    den   = D$den,
    ratio = D$ratio,
    cpts  = paste(cps, collapse=";")
  )
}

TB <- dplyr::bind_rows(tb_list)

# -------------------------
# Output TB (tagged)
# -------------------------
out_tb <- file.path(
  DIR_TABLE,
  sprintf("ae_ratio_trendbreaks_top%d_%s__%s.csv", TOPK_FOR_RATIO_CPT, CORPUS_TAG, DIC_TAG)
)
readr::write_csv(TB, out_tb)

# -------------------------
# QC plot: Top3 terms with changepoints (red dashed), point size = den
# -------------------------
if (nrow(TB) > 0) {

  # Choose top3 with a stable scoring rule
  score <- TB %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      max_ratio = max(ratio, na.rm = TRUE),
      n_years   = dplyr::n(),
      sc        = max_ratio * log1p(n_years),
      .groups="drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sc))

  sel3 <- head(score$term, 3)
  P <- list()

  for (t in sel3) {
    D <- TB %>% dplyr::filter(term == t) %>% dplyr::arrange(year)

    cps_chr <- D$cpts[!is.na(D$cpts)][1]
    cps_vec <- numeric(0)
    if (!is.na(cps_chr) && nzchar(cps_chr)) {
      cps_vec <- suppressWarnings(as.numeric(unlist(strsplit(cps_chr, ";"))))
      cps_vec <- cps_vec[is.finite(cps_vec)]
    }

    p <- ggplot2::ggplot(D, ggplot2::aes(x = year, y = ratio)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(size = den), alpha = 0.85) +
      { if (length(cps_vec)) ggplot2::geom_vline(xintercept = cps_vec, linetype="dashed", color="red") } +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::scale_size_area(max_size = 6) +
      ggplot2::labs(
        title = paste0(t, " × AE-ILD ratio (den-filtered, size=den)"),
        subtitle = paste0("Filter: den >= ", MIN_DEN_PER_YEAR, " | changepoints: red dashed"),
        x = "Year",
        y = "P(A & AE-ILD) / P(A)",
        size = "den"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    P[[length(P)+1]] <- p
  }

  if (length(P)) {
    p_all <- patchwork::wrap_plots(P, ncol = 1)

    f_png <- file.path(DIR_FIG2, sprintf("ae_ratio_trendbreak_top3_%s__%s.png", CORPUS_TAG, DIC_TAG))
    f_pdf <- file.path(DIR_FIG2, sprintf("ae_ratio_trendbreak_top3_%s__%s.pdf", CORPUS_TAG, DIC_TAG))

    ggplot2::ggsave(f_png, p_all, width = 9, height = 10, dpi = 300)
    ggplot2::ggsave(f_pdf, p_all, width = 9, height = 10, device = grDevices::cairo_pdf)

    log_msg("WROTE QC plot:", f_png)
    log_msg("WROTE QC plot:", f_pdf)
  }
} else {
  log_msg("QC plot skipped: TB is empty (nrow=0).")
}

log_msg("WROTE:", out_tb, " n=", nrow(TB))
log_msg("=== DONE 07_ae_ratio_trendbreak ===")
