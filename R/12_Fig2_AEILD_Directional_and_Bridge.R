# =========================================================
# 12_Fig2_SignedEffects_Summary_Final.R
# Fix: Fig2B now uses topA_support directly and fills missing balance/articles
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("dplyr","readr","stringr","forcats","ggplot2","scales","patchwork"))


# -------------------------
# Reproducibility log
# -------------------------
RUN_STAMP <- format(Sys.time(), "%Y%m%d-%H%M%S")
RUN_LOG <- file.path(DIR_LOG, paste0("run12_figure2_", RUN_STAMP, ".txt"))
writeLines(c(
  paste0("script=12_Fig2_SignedEffects_Summary_Final"),
  paste0("run_stamp=", RUN_STAMP),
  paste0("corpus_tag=", CORPUS_TAG),
  paste0("dic_tag=", DIC_TAG),
  paste0("root=", ROOT),
  paste0("dic_file=", DIC_FILE)
), con = RUN_LOG)
capture.output(sessionInfo(), file = RUN_LOG, append = TRUE)
log_msg("WROTE:", RUN_LOG)


C_MAIN <- "AE-ILD"
TOP_A_FIG5A <- 24L
TOP_A_FIG5B <- 20L      # you can set 10/12/15
A_FORCE <- NULL
TOP_B_PER_A <- 6L
A_BREAKDOWN <- c("smoking","silica","male_sex","MUC5B","TERT","glucocorticoid")

f_sum  <- file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_rank <- file.path(DIR_TABLE, sprintf("abc_rankings_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_sum))  stop("Missing: ", f_sum, " (run 06)")
if (!file.exists(f_rank)) stop("Missing: ", f_rank, " (run 04)")

Ssum <- readr::read_csv(f_sum,  show_col_types = FALSE)
rank <- readr::read_csv(f_rank, show_col_types = FALSE)

if (!(C_MAIN %in% Ssum$C)) C_MAIN <- Ssum$C[1]

theme_fig5 <- ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face="bold", size=12),
    axis.title = ggplot2::element_text(size=10),
    axis.text.y = ggplot2::element_text(size=8),
    legend.title = ggplot2::element_text(size=9),
    legend.text = ggplot2::element_text(size=8),
    legend.position = "right"
  )

# -------------------------
# Fig. 2A
# -------------------------
A_all <- Ssum |>
  dplyr::filter(C == C_MAIN) |>
  dplyr::mutate(
    articles = as.integer(articles),
    balance = as.numeric(balance),
    balance_low = as.numeric(balance_low),
    balance_high = as.numeric(balance_high)
  )

A_dat <- if (!is.null(A_FORCE) && length(A_FORCE)) {
  A_all |>
    dplyr::filter(A %in% A_FORCE) |>
    dplyr::mutate(A_f = factor(A, levels = A_FORCE))
} else {
  n_take <- min(TOP_A_FIG5A, nrow(A_all))
  A_all |>
    dplyr::arrange(dplyr::desc(articles), balance) |>
    dplyr::slice_head(n = n_take) |>
    dplyr::arrange(balance, dplyr::desc(articles)) |>
    dplyr::mutate(A_f = factor(A, levels = A))
}

pA <- ggplot2::ggplot(A_dat, ggplot2::aes(x = balance, y = A_f)) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = balance_low, xmax = balance_high),
    height = 0, linewidth = 0.5, color = "grey40"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = articles, fill = balance),
    shape = 21, color = "grey30", alpha = 0.95, stroke = 0.3
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  ggplot2::scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish
  ) +
  ggplot2::scale_size_area(max_size = 8) +
  ggplot2::labs(
    title = "AE-ILD risk (signed effects)",
    x = "Balance (risk_down \u2190 0 \u2192 risk_up)",
    y = NULL, fill = "Balance", size = "Articles"
  ) +
  theme_fig5

# -------------------------
# Fig. 2B (FIXED)
# -------------------------
support_tab <- rank |>
  dplyr::filter(C == C_MAIN) |>
  dplyr::mutate(score_q = as.numeric(score_q)) |>
  dplyr::group_by(A) |>
  dplyr::summarise(
    sum_score_q = sum(pmax(score_q, 0), na.rm = TRUE),
    triad_support = log10(pmax(sum_score_q, 1e-9)),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(triad_support))

topA_support <- support_tab |>
  dplyr::slice_head(n = min(TOP_A_FIG5B, nrow(support_tab))) |>
  dplyr::pull(A)

sum_for_B <- Ssum |>
  dplyr::filter(C == C_MAIN) |>
  dplyr::select(A, articles, balance) |>
  dplyr::mutate(
    articles = as.integer(articles),
    balance = as.numeric(balance)
  )

B_dat <- support_tab |>
  dplyr::filter(A %in% topA_support) |>
  dplyr::inner_join(sum_for_B, by="A") |>
  dplyr::mutate(A_f = factor(A, levels = rev(unique(A))))


pB <- ggplot2::ggplot(B_dat, ggplot2::aes(x = triad_support, y = A_f)) +
  ggplot2::geom_point(
    ggplot2::aes(size = articles, fill = balance),
    shape = 21, color = "grey30", alpha = 0.95, stroke = 0.3
  ) +
  ggplot2::scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish
  ) +
  ggplot2::scale_size_area(max_size = 8) +
  ggplot2::labs(
    title = "Triad-based AE-ILD support",
    x = "Triad-based AE-ILD LBD support (log10 \u03a3 score_q)",
    y = NULL, fill = "Balance", size = "Articles"
  ) +
  theme_fig5

# -------------------------
# Fig. 2C (6 panels, 2 rows)
# -------------------------
rank_bd <- rank |>
  dplyr::filter(C == C_MAIN, A %in% A_BREAKDOWN) |>
  dplyr::mutate(
    AB_n11 = as.integer(AB_n11),
    BC_n11 = as.integer(BC_n11),
    cooc = AB_n11 + BC_n11,
    ABBC_strength = AB_n11 + BC_n11
  ) |>
  dplyr::filter(cooc > 0) |>
  dplyr::group_by(A) |>
  dplyr::slice_max(order_by = cooc, n = TOP_B_PER_A, with_ties = FALSE) |>
  dplyr::ungroup()

reorder_within <- function(x, by, within, sep="___"){
  newx <- paste(x, within, sep=sep)
  ord <- tapply(by, newx, function(z) -mean(z, na.rm=TRUE))
  factor(newx, levels = names(sort(ord)))
}
scale_y_reordered <- function(sep="___"){
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0(sep, ".+$"), "", x))
}

rank_bd <- rank_bd |>
  dplyr::mutate(B_f = reorder_within(B, by = cooc, within = A))

pC <- ggplot2::ggplot(rank_bd, ggplot2::aes(x = cooc, y = B_f)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = cooc, yend = B_f),
    linewidth = 0.8, color = "steelblue", alpha = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = cooc, fill = ABBC_strength),
    shape = 21, color = "grey30", alpha = 0.95, stroke = 0.3
  ) +
  ggplot2::facet_wrap(~A, nrow = 2, scales = "free_y") +
  ggplot2::scale_fill_gradient(low = "grey90", high = "steelblue") +
  ggplot2::scale_size_area(max_size = 9) +
  scale_y_reordered() +
  ggplot2::labs(
    title = "Intermediate B pathways to AE-ILD",
    x = "Co-occurrence count [#(A\u2227B) + #(B\u2227C)] (C = AE-ILD)",
    y = "Intermediate B",
    fill = "ABBC_strength",
    size = "Bubble size:\n#(A\u2227B)+#(B\u2227C)"
  ) +
  theme_fig5

fig5 <- (pA + pB) / pC +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(heights = c(1, 1.15))

save_dual <- function(stub, plot, w, h, dpi=350){
  if (!dir.exists(dirname(stub))) dir.create(dirname(stub), recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(paste0(stub, ".pdf"), plot, width = w, height = h, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stub, ".png"), plot, width = w, height = h, dpi = dpi)
  log_msg("WROTE: ", paste0(stub, ".pdf"))
  log_msg("WROTE: ", paste0(stub, ".png"))
}

stub <- file.path(DIR_FIG2, sprintf("Figure2_%s__%s", CORPUS_TAG, DIC_TAG))
ggplot2::ggsave(paste0(stub, ".pdf"), fig5, width = 12, height = 7.8, device = grDevices::cairo_pdf)
ggplot2::ggsave(paste0(stub, ".png"), fig5, width = 12, height = 7.8, dpi = 350)

log_msg("WROTE: ", paste0(stub, ".pdf"))
log_msg("WROTE: ", paste0(stub, ".png"))