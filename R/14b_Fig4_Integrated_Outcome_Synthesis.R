# =========================================================
# 14_Fig4_Integrated_Outcome_Synthesis.R
# Refined manuscript-facing outcome-specific hot spot map (Figure 4)
# - tighter layout
# - equal-width cards
# - stronger emphasis on hot spots
# - shorter biomarker/exploratory layer
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("readr","dplyr","stringr","ggplot2","patchwork","grid","scales"))
library(grid)

# Optional reads to keep linkage with analysis outputs
f_cooc <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_sign <- file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_biom <- file.path(DIR_TABLE, sprintf("biomarker_metrics_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_drug <- file.path(DIR_TABLE, sprintf("nonreview_cooc_drug_outcomes_full_%s__%s.csv", CORPUS_TAG, DIC_TAG))

if (file.exists(f_cooc)) cooc <- readr::read_csv(f_cooc, show_col_types = FALSE)
if (file.exists(f_sign)) sign <- readr::read_csv(f_sign, show_col_types = FALSE)
if (file.exists(f_biom)) biom <- readr::read_csv(f_biom, show_col_types = FALSE)
if (file.exists(f_drug)) drug <- readr::read_csv(f_drug, show_col_types = FALSE)

cards <- tibble::tribble(
  ~outcome,       ~theme,                           ~hotspots,                                                                            ~lower,                                              ~fill,
  "AE-ILD",       "Fibrotic injury phenotype",     "Honeycombing\nTraction bronchiectasis\nUIP / CPFE context\nSmoking",              "KL-6 / SP-D\nPeriostin / Osteopontin",               "#C74B50",
  "Progression",  "Progressive fibrosing axis",    "Fibrosis extent\nFVC decline\nMale sex / ACPA-RF\nUIP-like remodelling",          "KL-6 / SP-D\nCXCL13 / Th17",                         "#E08E2B",
  "Mortality",    "Severity + fibrosis burden",    "CPI / GAP\nUIP\nMale sex / Smoking\nACPA-RF",                                      "KL-6 / MMP7\nCCL18 / YKL40",                         "#4E79A7"
)

make_card <- function(title, theme_txt, hotspots_txt, lower_txt, header_fill){
  ggplot2::ggplot() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +

    # Header
    annotate("label",
             x = 0.5, y = 0.955, label = title,
             size = 5.2, fontface = "bold",
             fill = header_fill, colour = "white",
             label.size = 0, label.padding = unit(0.38, "lines")) +

    # Theme box
    annotate("label",
             x = 0.5, y = 0.79,
             label = paste0("Dominant theme\n", theme_txt),
             size = 4.0,
             fill = scales::alpha(header_fill, 0.10), colour = "black",
             label.size = 0.25, label.padding = unit(0.34, "lines")) +

    # Hot spots box (slightly larger, visual main actor)
    annotate("label",
             x = 0.5, y = 0.49,
             label = paste0("Hot spots\n", hotspots_txt),
             size = 4.0, fontface = "plain",
             fill = "white", colour = "black",
             label.size = 0.35, label.padding = unit(0.45, "lines")) +

    # Lower box (shortened text)
    annotate("label",
             x = 0.5, y = 0.18,
             label = paste0("Biomarker / exploratory layer\n", lower_txt),
             size = 3.5,
             fill = scales::alpha(header_fill, 0.08), colour = "black",
             label.size = 0.25, label.padding = unit(0.32, "lines"))
}

plots <- lapply(seq_len(nrow(cards)), function(i){
  make_card(cards$outcome[i], cards$theme[i], cards$hotspots[i], cards$lower[i], cards$fill[i])
})

p <- plots[[1]] + plots[[2]] + plots[[3]] +
  patchwork::plot_layout(widths = c(1,1,1)) +
  patchwork::plot_annotation(
    title = "Outcome-specific hot spot map in RA-ILD",
    subtitle = "Interpretation-focused summary of outcome-centred signals from the LBD pipeline",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, margin = margin(b = 2)),
      plot.margin = margin(4, 6, 2, 6)
    )
  )

ggplot2::ggsave(
  file.path(if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG")) DIR_FIG else getwd(), "Figure4_outcome_specific_hotspot_map_v2.png"),
  p, width = 12.8, height = 5.6, dpi = 320, bg = "white"
)

ggplot2::ggsave(
  file.path(if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG")) DIR_FIG else getwd(), "Figure4_outcome_specific_hotspot_map_v2.pdf"),
  p, width = 12.8, height = 5.6, bg = "white"
)

readr::write_csv(cards, file.path(DIR_TABLE, sprintf("outcome_specific_hotspot_map_labels_v2_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

message("Saved refined figure:")
message(" - ", file.path(if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG")) DIR_FIG else getwd(), "Figure4_outcome_specific_hotspot_map_v2.png"))
message(" - ", file.path(if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG")) DIR_FIG else getwd(), "Figure4_outcome_specific_hotspot_map_v2.pdf"))
