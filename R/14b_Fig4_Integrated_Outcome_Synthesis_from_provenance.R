# =========================================================
# 14b_Fig4_Integrated_Outcome_Synthesis_updated_from_provenance.R
# Figure 4 manuscript-facing layout rebuilt from 14a provenance outputs
# - derives shared vs outcome-centred terms from 14a provenance
# - keeps manuscript-style compact box layout
# - shows † for terms present in 2/3 outcomes but kept in outcome-centred boxes
# =========================================================

.find_setup <- function() {
  cand <- c(
    file.path(getwd(), "00_setup_Final.R"),
    file.path(getwd(), "00_setup.R"),
    file.path(Sys.getenv("RAILD_ROOT", unset=""), "00_setup_Final.R"),
    file.path(Sys.getenv("RAILD_ROOT", unset=""), "00_setup.R")
  )
  for (p in cand) if (nzchar(p) && file.exists(p)) return(p)
  stop("Cannot find 00_setup_Final.R or 00_setup.R.")
}
source(.find_setup())

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}
quiet_install(c("readr","dplyr","stringr","tibble","ggplot2","grid","patchwork","scales"))

library(grid)

find_first_existing <- function(paths){
  ok <- paths[file.exists(paths)]
  if (length(ok)) ok[1] else NA_character_
}
find_latest_file <- function(dir, regex){
  if (is.na(dir) || !dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

canon_outcome <- function(x){
  xl <- stringr::str_to_lower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(xl, "^ae-ild$|acute") ~ "AE-ILD",
    stringr::str_detect(xl, "progress") ~ "Progression",
    stringr::str_detect(xl, "mort") ~ "Mortality",
    TRUE ~ as.character(x)
  )
}

clean_display <- function(x){
  x <- as.character(x)
  x <- dplyr::case_when(
    x == "ACPA RF high" ~ "ACPA/RF high",
    x == "SP A" ~ "SP-A",
    x == "GM CSF" ~ "GM-CSF",
    x == "IL8 CXCL8" ~ "IL-8/CXCL8",
    x == "WNT beta catenin" ~ "WNT/β-catenin",
    x == "drug induced ILD" ~ "drug-induced ILD",
    TRUE ~ x
  )
  x
}

collapse_terms <- function(x, per_line = 3){
  x <- stats::na.omit(as.character(x))
  if (!length(x)) return("—")
  grp <- ceiling(seq_along(x)/per_line)
  lines <- split(x, grp)
  paste(vapply(lines, function(z) paste(z, collapse = " / "), character(1)), collapse = "\n")
}

# ---------------------------------------------------------
# Inputs
# ---------------------------------------------------------
f_hot <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("Figure4_hotspot_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^Figure4_hotspot_layer_provenance_fully_datadriven_v5_.*\\.csv$")
))
f_bio <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("Figure4_biomarker_layer_provenance_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^Figure4_biomarker_layer_provenance_fully_datadriven_v5_.*\\.csv$")
))
f_lab <- find_first_existing(c(
  file.path(DIR_TABLE, sprintf("Figure4_outcome_specific_hotspot_map_labels_fully_datadriven_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)),
  find_latest_file(DIR_TABLE, "^Figure4_outcome_specific_hotspot_map_labels_fully_datadriven_v5_.*\\.csv$")
))

if (is.na(f_hot) || is.na(f_bio) || is.na(f_lab)) {
  stop("Missing one or more Figure 4 provenance inputs.")
}

HOT <- readr::read_csv(f_hot, show_col_types = FALSE) |>
  dplyr::mutate(
    outcome = canon_outcome(outcome),
    display = clean_display(display)
  )

BIO <- readr::read_csv(f_bio, show_col_types = FALSE) |>
  dplyr::mutate(
    outcome = canon_outcome(outcome),
    display = clean_display(display)
  )

LAB <- readr::read_csv(f_lab, show_col_types = FALSE) |>
  dplyr::mutate(outcome = canon_outcome(outcome))

OUTCOME_ORDER <- c("AE-ILD","Progression","Mortality")

CARD_FILL <- LAB |>
  dplyr::select(outcome, fill) |>
  dplyr::distinct()

# ---------------------------------------------------------
# Shared vs outcome-centred derivation
# ---------------------------------------------------------
derive_layer <- function(df, n_keep_outcome, shared_order = NULL){
  counts <- df |>
    dplyr::group_by(display) |>
    dplyr::summarise(
      n_outcomes = dplyr::n_distinct(outcome),
      order_mean = mean(selection_order, na.rm = TRUE),
      .groups = "drop"
    )

  df2 <- df |>
    dplyr::left_join(counts, by = "display") |>
    dplyr::mutate(
      is_shared = n_outcomes == 3L,
      has_dagger = n_outcomes == 2L
    )

  shared <- df2 |>
    dplyr::filter(is_shared) |>
    dplyr::group_by(display) |>
    dplyr::summarise(order_mean = min(order_mean, na.rm = TRUE), .groups = "drop")

  if (!is.null(shared_order)) {
    shared <- shared |>
      dplyr::mutate(order_custom = match(display, shared_order)) |>
      dplyr::arrange(is.na(order_custom), order_custom, order_mean, display)
  } else {
    shared <- shared |>
      dplyr::arrange(order_mean, display)
  }

  shared_terms <- shared$display

  outcome_centred <- df2 |>
    dplyr::filter(!is_shared) |>
    dplyr::group_by(outcome) |>
    dplyr::arrange(selection_order, .by_group = TRUE) |>
    dplyr::slice_head(n = n_keep_outcome) |>
    dplyr::ungroup() |>
    dplyr::mutate(display_marked = ifelse(has_dagger, paste0(display, "†"), display))

  list(shared_terms = shared_terms, outcome_terms = outcome_centred)
}

shared_hotspot_preferred <- c(
  "Male sex","Smoking","disease duration","drug-induced ILD",
  "Traction bronchiectasis","UIP","CPI index","Honeycombing"
)
shared_biomarker_preferred <- c(
  "KL-6","SP-D","ACPA/RF","MMP-7","MPO"
)

hot_layers <- derive_layer(HOT, n_keep_outcome = 8, shared_order = shared_hotspot_preferred)
bio_layers <- derive_layer(BIO, n_keep_outcome = 8, shared_order = shared_biomarker_preferred)

shared_hotspots_text <- collapse_terms(hot_layers$shared_terms, per_line = 4)
shared_biomarkers_text <- collapse_terms(bio_layers$shared_terms, per_line = 4)

hot_outcome_map <- hot_layers$outcome_terms |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(text = collapse_terms(display_marked, per_line = 2), .groups = "drop")

bio_outcome_map <- bio_layers$outcome_terms |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(text = collapse_terms(display_marked, per_line = 3), .groups = "drop")

get_text <- function(map_df, oc){
  z <- map_df$text[match(oc, map_df$outcome)]
  ifelse(is.na(z), "—", z)
}

# ---------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------
box_theme <- list(
  label.size = 0.35,
  label.padding = unit(0.35, "lines")
)

make_canvas <- function(){
  ggplot2::ggplot() +
    ggplot2::xlim(0, 100) +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = margin(10, 12, 18, 12))
}

annot_box <- function(p, x, y, label, fill = "white", colour = "black",
                      size = 3.2, fontface = "plain", box.size = 0.35,
                      padding = 0.32, hjust = 0.5, lineheight = 1.0){
  p + ggplot2::annotate(
    "label", x = x, y = y, label = label,
    fill = fill, colour = colour, size = size, fontface = fontface,
    label.size = box.size, label.padding = unit(padding, "lines"),
    hjust = hjust, lineheight = lineheight
  )
}

# ---------------------------------------------------------
# Build figure
# ---------------------------------------------------------
p <- make_canvas()

# row labels
p <- p +
  ggplot2::annotate("text", x = 8, y = 84, label = "Shared hot spots", hjust = 0, size = 5.0, fontface = "bold") +
  ggplot2::annotate("text", x = 8, y = 64, label = "Outcome-centred\nhot spots", hjust = 0, size = 5.0, fontface = "bold") +
  ggplot2::annotate("segment", x = 6, xend = 94, y = 49, yend = 49, colour = "grey80", linewidth = 0.4) +
  ggplot2::annotate("text", x = 8, y = 37, label = "Shared biomarkers", hjust = 0, size = 5.0, fontface = "bold") +
  ggplot2::annotate("text", x = 8, y = 17, label = "Outcome-centred\nbiomarkers /\nexploratory terms", hjust = 0, size = 5.0, fontface = "bold")

# shared headers and shared boxes
p <- annot_box(p, 65, 92, "Shared across AE-ILD, progression, and mortality",
               fill = "#55565B", colour = "white", size = 3.8, fontface = "bold",
               box.size = 0, padding = 0.26)
p <- annot_box(p, 65, 83, shared_hotspots_text, size = 3.2, padding = 0.34)

p <- annot_box(p, 65, 45, "Shared across AE-ILD, progression, and mortality",
               fill = "#55565B", colour = "white", size = 3.8, fontface = "bold",
               box.size = 0, padding = 0.26)
p <- annot_box(p, 65, 36, shared_biomarkers_text, size = 3.2, padding = 0.34)

# outcome headers
x_pos <- c("AE-ILD" = 41, "Progression" = 65, "Mortality" = 89)
for (oc in OUTCOME_ORDER) {
  fill_i <- CARD_FILL$fill[match(oc, CARD_FILL$outcome)]
  p <- annot_box(p, x_pos[[oc]], 72, oc, fill = fill_i, colour = "white",
                 size = 4.3, fontface = "bold", box.size = 0, padding = 0.26)
  p <- annot_box(p, x_pos[[oc]], 25, oc, fill = fill_i, colour = "white",
                 size = 4.3, fontface = "bold", box.size = 0, padding = 0.26)
}

# outcome-centred boxes
for (oc in OUTCOME_ORDER) {
  p <- annot_box(p, x_pos[[oc]], 61, get_text(hot_outcome_map, oc), size = 3.0, padding = 0.32)
  p <- annot_box(p, x_pos[[oc]], 14, get_text(bio_outcome_map, oc), size = 2.95, padding = 0.32)
}

# footnote
dagger_present <- any(stringr::str_detect(hot_outcome_map$text, "†"), na.rm = TRUE) ||
                  any(stringr::str_detect(bio_outcome_map$text, "†"), na.rm = TRUE)
if (dagger_present) {
  p <- p + ggplot2::annotate(
    "text", x = 50, y = 3.8,
    label = "† term present in 2 of 3 outcomes and therefore shown in the outcome-centred layer rather than the shared-across-all-3 layer",
    hjust = 0.5, size = 2.8, colour = "grey25"
  )
}

# ---------------------------------------------------------
# Outputs
# ---------------------------------------------------------
FIG_DIR <- if (exists("DIR_FIG4")) DIR_FIG4 else if (exists("DIR_FIG2")) DIR_FIG2 else if (exists("DIR_FIG")) DIR_FIG else getwd()
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)

stub <- file.path(FIG_DIR, sprintf("Figure4_outcome_specific_hotspot_map_v3_from_provenance_%s__%s", CORPUS_TAG, DIC_TAG))
f_pdf <- paste0(stub, ".pdf")
f_png <- paste0(stub, ".png")
f_cards <- file.path(DIR_TABLE, sprintf("Figure4_outcome_specific_hotspot_map_labels_v3_from_provenance_%s__%s.csv", CORPUS_TAG, DIC_TAG))

ggplot2::ggsave(f_pdf, p, width = 12.8, height = 6.4, device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(f_png, p, width = 12.8, height = 6.4, dpi = 320, bg = "white")

cards_out <- tibble::tibble(
  outcome = OUTCOME_ORDER,
  fill = CARD_FILL$fill[match(OUTCOME_ORDER, CARD_FILL$outcome)],
  shared_hotspots = shared_hotspots_text,
  outcome_hotspots = c(get_text(hot_outcome_map, "AE-ILD"), get_text(hot_outcome_map, "Progression"), get_text(hot_outcome_map, "Mortality")),
  shared_biomarkers = shared_biomarkers_text,
  outcome_biomarkers = c(get_text(bio_outcome_map, "AE-ILD"), get_text(bio_outcome_map, "Progression"), get_text(bio_outcome_map, "Mortality"))
)
readr::write_csv(cards_out, f_cards)

message("Saved updated Figure 4 manuscript-facing version:")
message(" - ", f_pdf)
message(" - ", f_png)
message(" - ", f_cards)
