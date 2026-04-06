# =========================================================
# 16_SuppFigS4_CaseReport_ABC_Only.R
# ---------------------------------------------------------
# Supplementary Figure S4: case-report ABC-only figure for AE-ILD in RA-ILD
#
# PURPOSE
#   Build a supplementary manuscript-facing case-report figure using ONLY the
#   already-generated case_report_abc_* outputs:
#     - Top case-report ABC triads (barplot)
#     - Case-report ABC network
#
# DESIGN
#   - No PubMed retrieval
#   - No re-running signed effects / biomarker atlas
#   - Reads existing case_report_abc_rankings + case_report_abc_edges
#   - Follows the visual logic of original SuppFigS3 as closely as practical
#   - Does NOT auto-run on source(); call the function explicitly
#
# USAGE
#   source("16_SuppFigS4_CaseReport_ABC_Only.R")
#   res <- run_case_report_suppfigS4_abc_only_figure()
# =========================================================

run_case_report_suppfigS4_abc_only_figure <- function(
  setup_file = NULL,
  root_dir   = NULL,
  verbose    = TRUE
) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  quiet_install <- function(pkgs) {
    to_install <- setdiff(pkgs, rownames(installed.packages()))
    if (length(to_install)) install.packages(to_install, dependencies = TRUE)
    invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
  }

  quiet_install(c("readr","dplyr","stringr","tibble","forcats","ggplot2","patchwork","igraph","ggraph"))

  coerce_scalar_character <- function(x) {
    if (inherits(x, "fs_path")) x <- as.character(x)
    if (is.factor(x)) x <- as.character(x)
    if (is.list(x) && length(x) == 1) x <- unlist(x, use.names = FALSE)
    if (is.null(x) || length(x) == 0) return(NULL)
    if (!is.atomic(x)) return(NULL)
    if (!is.character(x)) x <- tryCatch(as.character(x), error = function(e) NULL)
    if (is.null(x) || length(x) == 0) return(NULL)
    x <- x[[1]]
    if (is.na(x) || !nzchar(x)) return(NULL)
    x
  }

  normalize_path_safe <- function(path) {
    path <- coerce_scalar_character(path)
    if (is.null(path)) return(NULL)
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  get_scalar_from_env <- function(name, default = NULL, inherits = TRUE) {
    if (!exists(name, inherits = inherits)) return(default)
    coerce_scalar_character(get(name, inherits = inherits)) %||% default
  }

  find_first_existing_dir <- function(paths) {
    paths <- unique(Filter(Negate(is.null), lapply(paths, normalize_path_safe)))
    for (p in paths) if (dir.exists(p)) return(p)
    NULL
  }

  find_latest_file <- function(dir, regex, recursive = FALSE) {
    dir <- normalize_path_safe(dir)
    if (is.null(dir) || !dir.exists(dir)) return(NA_character_)
    xs <- list.files(dir, pattern = regex, full.names = TRUE, recursive = recursive)
    if (!length(xs)) return(NA_character_)
    xs[which.max(file.info(xs)$mtime)]
  }

  # -------------------------------------------------------
  # Load setup if available (00_setup.R or 00_setup_Final.R)
  # -------------------------------------------------------
  setup_candidates <- unique(Filter(Negate(is.null), c(
    coerce_scalar_character(setup_file),
    file.path(getwd(), "00_setup.R"),
    file.path(getwd(), "00_setup_Final.R"),
    file.path(get_scalar_from_env("RAILD_ROOT"), "00_setup.R"),
    file.path(get_scalar_from_env("RAILD_ROOT"), "00_setup_Final.R")
  )))

  setup_loaded <- FALSE
  for (sp in setup_candidates) {
    if (file.exists(sp)) {
      source(sp, local = FALSE)
      setup_loaded <- TRUE
      if (isTRUE(verbose)) message("Loaded setup: ", normalize_path_safe(sp))
      break
    }
  }

  if (!setup_loaded && isTRUE(verbose)) {
    message("No setup file was loaded. Using existing session objects and standard fallbacks.")
  }

  ROOT2 <- find_first_existing_dir(c(
    coerce_scalar_character(root_dir),
    get_scalar_from_env("ROOT"),
    get_scalar_from_env("RAILD_ROOT"),
    getwd(),
    path.expand("~")
  ))
  ROOT2 <- ROOT2 %||% normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  CORPUS_TAG2 <- get_scalar_from_env("CORPUS_TAG") %||% "pm_1980_20251231"
  DIC_TAG2    <- get_scalar_from_env("DIC_TAG") %||% "analysis_v1_genetics"

  DIR_TABLE2 <- normalize_path_safe(get_scalar_from_env("DIR_TABLE"))
  if (is.null(DIR_TABLE2)) DIR_TABLE2 <- file.path(ROOT2, "output", CORPUS_TAG2, DIC_TAG2, "table")

  DIR_FIG2 <- normalize_path_safe(get_scalar_from_env("DIR_FIG2"))
  if (is.null(DIR_FIG2)) DIR_FIG2 <- normalize_path_safe(get_scalar_from_env("DIR_FIG"))
  if (is.null(DIR_FIG2)) DIR_FIG2 <- file.path(ROOT2, "output", CORPUS_TAG2, DIC_TAG2, "fig")

  DIR_DIC2 <- normalize_path_safe(get_scalar_from_env("DIR_DIC"))
  if (is.null(DIR_DIC2)) DIR_DIC2 <- file.path(ROOT2, "dic")

  dir.create(DIR_TABLE2, recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_FIG2, recursive = TRUE, showWarnings = FALSE)

  # -------------------------------------------------------
  # Inputs
  # -------------------------------------------------------
  f_rank <- file.path(DIR_TABLE2, sprintf("case_report_abc_rankings_%s__%s.csv", CORPUS_TAG2, DIC_TAG2))
  if (!file.exists(f_rank)) {
    f_rank <- find_latest_file(DIR_TABLE2, "^case_report_abc_rankings_.*\\.csv$")
  }
  if (is.na(f_rank) || !file.exists(f_rank)) {
    stop("case_report_abc_rankings_*.csv not found under DIR_TABLE: ", DIR_TABLE2)
  }

  f_edges <- find_latest_file(DIR_TABLE2, "^case_report_abc_edges_top.*\\.csv$")
  if (!is.na(f_edges) && !file.exists(f_edges)) f_edges <- NA_character_

  dic_path <- file.path(DIR_DIC2, "ra_ild_dictionary_analysis_v1_genetics.csv")
  if (!file.exists(dic_path)) dic_path <- NA_character_

  if (isTRUE(verbose)) {
    message("Using rankings : ", f_rank)
    message("Using edges    : ", ifelse(is.na(f_edges), "<fallback from rankings>", f_edges))
    message("Using dictionary: ", ifelse(is.na(dic_path), "<not found>", dic_path))
    message("DIR_FIG2       : ", DIR_FIG2)
  }

  R0 <- readr::read_csv(f_rank, show_col_types = FALSE)
  if (!all(c("A","B","C") %in% names(R0))) stop("Rankings file must contain A, B, C columns.")

  score_col <- if ("score_q" %in% names(R0)) "score_q" else if ("score" %in% names(R0)) "score" else NULL
  if (is.null(score_col)) stop("No score column found in case_report_abc_rankings file.")

  E0 <- if (!is.na(f_edges) && file.exists(f_edges)) {
    readr::read_csv(f_edges, show_col_types = FALSE)
  } else {
    dplyr::bind_rows(
      R0 |> dplyr::transmute(from = A, to = B, w = if ("AB_n11" %in% names(R0)) AB_n11 else 1L, kind = "AB"),
      R0 |> dplyr::transmute(from = B, to = C, w = if ("BC_n11" %in% names(R0)) BC_n11 else 1L, kind = "BC")
    ) |>
      dplyr::group_by(from, to, kind) |>
      dplyr::summarise(w = max(w, na.rm = TRUE), .groups = "drop")
  }

  # -------------------------------------------------------
  # Parameters (easy to tweak)
  # -------------------------------------------------------
  TOP_TRIADS <- 16L
  CAP_PER_A  <- NA_integer_   # set e.g. 2L to diversify A terms
  NET_LAYOUT <- "fr"
  FIG_STUB_BASENAME <- "SupplementaryFigureS4_case_report_ABC_only"
  FIG_WIDTH  <- 12
  FIG_HEIGHT <- 10

  # -------------------------------------------------------
  # Dictionary-based role mapping (optional but preferred)
  # -------------------------------------------------------
  ROLE_A_CLASSES   <- c("drug","bio","gene","exposure","host","event","nonpharm","vaccine","population","system","strategy")
  ROLE_B_CLASSES   <- c("biomarker","molecular","microbiome","cell","activity","pft","phenotype","trend","complication")
  ROLE_BPR_CLASSES <- c("pattern","radiology","qct","airway")
  ROLE_C_CLASSES   <- c("outcome")

  class_lookup <- function(term_vec, dic = NULL) {
    if (is.null(dic)) return(rep(NA_character_, length(term_vec)))
    out <- dic$class[match(term_vec, dic$term)]
    ifelse(is.na(out), NA_character_, out)
  }

  role_of_class <- function(cls) {
    dplyr::case_when(
      cls %in% ROLE_C_CLASSES   ~ "C",
      cls %in% ROLE_BPR_CLASSES ~ "Bprime",
      cls %in% ROLE_B_CLASSES   ~ "B",
      cls %in% ROLE_A_CLASSES   ~ "A",
      TRUE ~ "Other"
    )
  }

  dic <- NULL
  if (!is.na(dic_path) && file.exists(dic_path)) {
    dic <- readr::read_csv(dic_path, show_col_types = FALSE) |>
      dplyr::mutate(
        term = as.character(term),
        class = tolower(trimws(as.character(class)))
      )
  }

  # -------------------------------------------------------
  # Select top triads
  # -------------------------------------------------------
  R1 <- R0 |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(score_col), as.numeric)
    ) |>
    dplyr::arrange(dplyr::desc(.data[[score_col]]))

  if (!is.na(CAP_PER_A)) {
    R1 <- R1 |>
      dplyr::group_by(A) |>
      dplyr::slice_head(n = CAP_PER_A) |>
      dplyr::ungroup() |>
      dplyr::arrange(dplyr::desc(.data[[score_col]]))
  }

  Rf <- R1 |> dplyr::slice_head(n = min(TOP_TRIADS, nrow(R1)))
  if (nrow(Rf) == 0) stop("No triads available in rankings file.")

  C_main <- Rf$C[1]

  # -------------------------------------------------------
  # Panel A: top triads barplot
  # -------------------------------------------------------
  B_role <- if (!is.null(dic)) {
    role_of_class(class_lookup(Rf$B, dic))
  } else {
    rep("B", nrow(Rf))
  }

  R_bar <- Rf |>
    dplyr::mutate(
      B_disp = ifelse(B_role == "Bprime", paste0(B, " [B']"), B),
      triad_label = paste0(A, " -> ", B_disp, " -> ", C),
      triad_label = forcats::fct_reorder(triad_label, .data[[score_col]])
    )

  pA <- ggplot2::ggplot(R_bar, ggplot2::aes(x = .data[[score_col]], y = triad_label)) +
    ggplot2::geom_col(fill = "#4E79A7") +
    ggplot2::labs(
      title = "Top case-report ABC triads (AE-ILD)",
      x = "Triad support score (score_q)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # -------------------------------------------------------
  # Panel B: ABC network
  # -------------------------------------------------------
  keep_nodes <- unique(c(Rf$A, Rf$B, Rf$C))

  E <- E0 |>
    dplyr::filter(from %in% keep_nodes, to %in% keep_nodes) |>
    dplyr::mutate(
      w = as.numeric(w),
      kind = as.character(kind)
    ) |>
    dplyr::filter(!is.na(from), from != "", !is.na(to), to != "")

  if (nrow(E) == 0) {
    pB <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No edges available for the selected case-report triads.", size = 4) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = "Case-report ABC network (AE-ILD)") +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  } else {
    nodes <- tibble::tibble(name = unique(c(E$from, E$to))) |>
      dplyr::mutate(
        class = class_lookup(name, dic),
        role  = role_of_class(class),
        role  = ifelse(name %in% Rf$A, "A", role),
        role  = ifelse(name %in% Rf$B & !(name %in% Rf$A), ifelse(role == "Bprime", "Bprime", "B"), role),
        role  = ifelse(name %in% Rf$C, "C", role),
        role2 = ifelse(is.na(role) | role == "Other", "Other", role)
      )

    g <- igraph::graph_from_data_frame(E, directed = TRUE, vertices = nodes)
    igraph::E(g)$kind   <- E$kind
    igraph::E(g)$weight <- as.numeric(E$w)
    igraph::E(g)$w_plot <- log1p(pmax(igraph::E(g)$weight, 1))
    igraph::V(g)$deg    <- igraph::degree(g, mode = "all")
    igraph::V(g)$label  <- igraph::V(g)$name

    col_edge_AB <- "#C74B50"
    col_edge_BC <- "#4E79A7"

    shape_map <- c(A = 22, B = 21, Bprime = 24, C = 23, Other = 25)
    fill_map  <- c(A = "#1b4f72", B = "#117a65", Bprime = "#d35400", C = "#b03a2e", Other = "#7d7d7d")

    set.seed(42)
    pB <- ggraph::ggraph(g, layout = NET_LAYOUT) +
      ggraph::geom_edge_link(
        ggplot2::aes(edge_colour = kind, edge_width = w_plot, linetype = kind),
        alpha = 0.65,
        lineend = "round"
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(shape = role2, size = deg, fill = role2),
        colour = "grey20", stroke = 0.7
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = label),
        size = 3.2,
        repel = TRUE
      ) +
      ggraph::scale_edge_colour_manual(
        values = c(AB = col_edge_AB, BC = col_edge_BC),
        labels = c(AB = "AB (A -> B/B')", BC = paste0("BC (B/B' -> ", C_main, ")")),
        name = "Edge type"
      ) +
      ggraph::scale_edge_linetype_manual(
        values = c(AB = "solid", BC = "longdash"),
        name = "Edge type"
      ) +
      ggraph::scale_edge_width_continuous(
        range = c(0.45, 2.8),
        name = "Evidence (log1p count)"
      ) +
      ggplot2::scale_shape_manual(values = shape_map, name = "Node role") +
      ggplot2::scale_fill_manual(values = fill_map, guide = "none") +
      ggplot2::scale_size_continuous(range = c(3, 9), guide = "none") +
      ggplot2::labs(
        title = "Case-report ABC network (AE-ILD)"
      ) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        legend.position = "bottom",
        panel.grid = ggplot2::element_blank()
      )
  }

  # -------------------------------------------------------
  # Assemble and save
  # -------------------------------------------------------
  fig <- (pA / pB) +
    patchwork::plot_layout(heights = c(0.78, 1.22)) +
    patchwork::plot_annotation(
      title = "Supplementary Figure S4. Exploratory case-report ABC bridge map for AE-ILD in RA-ILD",
      subtitle = "Case reports / case series only: top ABC triads and network view",
      tag_levels = "A"
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.99)
    )

  fig_stub <- file.path(DIR_FIG2, sprintf("%s_%s__%s", FIG_STUB_BASENAME, CORPUS_TAG2, DIC_TAG2))
  ggplot2::ggsave(paste0(fig_stub, ".pdf"), fig, width = FIG_WIDTH, height = FIG_HEIGHT, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(fig_stub, ".png"), fig, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 350)

  f_top <- file.path(DIR_TABLE2, sprintf("case_report_abc_toptriads_for_SupplementaryFigureS4_%s__%s.csv", CORPUS_TAG2, DIC_TAG2))
  readr::write_csv(R_bar, f_top)

  if (isTRUE(verbose)) {
    message("WROTE: ", paste0(fig_stub, ".pdf"))
    message("WROTE: ", paste0(fig_stub, ".png"))
    message("WROTE: ", f_top)
  }

  invisible(list(
    rankings = R0,
    selected_triads = Rf,
    edges = E,
    figure = fig,
    files = list(
      rankings = f_rank,
      edges = f_edges,
      toptriads_csv = f_top,
      fig_pdf = paste0(fig_stub, ".pdf"),
      fig_png = paste0(fig_stub, ".png")
    ),
    resolved = list(
      ROOT = ROOT2,
      DIR_TABLE = DIR_TABLE2,
      DIR_FIG2 = DIR_FIG2,
      DIR_DIC = DIR_DIC2,
      CORPUS_TAG = CORPUS_TAG2,
      DIC_TAG = DIC_TAG2
    )
  ))
}

# backward-compatible alias
run_case_report_abc_only_figure <- run_case_report_suppfigS4_abc_only_figure

res <- run_case_report_suppfigS4_abc_only_figure(verbose = TRUE)
