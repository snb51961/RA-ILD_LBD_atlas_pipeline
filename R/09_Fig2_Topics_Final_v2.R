# =========================================================
# 09_Fig2_Topics_Final.R
#  Figure 2. Topic structure and temporal evolution of RA-ILD literature
#   - Fig 2A: LDA topic–word distributions (K=8)  [Top terms panel]
#   - Fig 2B: Topic map (t-SNE of document-topic mixtures; shape=year_bin)
#
#  Input:
#   - data_proc/articles_YYYYMMDD.csv (from 01_fetch_pubmed.R)
#   - 00_setup_Final.R must exist (or adjust the source line) in project root
#
#  Output:
#   - {DIR_FIG2}/Fig2A_lda_topics_clean_{CORPUS_TAG}__{DIC_TAG}.pdf/png
#   - {DIR_FIG2}/Fig2B_topicmap_tsne_{CORPUS_TAG}__{DIC_TAG}.pdf/png
#   - {DIR_FIG2}/Fig2C_topic_year_counts_{CORPUS_TAG}__{DIC_TAG}.pdf/png
#   - {DIR_TABLE}/Fig2_topic_topwords_{CORPUS_TAG}__{DIC_TAG}.csv
#   - {DIR_TABLE}/Fig2_topic_labels_{CORPUS_TAG}__{DIC_TAG}.csv
#   - {DIR_TABLE}/Fig2C_topic_year_counts_{CORPUS_TAG}__{DIC_TAG}.csv
# =========================================================

source("00_setup_Final.R")

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c(
  "dplyr","readr","stringr","tibble","tidyr",
  "quanteda","topicmodels","Rtsne",
  "tidytext","ggplot2","ggrepel",
  "systemfonts","showtext"
))

# ------------------ Font (PDF safety) ------------------
cands  <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic","IPAexGothic","Noto Sans CJK JP")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

# ------------------ Helpers ------------------
latest_in_dir <- function(pattern, dir){
  xs <- Sys.glob(file.path(dir, pattern))
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

# Prefer today's legacy articles; else latest in data_proc
f_articles_csv <- file.path(DIR_PROC, sprintf("articles_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (!file.exists(f_articles_csv)) {
  f_articles_csv <- latest_in_dir("articles_*.csv", DIR_PROC)
}
if (is.na(f_articles_csv) || !file.exists(f_articles_csv)) {
  stop("articles_*.csv not found in data_proc. Run 01_fetch_pubmed.R first.")
}

articles <- readr::read_csv(f_articles_csv, show_col_types = FALSE) %>%
  dplyr::mutate(
    text = dplyr::if_else(!is.na(abstract) & !is.na(title), paste(title, abstract, sep=" "), ""),
    text = dplyr::if_else(is.na(text), "", text),
    year = suppressWarnings(as.integer(year))
  )

if (nrow(articles) < 30) stop("Documents too few (<30). Increase corpus or relax filters.")

# =========================================================
# Figure 2A: LDA topic–word distributions (K=8)
# =========================================================
set.seed(42)
K <- 8

# ---- Tokenize & DFM (domain stopwords consistent with your previous code) ----
sw_general <- quanteda::stopwords("en")
sw_domain_base <- c(
  "rheumatoid","arthritis","interstitial","lung","disease","ild","ip","pulmonary","fibrosis",
  "patients","patient","study","studies","review","case","cases","report","reports",
  "introduction","conclusion","methods","background","objective","aim","result","results",
  "purpose","mg","ml","day","days","week","weeks","year","years"
)

corp <- quanteda::corpus(articles, text_field = "text")
toks <- corp %>%
  quanteda::tokens(remove_punct=TRUE, remove_numbers=TRUE, remove_symbols=TRUE) %>%
  quanteda::tokens_tolower() %>%
  quanteda::tokens_remove(c(sw_general, sw_domain_base))

dfm_raw <- quanteda::dfm(toks)

# ---- staged trimming (robust) ----
dfm_trimmed <- quanteda::dfm_trim(dfm_raw, min_docfreq = 5, docfreq_type = "count")
dfm_trimmed <- quanteda::dfm_trim(dfm_trimmed, max_docfreq = 0.5, docfreq_type = "prop")

# relax if too small
if (quanteda::nfeat(dfm_trimmed) < 100) {
  dfm_trimmed <- quanteda::dfm_trim(dfm_raw, min_docfreq = 3, docfreq_type = "count")
  dfm_trimmed <- quanteda::dfm_trim(dfm_trimmed, max_docfreq = 0.8, docfreq_type = "prop")
}

stopifnot(quanteda::nfeat(dfm_trimmed) > 0, quanteda::ndoc(dfm_trimmed) > 0)

# ---- LDA ----
dtm <- quanteda::convert(dfm_trimmed, to = "topicmodels")
lda_fit <- topicmodels::LDA(
  dtm, k = K, method = "Gibbs",
  control = list(burnin = 500, iter = 2000, thin = 50, seed = 42)
)

# ---- beta ----
beta <- tidytext::tidy(lda_fit, matrix = "beta")  # topic, term, beta

# remove noisy tokens for plotting (based on your previous plot_stop idea)
plot_stop <- unique(c(
  tidytext::stop_words$word,
  "also","may","might","however","therefore","thus","suggest","suggests","suggested",
  "significant","significantly","associated","association","results","result",
  "use","used","using","methods","method","group","groups","level","levels",
  "clinical","evidence","management","systemic","diseases","connective","tissue",
  "month","months","ci","risk","factors","activity","primary","syndrome","pneumonia"
))

top_terms <- beta %>%
  dplyr::filter(!term %in% plot_stop) %>%
  dplyr::group_by(topic) %>%
  dplyr::slice_max(order_by = beta, n = 12, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(term_plot = tidytext::reorder_within(term, beta, topic))

# ---- save topwords table (for manuscript / supplement) ----
topic_topwords <- top_terms %>%
  group_by(topic) %>%
  summarise(top_words = paste(term[order(-beta)][1:min(10, n())], collapse=", "), .groups="drop") %>%
  arrange(topic)

f_topwords <- file.path(DIR_TABLE, sprintf("Fig2_topic_topwords_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(topic_topwords, f_topwords)
log_msg("WROTE:", f_topwords)

# Topic labels: do NOT use manual names (avoid mismatch after re-corpus)
facet_labels <- paste0("Topic ", 1:K)

# If K != 8, fallback to Topic 1..K
facet_labels <- if (length(manual_labels) == K) manual_labels else paste0("Topic ", 1:K)


# Topic color palette (consistent across Fig 2A/2B/2C)
topic_levels <- 1:K
topic_colors <- setNames(scales::hue_pal(l = 65, c = 100)(K), as.character(topic_levels))


lda_labeled <- top_terms %>%
  mutate(facet = factor(topic, levels = 1:K, labels = facet_labels),
         term_plot = tidytext::reorder_within(term, beta, facet))

p2a <- ggplot(lda_labeled, aes(x = beta, y = term_plot)) +
  geom_col(aes(fill = factor(topic)), color = "black", linewidth = 0.25) +
  scale_fill_manual(values = topic_colors, guide = "none") +
  tidytext::scale_y_reordered() +
  facet_wrap(~ facet, scales = "free_y", ncol = 2) +
  labs(
    title = "Figure 2A. LDA topic–word distributions (K=8)",
    x = "β (term | topic)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = jpfont),
        strip.text = element_text(face = "bold"))

f2a_pdf <- file.path(DIR_FIG2, sprintf("Fig2A_lda_topics_clean_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f2a_png <- file.path(DIR_FIG2, sprintf("Fig2A_lda_topics_clean_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggplot2::ggsave(f2a_pdf, p2a, device = cairo_pdf_device, width = 10, height = 8)
ggplot2::ggsave(f2a_png, p2a, width = 10, height = 8, dpi = 300)
log_msg("WROTE:", f2a_pdf)

# =========================================================
# Figure 2B: Topic map (t-SNE of document–topic mixtures)
# =========================================================
gamma <- topicmodels::posterior(lda_fit)$topics  # D x K
stopifnot(nrow(gamma) == quanteda::ndoc(dfm_trimmed))

# align doc_meta (safe)
doc_ids <- quanteda::docnames(dfm_trimmed)
idx <- suppressWarnings(as.integer(doc_ids))
if (any(is.na(idx)) || length(idx) != nrow(gamma)) idx <- seq_len(nrow(gamma))

doc_meta <- articles[idx, c("pmid","year","journal","title")]
year_bin <- cut(
  doc_meta$year,
  breaks = c(-Inf, 2004, 2009, 2014, 2019, Inf),
  labels = c("<=2004","2005-09","2010-14","2015-19","2020-"),
  right = TRUE
)

set.seed(42)
perp <- max(5, floor(nrow(gamma)/50))
tsne_res <- Rtsne::Rtsne(gamma, perplexity = perp, theta = 0.5,
                         check_duplicates = FALSE, pca = TRUE, max_iter = 1000)

tsne_df <- tibble(
  x = tsne_res$Y[,1],
  y = tsne_res$Y[,2],
  year = doc_meta$year,
  year_bin = year_bin,
  pmid = doc_meta$pmid,
  top_topic = max.col(gamma, ties.method = "first")
)

# label centers (median of each dominant-topic cluster)
centers <- tsne_df %>%
  group_by(top_topic) %>%
  summarise(cx = median(x), cy = median(y), n = n(), .groups = "drop") %>%
  mutate(label = facet_labels[top_topic])

# save labels
label_tbl <- tibble::tibble(topic = 1:K, label = facet_labels)
f_labels <- file.path(DIR_TABLE, sprintf("Fig2_topic_labels_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(label_tbl, f_labels)
log_msg("WROTE:", f_labels)

p2b <- ggplot(tsne_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(top_topic), shape = year_bin), size = 2, alpha = 0.85) +
  scale_color_manual(values = topic_colors) +
  ggrepel::geom_label_repel(
    data = centers,
    aes(x = cx, y = cy, label = label),
    size = 3.2, fontface = "bold", fill = "white", label.size = 0.2,
    box.padding = 0.5, point.padding = 0.5, segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(
    title = "Figure 2B. Topic map (t-SNE of document–topic mixtures)",
    x = "t-SNE 1", y = "t-SNE 2",
    color = "Dominant topic", shape = "Year bin"
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = jpfont))

f2b_pdf <- file.path(DIR_FIG2, sprintf("Fig2B_topicmap_tsne_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f2b_png <- file.path(DIR_FIG2, sprintf("Fig2B_topicmap_tsne_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggplot2::ggsave(f2b_pdf, p2b, device = cairo_pdf_device, width = 10, height = 8)
ggplot2::ggsave(f2b_png, p2b, width = 10, height = 8, dpi = 300)
log_msg("WROTE:", f2b_pdf)



# -------------------------
# Fig 2C: Yearly document counts by dominant topic
#  - Dominant topic for each document is defined as argmax(gamma)
# -------------------------
dominant_topic <- max.col(gamma, ties.method = "first")
doc_year <- suppressWarnings(as.integer(doc_meta$year))

topic_year_counts <- tibble::tibble(
  year = doc_year,
  topic = dominant_topic
) |>
  dplyr::filter(!is.na(year)) |>
  dplyr::count(year, topic, name = "n") |>
  tidyr::complete(
    year = seq(min(year, na.rm = TRUE), max(year, na.rm = TRUE)),
    topic = sort(unique(topic)),
    fill = list(n = 0)
  ) |>
  dplyr::arrange(topic, year)

p2c <- ggplot2::ggplot(topic_year_counts, ggplot2::aes(x = year, y = n, color = factor(topic))) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::geom_point(size = 0.6) +
  ggplot2::scale_color_manual(values = topic_colors, guide = "none") +
  ggplot2::facet_wrap(~ paste0("T", topic), ncol = 4, scales = "free_y") +
  ggplot2::labs(
    title = "Figure 2C. Yearly number of documents by dominant topic",
    x = "Year", y = "Number of documents"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    text = ggplot2::element_text(family = jpfont),
    strip.text = ggplot2::element_text(face = "bold")
  )

f2c_pdf <- file.path(DIR_FIG2, sprintf("Fig2C_topic_year_counts_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f2c_png <- file.path(DIR_FIG2, sprintf("Fig2C_topic_year_counts_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggplot2::ggsave(f2c_pdf, p2c, device = cairo_pdf_device, width = 12, height = 7)
ggplot2::ggsave(f2c_png, p2c, width = 12, height = 7, dpi = 300)
log_msg("WROTE:", f2c_pdf)


# ===============================
# Combined Figure 2 layout (A left-top, B right-top, C bottom spanning)
# ===============================
p2_combined <- (p2a + p2b) / p2c +
  patchwork::plot_layout(heights = c(1, 0.85))

f2_pdf <- file.path(DIR_FIG2, sprintf("Figure2_combined_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f2_png <- file.path(DIR_FIG2, sprintf("Figure2_combined_%s__%s.png", CORPUS_TAG, DIC_TAG))
ggplot2::ggsave(f2_pdf, p2_combined, device = cairo_pdf_device, width = 12.5, height = 10)
ggplot2::ggsave(f2_png, p2_combined, width = 12.5, height = 10, dpi = 300)
log_msg("WROTE:", f2_pdf)


f2c_csv <- file.path(DIR_TABLE, sprintf("Fig2C_topic_year_counts_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(topic_year_counts, f2c_csv)
log_msg("WROTE:", f2c_csv)

log_msg("=== DONE 09_Fig2_Topics_Final ===")


