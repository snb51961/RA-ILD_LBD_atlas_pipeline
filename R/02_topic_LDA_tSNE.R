############################################################
# 02) Topic Model (LDA) + t-SNE + Topic Trending
# ----------------------------------------------------------
#  - Replicates all topic-map figures used in the manuscript:
#       * lda_topics_*.pdf
#       * topicmap_*.pdf
#       * topicmap_labeled_manual_*.pdf
#       * topic_year_trends_*.pdf
#       * topic_year_heatmap_*.pdf
#       * topic_dominant_counts_*.pdf
#
#  - Fully faithful to "図の作成コードtopicmap.txt"
#  - ROOT uses here::here() for GitHub reproducibility
############################################################

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(quanteda)
  library(topicmodels)
  library(Rtsne)
  library(tidytext)
  library(ggrepel)
  library(systemfonts)
  library(showtext)
})

############################################################
# 0) PATHS
############################################################

ROOT     <- here::here()
DIR_PROC <- file.path(ROOT, "data_proc")
DIR_OUT  <- file.path(ROOT, "output")
DIR_FIG  <- file.path(ROOT, "fig_pub")

dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

############################################################
# 1) Load latest articles
############################################################

pick_latest <- function(pattern, dirs){
  xs <- unlist(lapply(dirs, function(d) Sys.glob(file.path(d, pattern))))
  if (!length(xs)) stop("articles_*.csv not found.")
  xs[which.max(file.info(xs)$mtime)]
}

f_articles <- pick_latest("articles_*.csv", c(DIR_PROC))
articles <- read_csv(f_articles, show_col_types = FALSE) %>%
  mutate(text = coalesce(text, paste(title, abstract, sep=" "))) %>%
  mutate(text = ifelse(is.na(text), "", text))

if (nrow(articles) < 30) stop("Too few documents (<30). LDA unsuitable.")

############################################################
# 2) Tokenization → DFM
############################################################

sw_general <- quanteda::stopwords("en")
sw_domain  <- c(
  "rheumatoid","arthritis","interstitial","lung","disease","ild","ip","pulmonary","fibrosis",
  "patients","patient","study","studies","review","case","cases","report","reports",
  "introduction","conclusion","methods","background","objective","aim","result","results",
  "purpose","mg","ml","day","days","week","weeks","year","years"
)

corp <- quanteda::corpus(articles, text_field = "text")
toks <- corp %>%
  tokens(remove_punct=TRUE, remove_numbers=TRUE, remove_symbols=TRUE) %>%
  tokens_tolower() %>%
  tokens_remove(c(sw_general, sw_domain))

dfm_raw <- dfm(toks)

# 2段階 trimming
dfm_trimmed <- dfm_trim(dfm_raw, min_docfreq=5, docfreq_type="count")
dfm_trimmed <- dfm_trim(dfm_trimmed, max_docfreq=0.5, docfreq_type="prop")

if (nfeat(dfm_trimmed) < 100){
  dfm_trimmed <- dfm_trim(dfm_raw, min_docfreq=3, docfreq_type="count")
  dfm_trimmed <- dfm_trim(dfm_trimmed, max_docfreq=0.8, docfreq_type="prop")
}

stopifnot(nfeat(dfm_trimmed)>0, ndoc(dfm_trimmed)>0)

############################################################
# 3) LDA (K=8, Gibbs)
############################################################

set.seed(42)
K <- 8
dtm <- quanteda::convert(dfm_trimmed, to="topicmodels")
lda_fit <- LDA(dtm, k=K, method="Gibbs",
               control=list(burnin=500, iter=2000, thin=50, seed=42))

############################################################
# 4) β (top terms per topic) → Fig2A
############################################################

beta <- tidy(lda_fit, matrix="beta")
TOP_N <- 12

top_terms <- beta %>%
  group_by(topic) %>%
  slice_max(beta, n = TOP_N, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(term_plot = tidytext::reorder_within(term, beta, topic))

# PDF
p_topics <- ggplot(top_terms, aes(beta, term_plot)) +
  geom_col() +
  tidytext::scale_y_reordered() +
  facet_wrap(~ paste0("Topic ", topic), scales="free_y", ncol=2) +
  labs(title = "LDA Top Terms", x="β", y=NULL) +
  theme_minimal(base_size=12)

f_topics <- file.path(DIR_FIG, paste0("lda_topics_", STAMP, ".pdf"))
ggsave(f_topics, p_topics, width=10, height=8)
message("Saved: ", f_topics)

# Save clean word table (再現性資料)
topic_words_table <- top_terms %>%
  group_by(topic) %>%
  summarise(top_words = paste(term[order(-beta)][1:10], collapse=", "),
            .groups="drop")
write_csv(topic_words_table, file.path(DIR_OUT, paste0("topic_topwords_", STAMP, ".csv")))

############################################################
# 5) γ (document-topic) → t-SNE → Fig2B
############################################################

gamma <- posterior(lda_fit)$topics
stopifnot(nrow(gamma) == ndoc(dfm_trimmed))

doc_ids <- docnames(dfm_trimmed)
idx <- suppressWarnings(as.integer(doc_ids))
if (any(is.na(idx)) || length(idx)!=nrow(gamma)) idx <- seq_len(nrow(gamma))

doc_meta <- articles[idx, c("pmid","year","journal","title")]

year_bin <- cut(doc_meta$year,
                breaks=c(-Inf,2004,2009,2014,2019,Inf),
                labels=c("<=2004","2005-09","2010-14","2015-19","2020-"))

set.seed(42)
perp <- max(5, floor(nrow(gamma)/50))
tsne_res <- Rtsne(gamma, perplexity=perp, theta=0.5,
                  check_duplicates=FALSE, pca=TRUE, max_iter=1000)

tsne_df <- tibble(
  x = tsne_res$Y[,1],
  y = tsne_res$Y[,2],
  year = doc_meta$year,
  year_bin = year_bin,
  pmid = doc_meta$pmid,
  top_topic = max.col(gamma, ties.method="first")
)

p_tsne <- ggplot(tsne_df, aes(x,y)) +
  geom_point(aes(shape=year_bin, alpha=0.9), size=2) +
  scale_alpha(guide="none") +
  labs(title="t-SNE of Document–Topic Mixtures",
       x="t-SNE 1", y="t-SNE 2", shape="Year bin") +
  theme_minimal(base_size=12)

f_tsne <- file.path(DIR_FIG, paste0("topicmap_", STAMP, ".pdf"))
ggsave(f_tsne, p_tsne, width=9, height=7)
message("Saved: ", f_tsne)

############################################################
# 6) 自動ラベル → ラベル付き t-SNE（Fig2B 拡張版）
############################################################

label_df <- top_terms %>%
  group_by(topic) %>%
  slice_max(beta, n=3, with_ties=FALSE) %>%
  summarise(label = paste(term, collapse=", "), .groups="drop") %>%
  mutate(topic = as.integer(topic))

centers <- tsne_df %>%
  group_by(top_topic) %>%
  summarise(cx=median(x), cy=median(y), .groups="drop") %>%
  rename(topic = top_topic) %>%
  left_join(label_df, by="topic")

p_tsne_lab <- ggplot(tsne_df, aes(x,y)) +
  geom_point(aes(color=factor(top_topic), shape=year_bin), size=2, alpha=0.8) +
  ggrepel::geom_label_repel(
    data=centers,
    aes(x=cx, y=cy, label=label),
    size=3.5, fontface="bold",
    fill="white", label.size=0.2,
    box.padding=0.5
  ) +
  labs(title="t-SNE with Topic Labels") +
  theme_minimal(base_size=12)

f_tsne_lab <- file.path(DIR_FIG, paste0("topicmap_labeled_", STAMP, ".pdf"))
ggsave(f_tsne_lab, p_tsne_lab, width=10, height=8)
message("Saved: ", f_tsne_lab)

############################################################
# 7) 年次トレンド（Fig2の補助図）
############################################################

K <- ncol(gamma)
topic_year <- as.data.frame(gamma) %>%
  mutate(year = doc_meta$year) %>%
  filter(!is.na(year)) %>%
  pivot_longer(starts_with("V"), names_to="topic", values_to="gamma") %>%
  mutate(topic = as.integer(gsub("\\D","",topic))) %>%
  group_by(year, topic) %>%
  summarise(mean_gamma = mean(gamma), n_docs=n(), .groups="drop")

# 1) line
p_trend <- ggplot(topic_year, aes(year, mean_gamma)) +
  geom_line() + geom_point(size=0.8) +
  facet_wrap(~ paste0("T", topic), ncol=4, scales="free_y") +
  labs(title="Yearly mean topic proportion (γ)") +
  theme_minimal(base_size=12)
f_trend <- file.path(DIR_FIG, paste0("topic_year_trends_", STAMP, ".pdf"))
ggsave(f_trend, p_trend, width=12, height=7)

# 2) heatmap
p_heat <- ggplot(topic_year, aes(year, factor(topic), fill=mean_gamma)) +
  geom_tile() +
  scale_y_discrete(labels=function(v) paste0("T",v)) +
  labs(title="Heatmap of mean γ", x="Year", y="Topic") +
  theme_minimal(base_size=12)
f_heat <- file.path(DIR_FIG, paste0("topic_year_heatmap_", STAMP, ".pdf"))
ggsave(f_heat, p_heat, width=12, height=6)

# 3) dominant topic count
dominant_topic <- max.col(gamma, ties.method="first")
dom_df <- tibble(year=doc_meta$year, topic=dominant_topic) %>%
  filter(!is.na(year)) %>%
  count(year,topic, name="n")

p_dom <- ggplot(dom_df, aes(year,n,color=factor(topic))) +
  geom_line() + geom_point(size=0.9) +
  labs(title="Counts of dominant topic per year") +
  theme_minimal(base_size=12)
f_dom <- file.path(DIR_FIG, paste0("topic_dominant_counts_", STAMP, ".pdf"))
ggsave(f_dom, p_dom, width=12, height=6)

############################################################
# END
############################################################
message("=== DONE (02) ===")
