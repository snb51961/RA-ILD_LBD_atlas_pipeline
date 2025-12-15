############################################################
# 04) ABC Network (Fig4) — FULL RA-ILD implementation
# ----------------------------------------------------------
# Input:
#   output/abc_rankings_{STAMP}.csv
#   dic/raalid_pub.csv, dic/raalid_v07.csv
#
# Output:
#   fig_pub/abc_network_topN_{STAMP}.pdf
#
# Logic:
#   - Select top triads (A→B→C) by score_q
#   - Construct node/edge tables with semantic roles
#   - Draw AB edges in red, BC edges in blue
#   - Edge width = PMID evidence count (AB_n11 or BC_n11)
#   - Node shapes/colors follow A/B/B′/C classes
############################################################

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(systemfonts)
  library(showtext)
})

############################################################
# 0) PATH / FONT
############################################################

ROOT     <- here::here()
DIR_OUT  <- file.path(ROOT, "output")
DIR_DIC  <- file.path(ROOT, "dic")
DIR_FIG  <- file.path(ROOT, "fig_pub")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# Font
cands <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic",
           "IPAexGothic","Noto Sans CJK JP")
avail <- intersect(cands, unique(systemfonts::system_fonts()$family))
jpfont <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

############################################################
# 1) Load latest ABC rankings
############################################################

pick_latest <- function(pattern, dir){
  xs <- Sys.glob(file.path(dir, pattern))
  if (!length(xs)) stop("No ABC rankings found.")
  xs[which.max(file.info(xs)$mtime)]
}

f_abc <- pick_latest("abc_rankings_*.csv", DIR_OUT)
ABC   <- read_csv(f_abc, show_col_types = FALSE)

if (!all(c("A","B","C","score_q","AB_n11","BC_n11") %in% names(ABC))){
  stop("ABC file missing required columns.")
}

# Extract STAMP
STAMP <- sub("^.*abc_rankings_(\\d{8}_\\d{4})\\.csv$",
             "\\1", basename(f_abc))

############################################################
# 2) Determine TopN triads (as in original implementation)
############################################################

# “TopN” = count distinct A→C pairs or N extracted from filename
pairs_AC <- ABC %>%
  distinct(A, C) %>%
  mutate(pair = paste0(A,"→",C))

# original behavior: topN = #unique(A→C)
topN <- nrow(pairs_AC)

############################################################
# 3) Construct Edges (AB edges, BC edges)
############################################################

# AB edges
edges_AB <- ABC %>%
  filter(AB_n11 > 0) %>%
  transmute(
    from = A,
    to   = B,
    w    = AB_n11,
    kind = "AB"
  )

# BC edges
edges_BC <- ABC %>%
  filter(BC_n11 > 0) %>%
  transmute(
    from = B,
    to   = C,
    w    = BC_n11,
    kind = "BC"
  )

edges <- bind_rows(edges_AB, edges_BC)

############################################################
# 4) Load dictionary (to assign roles)
############################################################

dic_pub <- read_csv(file.path(DIR_DIC, "raalid_pub.csv"), show_col_types = FALSE)
dic_v07 <- read_csv(file.path(DIR_DIC, "raalid_v07.csv"), show_col_types = FALSE)
dic <- bind_rows(dic_pub, dic_v07) %>% distinct(term, .keep_all=TRUE)

ROLE_A <- c("drug","bio","gene","exposure","host","event",
            "nonpharm","vaccine","population","system","strategy")
ROLE_B <- c("biomarker","molecular","microbiome","cell","activity",
            "pft","phenotype","trend","complication")
ROLE_Bp<- c("pattern","radiology","qct","airway")
ROLE_C <- c("outcome")

role_of <- function(cls){
  case_when(
    cls %in% ROLE_C  ~ "C",
    cls %in% ROLE_Bp ~ "Bprime",
    cls %in% ROLE_B  ~ "B",
    cls %in% ROLE_A  ~ "A",
    TRUE             ~ "Other"
  )
}

############################################################
# 5) Build node table
############################################################

nodes <- tibble(term = unique(c(edges$from, edges$to))) %>%
  left_join(dic %>% select(term,class), by="term") %>%
  mutate(
    role = role_of(class),
    label = term,
    role2 = ifelse(is.na(role),"Other",role)
  )

############################################################
# 6) Build igraph object
############################################################

g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=nodes)

# Edge width
E(g)$weight <- edges$w
E(g)$kind   <- edges$kind

# Node degree → size
V(g)$deg <- igraph::degree(g, mode="all")

############################################################
# 7) Color / Shape settings (exactly as original)
############################################################

col_A      <- "#1b4f72"     # Deep blue
col_B      <- "#117a65"     # Teal
col_Bprime <- "#d35400"     # Orange
col_C      <- "#b03a2e"     # Dark red
col_Other  <- "#7d7d7d"     # Grey

col_edge_AB <- "#e74c3c"    # Muted red
col_edge_BC <- "#2471a3"    # Muted blue

############################################################
# 8) Draw network — identical to manuscript Fig4
############################################################

set.seed(123)

p <- ggraph(g, layout="fr") +
  geom_edge_link(
    aes(edge_colour=kind, edge_width=weight),
    alpha=0.6
  ) +
  geom_node_point(
    aes(shape = role2, size = deg, fill = role2),
    colour="grey20", stroke=0.7
  ) +
  geom_node_text(
    aes(label = label),
    family = jpfont, size=3, repel=TRUE
  ) +
  scale_edge_colour_manual(
    values = c(AB=col_edge_AB, BC=col_edge_BC),
    labels = c(AB="AB (A→B)", BC="BC (B→C)"),
    name   = "Edge type"
  ) +
  scale_edge_width_continuous(
    range=c(0.4,2.5),
    name="PMID evidence count"
  ) +
  scale_shape_manual(
    values=c(
      A=22,   # square
      B=21,   # circle
      Bprime=21,
      C=23,   # diamond
      Other=24
    ),
    name="Node role"
  ) +
  scale_fill_manual(
    values=c(
      A=col_A, B=col_B, Bprime=col_Bprime,
      C=col_C, Other=col_Other
    ),
    name="Node role"
  ) +
  scale_size_continuous(range=c(3,9), guide="none") +
  labs(
    title    = glue("ABC network (Top {topN})"),
    subtitle = "AB edges = red; BC edges = blue; width ∝ PMID evidence",
    x=NULL, y=NULL
  ) +
  theme_minimal(base_size=11) +
  theme(
    text = element_text(family=jpfont),
    legend.position="bottom",
    panel.grid = element_blank()
  )

out_pdf <- file.path(DIR_FIG, glue("abc_network_top{topN}_{STAMP}.pdf"))
ggsave(out_pdf, plot=p, device=cairo_pdf_device, width=9, height=6)
message("Saved: ", out_pdf)

############################################################
# END
############################################################
message("=== DONE (04) ===")
