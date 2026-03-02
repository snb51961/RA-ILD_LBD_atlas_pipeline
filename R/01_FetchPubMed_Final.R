# =========================================================
# 01_fetch_pubmed.R  (based on your hardened code)
# - Retrieval window fixed by DATE_MAX (<= 2025/12/31)
# - Year counts + changepoint red dashed lines
# =========================================================

source(file.path(getwd(), "00_setup_Final.R"))

quiet_install <- function(pkgs){
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

quiet_install(c("rentrez","xml2","tibble","dplyr","readr","ggplot2","changepoint","systemfonts","showtext"))

# Font selection for PDF rendering (choose an available system font)
cands  <- c("Hiragino Sans","Hiragino Kaku Gothic ProN","Yu Gothic","IPAexGothic","Noto Sans CJK JP","Arial","Helvetica","DejaVu Sans")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
base_font <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

STAMP_DATE <- format(Sys.Date(), "%Y%m%d")

# Upper bound year (derived from DATE_MAX)
YEAR_MAX <- suppressWarnings(as.integer(substr(DATE_MAX, 1, 4)))
if (!is.finite(YEAR_MAX)) stop("DATE_MAX is invalid: ", DATE_MAX)

# ------------------ 1) PubMed search -----------------------
QUERY <- paste0(
  '(',
    '"Arthritis, Rheumatoid"[MH] OR "rheumatoid arthritis"[Title/Abstract]',
  ') AND (',
    '"interstitial lung disease"[Title/Abstract] OR ILD[Title/Abstract] OR ',
    '"interstitial pneumonia"[Title/Abstract] OR "pulmonary fibrosis"[Title/Abstract] OR ',
    '"diffuse parenchymal lung disease"[Title/Abstract]',
  ') AND hasabstract[text] NOT (editorial[Publication Type] OR letter[Publication Type])',
  ' AND ("', YEAR_MIN, '/01/01"[pdat] : "', DATE_MAX, '"[pdat])'
)

log_msg("Searching PubMed ...")
s0 <- rentrez::entrez_search(db="pubmed", term=QUERY, use_history=TRUE, retmax=0)
TOTAL <- as.integer(s0$count); log_msg("TOTAL:", TOTAL)

pmids <- character(0); step <- 5000L
if (TOTAL > 0) {
  for (start in seq(0L, TOTAL-1L, by=step)) {
    ids <- rentrez::entrez_search(db="pubmed", term=QUERY, retstart=start, retmax=step)$ids
    pmids <- c(pmids, ids)
    log_msg(sprintf("PMID: %d / %d", length(pmids), TOTAL))
    Sys.sleep(0.34)
  }
}
pmids <- unique(pmids)

# Save PMIDs
readr::write_csv(tibble::tibble(pmid=pmids),
                 file.path(DIR_RAW, sprintf("pmids_%s.csv", CORPUS_TAG)))

# ---- EFetch (XML) ----
extract_one <- function(n){
  z <- function(xpath) xml2::xml_text(xml2::xml_find_first(n, xpath))
  zx <- function(xpath) paste(xml2::xml_text(xml2::xml_find_all(n, xpath)), collapse = " ")
  zxs<- function(xpath) paste(xml2::xml_text(xml2::xml_find_all(n, xpath)), collapse = ";")
  ti  <- z(".//ArticleTitle")
  ab  <- zx(".//Abstract/AbstractText")
  yr  <- z(".//PubDate/Year"); if (is.na(yr) || yr=="") yr <- z(".//ArticleDate/Year")
  pt  <- zxs(".//PublicationType")
  mesh<- zxs(".//MeshHeading/DescriptorName")
  jr  <- z(".//Journal/Title")
  aff <- z(".//AuthorList/Author[last()]/AffiliationInfo/Affiliation")
  tibble::tibble(
    pmid = z(".//PMID"), title=ti, abstract=ab,
    year=suppressWarnings(as.integer(yr)), pubtype=pt, mesh=mesh, journal=jr,
    affiliation_last=aff
  )
}

BATCH <- 200L
rows <- list()
if (length(pmids) > 0) {
  for (i in seq(1, length(pmids), by=BATCH)) {
    j <- min(i+BATCH-1, length(pmids))
    xmltxt <- rentrez::entrez_fetch(db="pubmed", id=pmids[i:j], rettype="xml", parsed=FALSE)
    doc    <- xml2::read_xml(xmltxt)
    nodes  <- xml2::xml_find_all(doc, ".//PubmedArticle")
    rows[[length(rows)+1]] <- dplyr::bind_rows(lapply(nodes, extract_one))
    log_msg(sprintf("EFetch: %d / %d", j, length(pmids)))
    Sys.sleep(0.34)
  }
}

articles_raw <- dplyr::bind_rows(rows)

# Guardrail: exclude records beyond YEAR_MAX (double-check)
max_year_raw <- suppressWarnings(max(articles_raw$year, na.rm = TRUE))
if (is.finite(max_year_raw) && max_year_raw > YEAR_MAX) {
  log_msg("WARNING: raw records include year >", YEAR_MAX, " (max=", max_year_raw, "). They will be excluded.")
}

articles <- articles_raw |>
  dplyr::filter(!is.na(year), year >= YEAR_MIN, year <= YEAR_MAX) |>
  dplyr::mutate(text = paste(title, abstract, sep = " "))

if (nrow(articles) && any(articles$year > YEAR_MAX, na.rm=TRUE)) {
  stop("ERROR: year >", YEAR_MAX, " remains after filtering.")
}

# Save (backward-compatible file names)
saveRDS(articles, file.path(DIR_RAW,  sprintf("articles_%s.rds", STAMP_DATE)))
readr::write_csv(articles, file.path(DIR_PROC, sprintf("articles_%s.csv", STAMP_DATE)))
log_msg("Saved articles n=", nrow(articles))

# Annual publication counts + changepoint detection
year_counts <- articles |> dplyr::count(year, name="n") |> dplyr::arrange(year)
readr::write_csv(year_counts, file.path(DIR_PROC, sprintf("year_counts_%s.csv", STAMP_DATE)))

p_year <- ggplot2::ggplot(year_counts, ggplot2::aes(year, n)) +
  ggplot2::geom_col() +
  ggplot2::labs(title="Annual publication counts (RA-ILD corpus)", x="Year", y="Count") +
  ggplot2::theme_minimal(base_size=12) +
  ggplot2::theme(text=ggplot2::element_text(family=base_font))

ggplot2::ggsave(file.path(DIR_FIG, "year_counts.pdf"), p_year, device=cairo_pdf_device, width=8, height=5)

cpt <- tryCatch(changepoint::cpt.meanvar(year_counts$n, method="PELT", penalty="BIC"), error=function(e) NULL)
if (!is.null(cpt)) {
  idx <- changepoint::cpts(cpt); change_years <- year_counts$year[idx]
  p_year_cp <- p_year + ggplot2::geom_vline(xintercept=change_years, linetype="dashed", color="red")
  ggplot2::ggsave(file.path(DIR_FIG, "year_counts_changepoint.pdf"), p_year_cp, device=cairo_pdf_device, width=8, height=5)
  log_msg("Changepoints:", paste(change_years, collapse=", "))
}

log_msg("=== DONE 01_fetch_pubmed ===")
