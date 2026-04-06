# =========================================================
# 00_setup.R  (RA-ILD text-mining pipeline: analysis_v1_genetics)
# - ROOT auto-detect + override
# - corpus tag fixed by DATE_MAX (pdat upper bound)
# - shared dirs + logging
# =========================================================

suppressPackageStartupMessages({
  library(fs)
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
})

# -------------------------
# 0) ROOT (auto-detect + override)
# -------------------------
# Edit only this if you need to set a custom project root.
ROOT_OVERRIDE <- Sys.getenv("RAILD_ROOT", unset = NA_character_)

guess_root <- function() {
  # Prefer directories that already look like this project.
  candidates <- c(
    ROOT_OVERRIDE,
    getwd()
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) {
    if (dir_exists(p) && (dir_exists(path(p, "dic")) || dir_exists(path(p, "data_raw")))) return(p)
  }
  # Last resort
  getwd()
}

ROOT <- guess_root()

# -------------------------
# 1) Corpus window (pdat upper bound fixed)
# -------------------------
YEAR_MIN <- 1980L
DATE_MAX <- "2025/12/31"  # <= fixed upper bound for this release

# Tags encode the exact search window; include DATE_MAX to keep analyses comparable.
CORPUS_TAG <- sprintf("pm_%d_%s", YEAR_MIN, gsub("[^0-9]", "", DATE_MAX))  # e.g., pm_1980_20251231
DIC_TAG    <- "analysis_v1_genetics"

# -------------------------
# 2) Directories (legacy + new tagged output)
# -------------------------
DIR_RAW   <- path(ROOT, "data_raw")
DIR_PROC  <- path(ROOT, "data_proc")
DIR_DIC   <- path(ROOT, "dic")
DIR_FIG   <- path(ROOT, "fig")
DIR_OUT   <- path(ROOT, "output", CORPUS_TAG, DIC_TAG)
DIR_LOG   <- path(DIR_OUT, "log")
DIR_TABLE <- path(DIR_OUT, "table")
DIR_FIG2  <- path(DIR_OUT, "fig")

dir_create(c(DIR_RAW, DIR_PROC, DIR_DIC, DIR_FIG, DIR_OUT, DIR_LOG, DIR_TABLE, DIR_FIG2), recurse = TRUE)

# Dictionary file (one true source)
DIC_FILE <- path(DIR_DIC, "ra_ild_dictionary_analysis_v1_genetics.csv")

# -------------------------
# 3) Logging
# -------------------------
log_file <- path(DIR_LOG, "pipeline.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  message(msg)
  write_lines(msg, log_file, append = TRUE)
}

# -------------------------
# 4) Common helpers
# -------------------------
stamp_date <- function() format(Sys.Date(), "%Y%m%d")
stamp_time <- function() format(Sys.time(), "%Y%m%d_%H%M")

write_csv2 <- function(df, file) {
  readr::write_csv(df, file)
  log_msg("WROTE:", file)
  invisible(file)
}

# -------------------------
# 5) Session header
# -------------------------
set.seed(42)
log_msg("=== SETUP ===")
log_msg("ROOT:", ROOT)
log_msg("YEAR_MIN:", YEAR_MIN, "DATE_MAX:", DATE_MAX)
log_msg("CORPUS_TAG:", CORPUS_TAG, "DIC_TAG:", DIC_TAG)
log_msg("DIC_FILE:", DIC_FILE)