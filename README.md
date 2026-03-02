# RA-ILD Literature-Based Discovery (ABC + Signed Effects) — Public R Code

This repository provides a reproducible R pipeline for literature-based discovery (LBD) in rheumatoid arthritis–associated interstitial lung disease (RA-ILD).

The pipeline retrieves PubMed records, builds document-level term-hit matrices using a curated dictionary, derives ABC triads and association statistics, and generates publication figures (Fig 2–6), including a rule-based signed-effect summary and a biomarker evidence atlas.

---

## Repository Contents

- **R/** — scripts (00–13) to run the pipeline and generate figures  
- **dic/** — curated dictionary  
  - `ra_ild_dictionary_analysis_v1_genetics.csv`  
- **out/** — outputs (figures/tables/logs; created automatically)  
- **data_raw/**, **data_proc/** — intermediate files (created automatically; typically git-ignored)

---

## Requirements

- **R (>= 4.2 recommended)**  

Required R packages are installed automatically by scripts where possible (or see `renv.lock` if provided).

---

## How to Run

### Option A (recommended): Run from the project root

1. Download/clone this repository to any location.
2. Set your working directory to the project root (the folder that contains `R/` and `dic/`).
3. Run scripts in order.

---

### Option B: Run from anywhere using `RAILD_ROOT`

Set an environment variable pointing to the project root.

In R:

```r
Sys.setenv(RAILD_ROOT = "/path/to/this/repo")
```

---

## Typical Run Order

1. `00_setup_Final.R` — global configuration (directories/tags; corpus window)  
2. `01_FetchPubMed_Final.R` — retrieve PubMed records and write `articles_*.csv` + PMID list  
3. `02_BuildHitsMatrix_Final.R` — build document × term hit matrix  
4. `03_CoocAndCollocation_Final.R` — co-mention / collocation preprocessing  
5. `04_ABC_Rankings_Final.R` — construct ABC triads and rankings  
6. `05_AC_NPMI_Final.R` — association statistics (NPMI, lift, OR, q-values)  
7. `06_Fig5_SignedEffects_Final.R` — rule-based signed-effect extraction (AE-ILD / progression / mortality)

---

### Optional Analyses

- `07_AE_Ratio_TrendBreak_Final.R` — trend-break analysis  
- `08_Sensitivity_Nonreview_Final.R` — sensitivity analysis excluding review / case-report  

---

## Figure Scripts (Final Versions Used in the Manuscript)

- **Fig 2:** `09_Fig2_Topics_Final_v2.R` (consistent topic colors + combined layout)  
- **Fig 3:** `10_Fig3_AC_Coherence_Final.R`  
- **Fig 4:** `11_Fig4_ABC_Final_v3_moreTriads.R`  
- **Fig 5:** `12_Fig5_SignedEffects_Summary_Final.R`  
- **Fig 6:** `13_Fig6_BiomarkerEvidenceAtlas_Final_FINAL.R`  

---

## Outputs

Figures and tables are written under `out/` (exact subfolders depend on tags).

Some scripts also write reproducibility logs (including `sessionInfo()` and key input file metadata) under:

```
out/.../log/
```

---

## Supplementary Files (Manuscript Submission + GitHub)

We provide the following supplementary items used for reproducibility:

### Supplementary Table S1  
`supplementary/Supplementary_Table_S1_PubMed_Retrieval_Metadata.docx`  
(PubMed query, fixed time window, and retrieval date)

### Supplementary Data S1 (PMID list; exact corpus definition)  
`supplementary/Supplementary_Data_S1_pmids.csv`  
(originally generated as `pmids_pm_1980_20251231.csv`)

### Supplementary Data S2 (Signed-Effects Summary Used for Fig 5–6)  
`supplementary/Supplementary_Data_S2_signed_effects_summary_withCI.csv`  
(originally generated as  
`signed_effects_summary_withCI_pm_1980_20251231__analysis_v1_genetics.csv`)

---

## Notes

This pipeline is intended for hypothesis generation (LBD).  
Ranking scores and signed effects are heuristic summaries and are not formal confirmatory inference.

PubMed contents can change over time; for exact reproducibility of the corpus, use the provided PMID list (Supplementary Data S1).

---

## Citation

If you use this code or the released corpus definition, please cite the accompanying manuscript and the repository release (e.g., Zenodo DOI if available).