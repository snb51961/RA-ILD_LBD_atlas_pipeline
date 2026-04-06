# RA-ILD literature-derived atlas pipeline

This repository contains the R scripts, curated dictionary, and released datasets used for the RA-ILD literature-derived atlas manuscript.

## Repository structure

```text
R/
  00_setup.R
  01_FetchPubMed.R
  02_BuildHitsMatrixl.R
  03_CoocAndCollocation.R
  04_ABC_Rankings.R
  05_AC_NPMI.R
  06_SignedEffectsl.R
  07_AE_Ratio_TrendBreakl.R
  08_Supp_Nonreview_Sensitivity.R
  09_SuppFigS1_Topics_TemporalEvolution.R
  10_SuppFigS2_AC_Coherence.R
  11_SuppFigS3_ABC_Exhaustive.R
  12_Fig2_AEILD_Directional_and_Bridge.R
  13_Fig3_Biomarker_Atlas.R
  13_Fig3_Biomarker_Atlas_renumbered.R
  14a_Fig4_Integrated_Outcome_Synthesis_fully_datadriven.R
  14b_Fig4_Integrated_Outcome_Synthesis.R
  15_Supp_DesignTier_Weighted_Sensitivity.R
  16_SuppFigS4_CaseReport_ABC_Only.R
  17_TimeSlice_BuildAndAnalyze.R
  18_TimeSlice_CompareAndSummarize.R

dic/
  ra_ild_dictionary_analysis_v1_genetics.csv

Dataset/
  Dataset_S1_pmids1.csv
  Dataset_S2_signed_effects_summary_withCI.csv
```

## Included released datasets

### Dataset S1
**File name:** `Dataset_S1_pmids1.csv`  
PMID-level frozen manuscript corpus file corresponding to the fixed PubMed retrieval window and post-retrieval original-article corpus definition.

This publication-facing dataset corresponds to the internal pipeline output:
- `pmids_pm_1980_20251231.csv`

### Dataset S2
**File name:** `Dataset_S2_signed_effects_summary_withCI.csv`  
Signed-effect summary dataset underlying the outcome-centred directional summaries reported in the manuscript.

This publication-facing dataset corresponds to the internal pipeline output:
- `signed_effects_summary_withCI_pm_1980_20251231__analysis_v1_genetics.csv`

## Curated dictionary

### Dictionary file
**File name:** `dic/ra_ild_dictionary_analysis_v1_genetics.csv`

This is the fixed analysis dictionary used throughout the released pipeline. It contains:
- canonical terms
- regular expressions
- semantic classes
- ABC roles used in downstream interpretation

The released analyses are tied to this dictionary version. For reproducibility, use this file unchanged.

## Script overview

### Core pipeline
- **00_setup.R**  
  Defines project root, corpus tag, dictionary tag, output directories, and shared logging helpers.

- **01_FetchPubMed.R**  
  Performs broad PubMed retrieval, applies post-retrieval publication-type classification, and writes the main original-article corpus plus auxiliary corpora.

- **02_BuildHitsMatrixl.R**  
  Builds the dictionary-based `hit__` matrix from the original-article corpus, including `_var` parent aggregation and consistency fixes.

- **03_CoocAndCollocation.R**  
  Computes class-vs-outcome co-occurrence tables and optional bigram/trigram collocations.

- **04_ABC_Rankings.R**  
  Computes AE-ILD-centred ABC triad rankings, evidence PMIDs, and network edges.

- **05_AC_NPMI.R**  
  Computes A↔C co-occurrence probabilities, lift, and NPMI for downstream figures.

- **06_SignedEffectsl.R**  
  Performs rule-based signed-effect extraction from title/abstract contexts and writes sentence-level, article-level, summary, numeric, and biomarker outputs.

- **07_AE_Ratio_TrendBreakl.R**  
  Detects changepoints in yearly AE-ILD ratios for top A terms.

- **08_Supp_Nonreview_Sensitivity.R**  
  Computes non-review / non-case-report drug–outcome co-occurrence summaries.

### Figure and supplementary scripts
- **09_SuppFigS1_Topics_TemporalEvolution.R**  
  Generates topic-model outputs and the temporal evolution figure.

- **10_SuppFigS2_AC_Coherence.R**  
  Generates supplementary A↔C coherence heatmaps and scatterplots.

- **11_SuppFigS3_ABC_Exhaustive.R**  
  Generates the full AE-ILD-centred ABC supplementary figure.

- **12_Fig2_AEILD_Directional_and_Bridge.R**  
  Generates the main Figure 2 panels combining signed effects and AE-ILD bridge summaries.

- **13_Fig3_Biomarker_Atlas.R**  
  Main biomarker atlas script and atlas export.

- **13_Fig3_Biomarker_Atlas_renumbered.R**  
  Renumbered figure script version used for manuscript figure ordering.

- **14a_Fig4_Integrated_Outcome_Synthesis_fully_datadriven.R**  
  Primary Figure 4 synthesis script. Selects hotspot and biomarker/exploratory terms from analysis outputs and writes provenance tables.

- **14b_Fig4_Integrated_Outcome_Synthesis.R**  
  Refined manuscript-facing Figure 4 layout script. This is a presentation-oriented version derived from the analytical outputs and is not the main provenance-generating script.

- **15_Supp_DesignTier_Weighted_Sensitivity.R**  
  Generates design-tier weighted sensitivity outputs and Supplementary Tables S2–S5.

- **16_SuppFigS4_CaseReport_ABC_Only.R**  
  Generates the exploratory case-report AE-ILD ABC supplementary figure from case-report-specific outputs.

- **17_TimeSlice_BuildAndAnalyze.R**  
  Builds discovery/holdout/full temporal slices and computes grouped slice outputs.

- **18_TimeSlice_CompareAndSummarize.R**  
  Compares discovery vs holdout outputs and writes the temporal holdout summary and revised Supplementary Figure S5.

## Recommended run order

For the main released pipeline, the practical execution order is:

1. `00_setup.R`
2. `01_FetchPubMed.R`
3. `02_BuildHitsMatrixl.R`
4. `03_CoocAndCollocation.R`
5. `04_ABC_Rankings.R`
6. `05_AC_NPMI.R`
7. `06_SignedEffectsl.R`
8. `07_AE_Ratio_TrendBreakl.R`
9. `08_Supp_Nonreview_Sensitivity.R`
10. `09_SuppFigS1_Topics_TemporalEvolution.R`
11. `10_SuppFigS2_AC_Coherence.R`
12. `11_SuppFigS3_ABC_Exhaustive.R`
13. `12_Fig2_AEILD_Directional_and_Bridge.R`
14. `13_Fig3_Biomarker_Atlas.R`
15. `14a_Fig4_Integrated_Outcome_Synthesis_fully_datadriven.R`
16. `14b_Fig4_Integrated_Outcome_Synthesis.R`
17. `15_Supp_DesignTier_Weighted_Sensitivity.R`
18. `16_SuppFigS4_CaseReport_ABC_Only.R`

Temporal holdout scripts are run separately:
- `17_TimeSlice_BuildAndAnalyze.R`
- `18_TimeSlice_CompareAndSummarize.R`

## Reproducibility notes

- The released corpus window is fixed by `DATE_MAX = "2025/12/31"` in `00_setup.R`.
- The released analysis dictionary is fixed as `ra_ild_dictionary_analysis_v1_genetics.csv`.
- Figure 4 reproducibility is anchored to **14a**, which writes term-level provenance tables.  
  **14b** should be read as the refined manuscript-facing layout version.
- Dataset S1 and Dataset S2 are publication-facing names used for public release.  
  Internal pipeline outputs may have longer tag-stamped file names, as noted above.
- PubMed content changes over time. For exact manuscript reproduction, use the released `Dataset_S1_pmids1.csv`.

## Data redistribution note

This repository releases PMID-level identifiers and derived analysis tables. PubMed abstract text itself should not be redistributed here unless licensing and source terms clearly allow it.

## Contact

**Shinji Maeda, MD, PhD**  
Department of Respiratory Medicine, Allergy and Clinical Immunology,  
Nagoya City University Graduate School of Medical Sciences, Nagoya, Japan  
Email: `snb51961@med.nagoya-cu.jp`
