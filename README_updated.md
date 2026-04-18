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
  14b_Fig4_Integrated_Outcome_Synthesis_from_provenance.R
  14c_Build_Supplementary_Table_S7_from_Fig4_Provenance.R
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

The released analyses are tied to this dictionary file name. Minor internal refinements of regex content were incorporated into the final released file without changing the dictionary tag or file name; for reproducibility, use the released file in `dic/` unchanged.

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

- **14b_Fig4_Integrated_Outcome_Synthesis_from_provenance.R**  
  Manuscript-facing Figure 4 layout rebuilt from the 14a provenance outputs. This script derives shared versus outcome-centred groupings and writes the final manuscript-style Figure 4 panel.

- **14c_Build_Supplementary_Table_S7_from_Fig4_Provenance.R**  
  Builds Supplementary Table S7 from the Figure 4 provenance outputs by applying the manuscript-facing sharedness logic to the 14a exports.

- **15_Supp_DesignTier_Weighted_Sensitivity.R**  
  Generates design-tier weighted sensitivity outputs and Supplementary Tables S2–S5.

- **16_SuppFigS4_CaseReport_ABC_Only.R**  
  Generates the exploratory case-report AE-ILD ABC supplementary figure from case-report-specific outputs.

- **17_TimeSlice_BuildAndAnalyze.R**  
  Self-contained temporal holdout build script. It reruns grouped discovery/full/holdout analyses from a frozen original-article corpus and fixed dictionary placed under `temporal_holdout/input/`.

- **18_TimeSlice_CompareAndSummarize.R**  
  Compares discovery versus holdout outputs from script 17 and writes the temporal holdout summary and revised Supplementary Figure S5.

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
16. `14b_Fig4_Integrated_Outcome_Synthesis_from_provenance.R`
17. `14c_Build_Supplementary_Table_S7_from_Fig4_Provenance.R`
18. `15_Supp_DesignTier_Weighted_Sensitivity.R`
19. `16_SuppFigS4_CaseReport_ABC_Only.R`

Temporal holdout scripts are run separately:

1. `17_TimeSlice_BuildAndAnalyze.R`
2. `18_TimeSlice_CompareAndSummarize.R`

## Temporal holdout notes

The temporal holdout scripts are intentionally self-contained and separate from the main manuscript pipeline. They expect the following fixed input layout:

```text
temporal_holdout/
  input/
    articles_main_original_frozen.csv
    ra_ild_dictionary_analysis_v1_genetics.csv
  output/
  R/
    17_TimeSlice_BuildAndAnalyze.R
    18_TimeSlice_CompareAndSummarize.R
```

Script 17 writes the discovery/full/holdout slice outputs and the manifest. Script 18 reads that manifest and writes Supplementary Table S6 and the revised Supplementary Figure S5 outputs.

## Reproducibility notes

- The released corpus window is fixed by `DATE_MAX = "2025/12/31"` in `00_setup.R`.
- The released analysis dictionary is fixed as `ra_ild_dictionary_analysis_v1_genetics.csv`.
- Figure 4 reproducibility is anchored to **14a**, which writes the term-level provenance tables.  
  **14b** is the manuscript-facing regrouping/layout script, and **14c** generates Supplementary Table S7 from those provenance outputs.
- Dataset S1 and Dataset S2 are publication-facing names used for public release.  
  Internal pipeline outputs may have longer tag-stamped file names, as noted above.
- PubMed content changes over time. For exact manuscript reproduction, use the released `Dataset_S1_pmids1.csv`.
- Design-tier weighted analyses (Supplementary Tables S2–S5) and temporal holdout outputs (Supplementary Table S6) are sensitivity / robustness analyses and do not replace the main manuscript analyses.

## Data redistribution note

This repository releases PMID-level identifiers and derived analysis tables. PubMed abstract text itself should not be redistributed here unless licensing and source terms clearly allow it.

## Contact

**Shinji Maeda, MD, PhD**  
Department of Respiratory Medicine, Allergy and Clinical Immunology,  
Nagoya City University Graduate School of Medical Sciences, Nagoya, Japan  
Email: `snb51961@med.nagoya-cu.jp`
