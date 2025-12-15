# RA-ILD LBD Atlas Pipeline
## Reproducible Text-Mining and Biomarker Evidence Atlas

This repository provides fully reproducible R code for the literature-based discovery (LBD)
analysis of rheumatoid arthritis–associated interstitial lung disease (RA-ILD),
as reported in the accompanying manuscript.


The code reproduces all main and supplementary figures, including:

- Topic modelling and t-SNE maps (Figure 2)  
- A–C coherence analysis using NPMI and lift (Figure 3)  
- ABC triad network analysis (Figure 4)  
- Signed-effect analysis for AE-ILD (Figure 5A–B)  
- White-map decomposition and candidate prioritization (Figure 5C)  
- Biomarker Evidence Atlas (Figure 6)

The repository is structured to separate **full reproducible analysis code** from **minimal Supplementary Codes** provided for transparency.

---

## 📁 Directory Structure

RA-ILD_LBD_atlas_pipeline/
├── core/ # Supplementary Codes (S1)
│ ├── S1_AC_coherence.R
│ ├── S1_white_map_core.R
│ └── S1_white_top10.R
│
├── dic/ # Dictionaries and curated knowledge tables
│ ├── raalid_v07.csv
│ ├── raalid_pub.csv
│ └── external_AC_evidence.csv
│
├── R/ # Full analysis pipeline (01–07)
│ ├── 01_corpus_build.R
│ ├── 02_topic_LDA_tSNE.R
│ ├── 03_coherence_AC_npmi_lift.R
│ ├── 04_ABC_triad_score.R
│ ├── 05_signed_effects_AE_ILD.R
│ ├── 06_whitemap_and_Bbreakdown.R
│ └── 07_biomarker_atlas.R
│
└── README.md


---

## ▶️ How to Reproduce All Manuscript Figures

### 1. Prepare the project directory

Clone this repository and set the project root as the working directory.  
All scripts use `here::here()` and assume the following structure:

data_raw/
data_proc/
output/
fig_pub/
fig_atlas/


### 2. Place input files

- PubMed-derived article files (`articles_*.csv`) → `data_proc/`
- Dictionary and curated knowledge tables → `dic/`

### 3. Run the full pipeline scripts in order

R/01_corpus_build.R
R/02_topic_LDA_tSNE.R
R/03_coherence_AC_npmi_lift.R
R/04_ABC_triad_score.R
R/05_signed_effects_AE_ILD.R
R/06_whitemap_and_Bbreakdown.R
R/07_biomarker_atlas.R


Figures will be generated automatically in `fig_pub/` and `fig_atlas/`.

---

## 🧩 Supplementary Codes (S1)

The `core/` directory contains **minimal, self-contained R scripts** corresponding to the Supplementary Codes in the manuscript:

- **S1_AC_coherence.R**  
  Core definition and visualization of A–C coherence (NPMI / lift).

- **S1_white_map_core.R**  
  Core construction of the white map defining A–C unknown candidates.

- **S1_white_top10.R**  
  Extraction and visualization of the top-ranked white map candidates.

These scripts are intended for conceptual transparency and do not require execution of the full pipeline.

---

## 📊 Dictionary Lifecycle and Supplementary Dataset

The initial dictionary (v07) used in this study was intentionally broad and manually curated to define the vocabulary space for RA-ILD–related literature mining.

The complete lifecycle of all initial dictionary terms—including:
- initial inclusion in the dictionary,
- occurrence in the literature corpus, and
- retention in the curated public dictionary

is provided as **Supplementary Dataset S1**.

**Supplementary Dataset S1 is not included in this GitHub repository** and is distributed separately via Zenodo together with the code release (see manuscript for DOI).

---

## 📌 Notes

- All scripts are designed for cross-platform compatibility.
- Font settings default to `family = "sans"`.
- Files in `dic/` are curated reference tables and are not modified by the pipeline.

---

## 📚 Citation

Please cite the associated manuscript when using this code.

A DOI for this repository and Supplementary Dataset S1 is issued via Zenodo.

---

## ❓ Contact

For questions regarding the analysis pipeline, dictionary design, or reproducibility,  
please contact the authors.
