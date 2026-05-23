# Trigeminal Chemosensitivity Analysis

Computational pipeline accompanying a scientific manuscript on trigeminal chemosensitivity phenotypes and their candidate modulators.

This repository is intended primarily for reproducibility and review. It contains the scripts used for data preparation, exploratory analyses, phenotype derivation, and predictive modeling. Detailed methodological rationale and interpretation are described in the manuscript and are not repeated here.

---

## Overview

The analyses are organized into three stages:

1. **Preliminary characterization of observed data**
   Characterization of psychophysical trigeminal measures in the non-imputed dataset (`n = 1,001`).

2. **Exploratory analysis**
   Characterization of the trigeminal data space and screening of candidate modulators in the complete imputed dataset (`n = 1,001`).

3. **Machine learning framework**
   Derivation of a composite trigeminal phenotype and evaluation of predictive models using a training/validation split (`n = 800/201`).

A key design feature is that the **training/validation split was performed before imputation**, and imputation was carried out independently in both subsets to prevent information leakage.

---

## Repository Structure
```

R/
├── Core infrastructure
│   ├── utils.R
│   ├── globals.R
│   └── ProjectionsBiomed_MainFunctions_6_1core.R
│
├── Data preparation
│   ├── read_data_and_basic_corrections.R
│   └── build_analysis_dataset.R
│
├── Reported analyses
│   ├── 01_stage1_2_trigeminal_data_space_exploration.R
│   ├── 02_stage3_trigeminal_clustering_phenotype_derivation.R
│   ├── 03_stage3_trigeminal_cluster_interpretation.R
│   ├── 04_stage3_psychophysical_regression_analysis.R
│   ├── 05_stage3_modulator_identification_psychophysical.R
│   ├── 06_stage3_modulator_identification_clusters.R
│   └── 07_stage3_trigeminal_cluster_modulators_interpretation.R
│
└── Supporting / exploratory scripts
    ├── descriptive_statistics_overview.R
    ├── smoking_overview.R
    ├── facial_pain_overview.R
    ├── chronic_diseases_overview.R
    ├── covid_overview.R
    └── nasal_breathing_and_ENT_surgery_overview.R

Python/
├── trigeminal_measures_PCA.py
├── compute_PCA_altColors.py
├── cABC_analysis.py
├── trigeminal_measures_exploreTukey.py
└── explore_tukey_lop.py
```

Raw and processed data files are stored locally and are **not version-controlled** in this repository.

---

## Analysis Summary

### Stage 1: Preliminary characterization
**File:** `R/read_data_and_basic_corrections.R`

Performed on observed, non-imputed data:
- Fisher’s exact tests for overlap of high-sensitivity classifications across psychophysical measures
- Spearman correlations among psychophysical measures
- Sex differences assessed using Kruskal–Wallis tests with effect sizes

### Stage 2: Exploratory analysis
**File:** `R/01_stage1_2_trigeminal_data_space_exploration.R`

Performed on the complete imputed dataset obtained by merging the independently imputed training and validation subsets:
- Correlation analysis of 21 trigeminal variables comprising 3 psychophysical measures and 18 TriFunQ items
- Principal component analysis
- Variable-importance ranking
- Exploratory screening of 220 candidate modulators against psychophysical measures

These analyses are descriptive and hypothesis-generating.

### Stage 3: Machine learning framework
**Files:**
- `R/02_stage3_trigeminal_clustering_phenotype_derivation.R`
- `R/03_stage3_trigeminal_cluster_interpretation.R`
- `R/04_stage3_psychophysical_regression_analysis.R`
- `R/05_stage3_modulator_identification_psychophysical.R`
- `R/06_stage3_modulator_identification_clusters.R`
- `R/07_stage3_trigeminal_cluster_modulators_interpretation.R`

Performed using training/validation separation:
- Unsupervised derivation of a composite trigeminal phenotype (clustering)
- Factor analysis for phenotype interpretation and hypothesis testing
- Supervised prediction of psychophysical measures
- Supervised classification of cluster membership
- Identification of non-trigeminal candidate modulators (psychophysical and cluster-based)

All model development, including clustering, feature selection, and hyperparameter tuning, is confined to the training set; generalization is evaluated in the held-out validation set.

---

## Environment

Tested under:
- **R 4.5.x / 4.6** on Linux
- **Python 3.12.x** on Linux

Main R packages used include:
`dplyr`, `tidyr`, `ggplot2`, `missRanger`, `FactoMineR`, `factoextra`, `cluster`, `NbClust`, `psych`, `glmnet`, `randomForest`, `Boruta`, `nnet`, `boot`, `car`, `ComplexHeatmap`, `uwot`, `Rtsne`, `fastICA`.

Additional packages are cited in the scripts and manuscript.

---

## Data and Processing Notes

- Analyses are based on `n = 1,001` participants.
- The recoded analytic dataset comprised **241 variables**.
- The reported analyses use **21 trigeminal variables** and **220 non-trigeminal candidate modulators**.
- The training/validation split (`80%/20%`) was created **before imputation**.
- Missing data were imputed independently in the training and validation sets using `missRanger`.
- Additional preprocessing and variable transformations are implemented in the data-preparation scripts and helper functions in `utils.R`.

For detailed methodological justification, see the manuscript.

---

## Usage

All analysis scripts load required helper functions through the internal dependency chain:

`utils.R -> globals.R -> analysis scripts`

### 1. Data preparation

`source("R/read_data_and_basic_corrections.R")`
`source("R/build_analysis_dataset.R")`

### 2. Preliminary characterization in observed data

`source("R/read_data_and_basic_corrections.R")`

### 3. Exploratory analysis in the complete imputed dataset

`source("R/01_stage1_trigeminal_data_space_exploration.R")`

### 4. Machine learning analyses

`source("R/01_stage1_2_trigeminal_data_space_exploration.R")`
`source("R/03_stage3_trigeminal_cluster_classification.R")`
`source("R/04_stage3_trigeminal_cluster_interpretation.R")`
`source("R/05_stage3_psychophysical_regression_analysis.R")`
`source("R/06_stage3_modulator_identification_psychophysical.R")`
`source("R/07_stage3_modulator_identification_clusters.R")`

### 5. Optional supporting analyses

`source("R/descriptive_statistics_overview.R")`
`source("R/smoking_overview.R")`
`source("R/facial_pain_overview.R")`
`source("R/chronic_diseases_overview.R")`
`source("R/covid_overview.R")`
`source("R/nasal_breathing_and_ENT_surgery_overview.R")`

---

## Citation

Lötsch J, Weise S, Himmelspach A, Bormann A, Hummel T.
**Phenotypes of trigeminal chemosensitivity in a general-population sample**.
2026, in preparation.

---

## License

[Code and documentation](https://github.com/JornLotsch/TrigeminalChemosensitivity)

CC BY 4.0