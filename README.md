# Trigeminal Chemosensitivity Analysis

Computational pipeline for characterizing trigeminal sensitivity phenotypes and identifying their modulatory factors. This repository accompanies a peer-reviewed scientific manuscript describing a two-stage analysis: (1) characterization of the trigeminal data space and derivation of a composite trigeminal phenotype via unsupervised clustering, and (2) identification of modulators of both the phenotype and three psychophysical measures (AmmoLa intensity, CO₂ threshold, Lateralization).

**Development environment:** R 4.5.3 and 4.6 for Linux

---

## Analysis Design

**Two-stage workflow:**
- **Stage 1:** Exploratory analysis (correlations, PCA) → Unsupervised clustering for phenotypes → Supervised cluster classification
- **Stage 2:** Regression models for psychophysical measures → Modulator identification

**Validation framework:** Train/validation split (80%/20%) before imputation. All models trained on training set, validated on held-out data.

See manuscript Methods for full details.

## Table of Contents

1. [Environment](#environment)
2. [Repository Structure](#repository-structure)
3. [Reported Analyses (Manuscript)](#reported-analyses-manuscript)
4. [Exploratory & Supporting Analyses](#exploratory--supporting-analyses)
5. [Data Processing](#data-processing)
6. [Usage](#usage)

---

## Environment

**Language & Platform:**
- R (version 4.x, Linux environment recommended)
- Python (optional, for supplementary analyses)

**Key R Packages:**
- Data handling: `dplyr`, `tidyr`, `reshape2`
- Visualization: `ggplot2`, `ggpubr`, `cowplot`, `patchwork`, `ComplexHeatmap`
- Distribution modeling: `opGMMassessment`, `DataVisualizations`
- Imputation: `missForest`
- Dimensionality reduction: `FactoMineR`, `factoextra`
- Clustering: `cluster`, `NbClust`, `dendextend`
- Statistical analysis: `psych`, `MASS`, `boot`, `effsize`

## Repository Structure

```
📂 R/
├── Core Infrastructure
│   ├── utils.R                                     [General utility functions: 25 helpers]
│   ├── globals.R                                   [Global config: libraries, colors, labels]
│   └── ProjectionsBiomed_MainFunctions_6_1core.R  [Helper library for projection methods]
│
├── Data Processing Pipeline (Required)
│   ├── read_data_and_basic_corrections.R          [Data import & validation]
│   └── build_analysis_dataset.R                   [Dataset construction: 242 variables]
│
├── *** REPORTED ANALYSES (Manuscript) ***
│   ├── Data Preparation
│   │   ├── read_data_and_basic_corrections.R       [Data import & validation]
│   │   └── build_analysis_dataset.R                [Dataset construction: 241 variables]
│   │
│   ├── Stage 1: Trigeminal Data Space Characterization & Phenotyping
│   │   ├── 01_stage1_trigeminal_data_space_exploration.R
│   │   │   └─ Exploratory analysis: correlations, PCA, variable importance
│   │   │
│   │   ├── 02_stage1_trigeminal_clustering_phenotype_derivation.R
│   │   │   └─ Unsupervised clustering: composite trigeminal phenotype derivation
│   │   │
│   │   ├── 03_stage1_trigeminal_cluster_classification.R
│   │   │   └─ Supervised classification of cluster membership
│   │   │
│   │   └── 04_stage1_trigeminal_cluster_interpretation.R
│   │       └─ Clinical and biological cluster interpretation
│   │
│   └── Stage 2: Psychophysical Regression & Modulator Identification
│       ├── 05_psychophysical_regression_analysis.R
│       │   └─ Regression models: AmmoLa, CO₂ threshold, Lateralization
│       │
│       ├── 06_stage2_modulator_identification_psychophysical.R
│       └── 07_stage2_modulator_identification_clusters.R
│           └─ Identification of demographic, medical, lifestyle modulators
│
├── *** EXPLORATORY ANALYSES (Not part of core reported pipeline) ***
│   ├── trigeminal_measures_clustering_y_trig_2.R
│   │   └─ Alternative clustering approach (kept for reproducibility)
│   │
│   ├── descriptive_statistics_overview.R             [Cohort characteristics]
│   ├── smoking_overview.R                            [Smoking patterns]
│   ├── facial_pain_overview.R                        [Facial pain]
│   ├── chronic_diseases_overview.R                   [Disease comorbidity]
│   ├── covid_overview.R                              [COVID-19 effects]
│   └── nasal_breathing_and_ENT_surgery_overview.R    [Nasal/ENT factors]
│
├── Data Files
│   └─ [Not version-controlled: raw and processed datasets stored locally]
│
└── TrigeminalSensitivity.Rproj

📂 Python/
├── trigeminal_measures_PCA.py                    [Principal Component Analysis of trigeminal variables]
├── compute_PCA_altColors.py                      [PCA visualization with alternative color scheme]
├── cABC_analysis.py                              [ABC analysis for variable importance ranking]
├── trigeminal_measures_exploreTukey.py           [Tukey's ladder of powers distribution exploration]
└── [explore_tukey_lop.py](https://github.com/JornLotsch/PythonDataExploration/blob/main/explore_tukey_lop.py) [External: distribution exploration utility]
```

## Reported Analyses (Manuscript)

### Data Preparation (Required)

**Files:** [read_data_and_basic_corrections.R](R/read_data_and_basic_corrections.R), [build_analysis_dataset.R](R/build_analysis_dataset.R)

- Raw data import and validation
- Variable recoding: 8 categories → 241 variables
- Train/validation split (80%/20%) before imputation
- Imputation: missRanger (random forest, k=5)
- **Output: 21 trigeminal variables** (18 TriFunQ items + 3 psychophysical measures)

---

### Stage 1: Characterization of Trigeminal Data Space

#### Trigeminal Data Space Exploration

**File:** [01_stage1_trigeminal_data_space_exploration.R](R/01_stage1_trigeminal_data_space_exploration.R)

Spearman correlations, PCA, and variable importance analysis of the 21 trigeminal variables.

#### Derivation of Composite Trigeminal Phenotype (Trigeminal Clusters)

**File:** [02_stage1_trigeminal_clustering_phenotype_derivation.R](R/02_stage1_trigeminal_clustering_phenotype_derivation.R)

Unsupervised clustering (48 combinations of 6 projections × 8 algorithms) to identify trigeminal phenotypes on training data (n=800). Validation and cross-dataset concordance assessment.

#### Supervised Classification of Trigeminal Cluster Membership

**File:** [03_stage1_trigeminal_cluster_classification.R](R/03_stage1_trigeminal_cluster_classification.R)

Supervised classification of cluster membership using multinomial logistic regression, penalized variants, and random forest. Evaluated on validation set.

#### Cluster Interpretation

**File:** [04_stage1_trigeminal_cluster_interpretation.R](R/04_stage1_trigeminal_cluster_interpretation.R)

Clinical and biological interpretation of derived clusters via TriFunQ psychometric dimensions.

---

### Stage 2: Regression Analysis of Psychophysical Measures & Modulator Identification

#### Regression Models: Predicting Psychophysical Measures

**File:** [05_psychophysical_regression_analysis.R](R/05_psychophysical_regression_analysis.R)

Regression models predicting each psychophysical measure (AmmoLa, CO₂, Lateralization) from remaining 20 trigeminal variables. OLS, penalized, and random forest approaches with validation on held-out data.

#### Identification of Modulators

**Files:** [06_stage2_modulator_identification_psychophysical.R](R/06_stage2_modulator_identification_psychophysical.R), [07_stage2_modulator_identification_clusters.R](R/07_stage2_modulator_identification_clusters.R)

Identification of demographic, medical, and lifestyle factors associated with trigeminal sensitivity and cluster membership using same regression pipeline (220 non-trigeminal candidate predictors).

---

## Exploratory & Supporting Analyses

### Exploratory Descriptive Context

**Files:** [descriptive_statistics_overview.R](R/descriptive_statistics_overview.R), [smoking_overview.R](R/smoking_overview.R), [facial_pain_overview.R](R/facial_pain_overview.R), [chronic_diseases_overview.R](R/chronic_diseases_overview.R), [covid_overview.R](R/covid_overview.R), [nasal_breathing_and_ENT_surgery_overview.R](R/nasal_breathing_and_ENT_surgery_overview.R)

Population-level descriptive analyses. Not part of core reported pipeline.

### Alternative Clustering Approach

**File:** [trigeminal_measures_clustering_y_trig_2.R](R/trigeminal_measures_clustering_y_trig_2.R)

Alternative clustering solution; not included in manuscript but available for methodological completeness.

---

## Data Processing

See **Data Preparation** section above and manuscript Methods for preprocessing details.

**Variable Transformations:** [utils.R](R/utils.R)
- **AmmoLa intensity:** Reflected sign-preserving log transformation → rescaled to [0, 3]
- **CO₂ threshold:** Sign-preserving log transformation (sign-inverted) → rescaled to [0, 3]
- **Lateralization:** Linear rescaling to [0, 3]

---

## Usage

**Code Organization:**
All analysis scripts automatically load required utilities through the following dependency chain:
```
utils.R → globals.R → analysis scripts
```

**To Reproduce the Manuscript Analyses (in order):**

**Step 0: Data Preparation (required)**
```r
source("R/read_data_and_basic_corrections.R")
source("R/build_analysis_dataset.R")
```
Creates: Train/validation split (n=800/201), 21-variable trigeminal matrix

**Optional: Exploratory Descriptive Context** (population characterization, not part of core pipeline)
```r
source("R/descriptive_statistics_overview.R")
source("R/smoking_overview.R")
source("R/facial_pain_overview.R")
source("R/chronic_diseases_overview.R")
source("R/covid_overview.R")
source("R/nasal_breathing_and_ENT_surgery_overview.R")
```
Produces: Supporting contextual characterization of cohort

---

**Stage 1: Trigeminal Data Space Characterization** (core reported analyses)

**Step 1a: Trigeminal Data Space Exploration**
```r
source("R/01_stage1_trigeminal_data_space_exploration.R")
```
Produces: Correlations, PCA results, variable importance

**Step 1b: Derive Composite Trigeminal Phenotype (Clusters)**
```r
source("R/02_stage1_trigeminal_clustering_phenotype_derivation.R")
```
Produces: Trigeminal cluster assignments, quality metrics

**Step 1c: Interpret & Validate Clusters**
```r
source("R/03_stage1_trigeminal_cluster_classification.R")
source("R/04_stage1_trigeminal_cluster_interpretation.R")
```
Produces: Cluster characterization, classification models, interpretation

**Stage 2: Regression & Modulator Analyses**

**Step 2a: Regression Models of Psychophysical Measures**
```r
source("R/05_psychophysical_regression_analysis.R")
```
Produces: Regression models for AmmoLa intensity, CO₂ threshold, Lateralization

**Step 2b: Identify Modulators of Trigeminal Sensitivity**
```r
source("R/06_stage2_modulator_identification_psychophysical.R")
source("R/07_stage2_modulator_identification_clusters.R")
```
Produces: Modulator identification, feature importance

**Note:** All scripts automatically load utilities via `globals.R`. No need to manually source `utils.R` or `ProjectionsBiomed_MainFunctions_6_1core.R`.

---

## Citation

- Lötsch J, Weise S, Himmelspach A, Bormann A, and Hummel T
Phenotypes of trigeminal chemosensitivity in a general-population sample 2026 (in preparation)


## License

[Code and documentation](https://github.com/JornLotsch/TrigeminalChemosensitivity)

CC-BY 4.0

*This repository contains analysis code for scientific research. Use restricted to authorized collaborators and publication in peer-reviewed venues.*
