# Trigeminal Chemosensitivity Analysis

Computational pipeline for analyzing trigeminal sensitivity measures and their phenotypic characterization through psychophysical testing, distribution modeling, and unsupervised clustering approaches.

## Table of Contents

1. [Environment](#environment)
2. [Repository Structure](#repository-structure)
3. [Core Analyses](#core-analyses)
   - [Distribution and Univariate Analysis](#distribution-and-univariate-analysis)
   - [Clustering and Phenotyping](#clustering-and-phenotyping)
4. [Statistical Methods](#statistical-methods)
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
├── Core Analysis Scripts
│   ├── read_data_and_basic_corrections.R          [Data import & validation]
│   ├── build_analysis_dataset.R                   [Dataset construction: 242 variables across 8 categories]
│   │
│   ├── trigeminal_measures_distributions_correlations_y_trig_1.R
│   │   └─ Exploratory validation: distributions, correlations, group comparisons on non-imputed data
│   │
│   ├── trigeminal_measures_clustering_y_trig_2_3.R
│   │   └─ Phenotype discovery: unsupervised clustering (6 projections × 8 algorithms)
│   │   └─ Training/validation evaluation of cluster stability
│   │
│   ├── ProjectionsBiomed_MainFunctions_6_1core.R  [Helper functions]
│   └── globals.R                                   [Global variables]
│
├── Exploratory Analyses (Independent)
│   ├── descriptive_statistics_overview.R
│   ├── smoking_overview.R
│   ├── facial_pain_overview.R
│   ├── chronic_diseases_overview.R
│   ├── covid_overview.R
│   └── nasal_breathing_and_ENT_surgery_overview.R
│
├── Data Files
│   ├── *.csv                     [Raw and processed datasets]
│   └── projectionsAndPlots_*.RData [Cached analysis results]
│
└── TrigeminalSensitivity.Rproj

📂 Python/
├── trigeminal_measures_exploreTukey.py
└── [explore_tukey_lop.py](https://github.com/JornLotsch/PythonDataExploration/blob/main/explore_tukey_lop.py)
```

---

## Core Analyses

### Exploratory Validation (Non-Imputed Data)

**File:** `trigeminal_measures_distributions_correlations_y_trig_1.R`

Exploratory validation on **non-imputed data** to document observed patterns without imputation assumptions:

- **Measurement completeness:** AmmoLa (100%), Lateralization (45.8%), CO₂ threshold (45.2%)
- **A priori high-sensitivity classification:** AmmoLa ≥90th percentile (≥90 on 0-100 scale); CO₂ ≤10th percentile
- **Distribution assessment:** Pareto Density Estimation, Gaussian Mixture Modeling (1-4 modes)
- **Agreement analysis:** Fisher's exact tests for concordance across high-sensitivity classifications
- **Univariate comparisons:** Wilcoxon-Mann-Whitney U tests with **rank-biserial correlation effect sizes** (95% bootstrap CI)
- **Correlation structure:** Spearman correlations with bootstrap confidence intervals (1000 replicates)
- **Age/sex effects:** Spearman correlation and Kruskal-Wallis tests

**Key outputs:** Distribution plots, agreement tables, effect size estimates, correlation matrices

### Structure Exploration (Full Imputed Data)

**File:** `trigeminal_measures_distributions_correlations_y_trig_1.R`

Structure exploration on **complete imputed dataset** (n=1001) to characterize population-level trigeminal dimensionality:

- **Input matrix:** 18 TriFunQ items + 3 transformed psychophysical measures
- **Correlation analysis:** Spearman correlations with bootstrap confidence intervals (1000 replicates) across all 153 variable pairs
- **Principal Component Analysis:** Kaiser criterion (eigenvalue > 1), variable loadings, variance explained
- **Effect size interpretation:** Cohen's guidelines (small <0.3, moderate 0.3-0.5, large ≥0.5)

**Key outputs:** PCA results, correlation heatmaps with bootstrap CI, scree plots, variable loading tables

### Cluster-Based Phenotypes (Training/Validation Data)

**File:** `trigeminal_measures_clustering_y_trig_2_3.R`

Data-driven phenotype discovery using unsupervised clustering on **training data only** (n=801) to prevent information leakage:

**Clustering Strategy:**
- **Input matrix:** 18 TriFunQ items + 3 transformed psychophysical measures, scaled to [0,3] and standardized
- **Projection methods** (6): No projection, PCA, ICA, MDS, t-SNE, UMAP
- **Clustering algorithms** (8): k-means, k-medoids, Ward, single/complete/average/median/centroid linkage
- **Total combinations:** 48 (6 projections × 8 algorithms)

**Model Selection:**
- **Optimal cluster number:** NbClust (30 validity indices with majority voting)
- **Quality metrics:** Silhouette Index, Dunn Index, DBCV, Calinski-Harabasz, Davies-Bouldin, inertia
- **Best solution:** Rank product aggregation across all metrics

**Validation Set Assignment:**
- Participants (n=200) assigned to training-derived clusters based on proximity to cluster centers
- No validation data used in cluster definition

**Cluster Characterization:**
- Wilcoxon-Mann-Whitney U tests (2 clusters) or Kruskal-Wallis tests (>2 clusters)
- Effect sizes: **Rank-biserial correlation** for all comparisons
- Bootstrap confidence intervals: 1000 iterations

**Cross-Dataset Validation:**
- Kendall's tau correlation of rank-biserial effect sizes between training and validation sets
- High positive τ indicates stable, reproducible cluster structure

**Phenotype Outcomes:**
- **Y_trig1 (discrete):** Cluster membership from optimal clustering solution
- **Y_trig2 (continuous):** Row means across all 21 trigeminal variables, assessed via Gaussian Mixture Models (1-4 modes)

**Key outputs:** Cluster assignments, quality metrics (dfClusterQuality.csv), silhouette/dendrogram plots, effect size estimates with bootstrap CI, heatmaps

---

## Statistical Methods

**Distribution Modeling:**
- Pareto Density Estimation (PDE) for non-parametric density visualization
- Gaussian Mixture Models (opGMMassessment) with likelihood ratio model selection
- Anderson-Darling test for distribution comparisons

**Hypothesis Testing:**
- Wilcoxon rank-sum test (2 groups), Kruskal-Wallis test (>2 groups)
- Fisher's exact test (contingency tables)
- Spearman rank correlation

**Effect Size Estimation:**
- **Rank-biserial correlation:** Non-parametric effect size for all group comparisons (interpretable as probability of superiority)
- Bootstrap confidence intervals: 1000 iterations using percentile method
- Cohen's guidelines for effect size interpretation

**Dimensionality Reduction:**
- Principal Component Analysis (PCA) with Kaiser criterion (eigenvalue > 1)
- Projection methods: ICA, MDS, t-SNE, UMAP for clustering visualization

**Clustering & Validation:**
- Multiple quality indices: Silhouette, Dunn, DBCV, Calinski-Harabasz, Davies-Bouldin, inertia
- Rank product aggregation for robust solution selection
- Cross-dataset validation: Kendall's tau correlation of effect sizes

---

## Data Processing

**Data Input & Corrections:** [read_data_and_basic_corrections.R](R/read_data_and_basic_corrections.R)
- Raw data import and format validation
- Identifier correction and standardization

**Dataset Construction:** [build_analysis_dataset.R](R/build_analysis_dataset.R)
- Variable selection and categorization across 8 categories:
  - Category 1: Demographics (age, weight, height, sex)
  - Category 2 - 8: Other, see paper
- **Result: 242 input variables**
- Missing data assessment: CO₂ threshold (66.5% missing), Lateralization (45.3% missing)
- Imputation: missRanger (random forest with predictive mean matching, k=5)
- Train/validation split: 80% training (n=801), 20% validation (n=200) using stratified downsampling (Anderson-Darling criterion)
- Post-imputation: Zero missing values across all variables

**Variable Transformations:**
- **AmmoLa intensity**: Reflected sign-preserving logarithmic transformation, rescaled to [0,3]
- **CO₂ threshold**: Sign-preserving logarithmic transformation (negative to preserve sensitivity direction), rescaled to [0,3]
- **Lateralization**: No transformation (approximately normal), rescaled to [0,3]

---

## Usage

**Step 1: Data Preparation (required)**
```r
source("R/read_data_and_basic_corrections.R")
source("R/build_analysis_dataset.R")
```

**Step 2: Univariate Analysis (Y_trig1)**
```r
source("R/trigeminal_measures_distributions_correlations_y_trig_1.R")
```

**Step 3: Clustering Analysis (Y_trig2, Y_trig3)**
```r
source("R/trigeminal_measures_clustering_y_trig_2_3.R")
```

**Optional: Exploratory Analyses**
```r
source("R/descriptive_statistics_overview.R")       # Cohort characteristics
source("R/smoking_overview.R")                      # Smoking-related analyses
source("R/facial_pain_overview.R")                  # Facial pain patterns
source("R/chronic_diseases_overview.R")             # Disease comorbidity
source("R/covid_overview.R")                        # COVID-19 effects
source("R/nasal_breathing_and_ENT_surgery_overview.R")  # Nasal factors
```

---

## Copyright

© 2025 Computational Analysis  
[Code and documentation](https://github.com/JornLotsch/TrigeminalChemosensitivity)

*This repository contains analysis code for scientific research. Use restricted to authorized collaborators and publication in peer-reviewed venues.*
