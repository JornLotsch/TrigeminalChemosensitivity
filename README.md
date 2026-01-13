# Trigeminal Chemosensitivity Analysis

Computational pipeline for analyzing trigeminal sensitivity phenotypes and their relationship to health factors, behavioral patterns, and demographics.

<!-- <img src="./wip-figure.png " width="30%"> -->

## Table of Contents

1. [Computational Setup](#computational-setup)
2. [Project Overview](#project-overview)
3. [Repository Structure](#repository-structure)
4. [Analysis Design](#analysis-design)
   - [Input and Output Spaces](#input-and-output-spaces)
5. [Construction of Trigeminal Sensitivity Phenotypes](#construction-of-trigeminal-sensitivity-phenotypes)
   - [Trigeminal Sensitivity Measurements](#trigeminal-sensitivity-measurements)
   - [Univariate Phenotype (Y_trig1)](#univariate-phenotype-y_trig1)
   - [Composite Phenotypes (Y_trig2 & Y_trig3)](#composite-phenotypes-y_trig2--y_trig3)
   - [Data Preprocessing](#data-preprocessing)
6. [Processing of Modulatory Factors](#processing-of-modulatory-factors)
   - [Feature Categories and Variable Construction](#feature-categories-and-variable-construction)
   - [Missing Data Handling](#missing-data-handling)
7. [Statistical Methods](#statistical-methods)
8. [Study Cohort Characteristics](#study-cohort-characteristics)
9. [Usage](#usage)

---

## Environment

**Language & Platform:**
- R (Linux environment)
- Python
- Tested on Linux; compatibility with other OS not guaranteed

**Key R Packages:**
- Data processing: `dplyr`, `tidyr`
- Visualization: `ggplot2`
- Imputation: `missForest`
- Dimensionality reduction: `FactoMineR`, `fastICA`, `Rtsne`, `umap`
- Clustering: `cluster`, `NbClust`
- Statistical analysis: `MASS`, `psych`

## Project Overview

This pipeline analyzes modulatory effects of health, behavioral, and demographic factors on trigeminal sensitivity. Three complementary phenotypes are constructed from psychophysical measurements and self-report data.

**Phenotypes:**
- **Y_trig1**: Binary classification based on sensitivity threshold
- **Y_trig2**: Continuous multivariate metric
- **Y_trig3**: Data-driven categorical phenotype

**Modulatory factors examined:**
- Health: chronic diseases, ENT surgeries, facial pain, nasal breathing
- Behavioral: smoking, alcohol use
- Temporal: COVID-19 infection history
- Demographics: age, sex, anthropometric measures
## Repository Structure

```
📂 R/
├── Data Preparation
│   ├── read_data_and_basic_corrections.R
│   ├── build_analysis_dataset.R
│   └── globals.R
│
├── Exploratory & Descriptive Analysis
│   ├── descriptive_statistics_overview.R
│   ├── smoking_overview.R
│   ├── facial_pain_overview.R
│   ├── chronic_diseases_overview.R
│   ├── covid_overview.R
│   └── nasal_breathing_and_ENT_surgery_overview.R
│
├── Statistical Analysis
│   ├── trigeminal_measures_distributions_correlations.R
│   ├── ammo_distribution_analysis.R
│   └── ProjectionsBiomed_MainFunctions_6_1core.R
│
└── TrigeminalSensitivity.Rproj

📂 Python/
├── trigeminal_measures_exploreTukey.py
└── explore_tukey_lop.py
```

## Analysis Design

**Output (Y):** Trigeminal sensitivity phenotypes constructed from psychophysical measurements and questionnaire responses

**Input (X):** Modulatory factors including demographics, health history, behavioral factors, and temporal variables

---

## Construction of Trigeminal Sensitivity Phenotypes

### Overview

Three phenotypes are derived from psychophysical measurements and questionnaire responses.

**Data Preprocessing:** Missing values are handled via random forest imputation. Distributions are tested for normality; transformations are applied where needed.

### Phenotype 1: Binary Classification (Y_trig1)

**[trigeminal_measures_distributions_correlations.R](R/trigeminal_measures_distributions_correlations.R)**

Binary classification based on sensitivity threshold. Validation includes correlation analysis, agreement testing, and group comparisons using appropriate statistical tests.

### Phenotype 2: Continuous Multivariate Metric (Y_trig2)

**[trigeminal_measures_distributions_correlations.R](R/trigeminal_measures_distributions_correlations.R) & [ammo_distribution_analysis.R](R/ammo_distribution_analysis.R)**

Derived from questionnaire responses and transformed psychophysical measurements. Distribution tested via Gaussian Mixture Modeling with model selection based on likelihood ratio tests. Principal Component Analysis applied for dimension examination.

---

## Processing of Modulatory Factors

**[read_data_and_basic_corrections.R](R/read_data_and_basic_corrections.R) & [build_analysis_dataset.R](R/build_analysis_dataset.R)**

Data processing pipeline:
1. Data import and validation
2. Error correction and standardization
3. Feature engineering and transformation
4. Handling of categorical and temporal variables

**Factor categories:**
- Demographics (age, sex, anthropometric measures)
- Health: chronic diseases, ENT surgeries, facial pain, nasal breathing
- COVID-19: infection status, temporal patterns
- Lifestyle: smoking status, pack-years, alcohol use
- Olfactory and nasal variables

**Missing data:** Random forest imputation for low-level missing data; variables with >20% missing excluded.

## Statistical Methods

**Univariate Testing:**
- Correlation: Spearman, Pearson
- Group comparison: Kruskal-Wallis, Wilcoxon-Mann-Whitney U, Fisher's exact
- Effect sizes: Cohen's d with bootstrap confidence intervals

**Dimensionality Reduction:**
- PCA, ICA, MDS, t-SNE, UMAP

**Clustering:**
- Partitioning: k-means, PAM
- Hierarchical: Ward, single, average, median, complete, centroid linkage
- Validation: NbClust (30 validity indices)

**Distribution Modeling:**
- Gaussian Mixture Modeling (1-4 modes)
- Likelihood ratio model selection

---

## Usage

**Data Preparation (required first):**
```r
source("R/read_data_and_basic_corrections.R")
source("R/build_analysis_dataset.R")
```

**Exploratory Analyses (independent):**
```r
source("R/descriptive_statistics_overview.R")
source("R/smoking_overview.R")
source("R/facial_pain_overview.R")
source("R/chronic_diseases_overview.R")
source("R/covid_overview.R")
source("R/nasal_breathing_and_ENT_surgery_overview.R")
```

**Statistical Analysis:**
```r
source("R/trigeminal_measures_distributions_correlations.R")
source("R/ammo_distribution_analysis.R")
```

---

## License

Analysis code and documentation available at: https://github.com/JornLotsch/TrigeminalChemosensitivity
