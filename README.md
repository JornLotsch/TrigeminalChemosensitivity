# Trigeminal Chemosensitivity Analysis

Computational pipeline for characterizing trigeminal sensitivity phenotypes and identifying their modulatory factors. This repository accompanies a peer-reviewed scientific manuscript describing a two-stage analysis: (1) characterization of the trigeminal data space and derivation of a composite trigeminal phenotype via unsupervised clustering, and (2) identification of modulators of both the phenotype and three psychophysical measures (AmmoLa intensity, CO₂ threshold, Lateralization).

**Development environment:** R 4.5.3 and 4.6 for Linux

---

## Analysis Design

**Two complementary analytical levels:**
1. **Internal trigeminal data structure**: 21 variables (3 psychophysical measures + 18 TriFunQ items)
2. **External modulatory factors**: 220 non-trigeminal candidate predictors (demographic, medical, lifestyle, olfactory)

**Two-part workflow:**

**Part 1: Exploratory Analysis (Complete Imputed Dataset, n=1,001)**
- Correlation analysis and PCA of trigeminal variables
- Regression screening of 220 candidate modulators against 3 psychophysical measures
- Descriptive, hypothesis-generating; no held-out validation
- *Files:* [01_stage1_trigeminal_data_space_exploration.R](R/01_stage1_trigeminal_data_space_exploration.R)

**Part 2: Machine Learning Framework (Train/Validation Split, 80%/20%)**
- All model development confined to training set (n=800)
- Derivation of composite trigeminal phenotype via unsupervised clustering
- Supervised models: Predict psychophysical measures and cluster membership
- Generalization: Evaluated exclusively on validation set (n=201)
- *Files:* [02_stage1_trigeminal_clustering_phenotype_derivation.R](R/02_stage1_trigeminal_clustering_phenotype_derivation.R), [03_stage1_trigeminal_cluster_classification.R](R/03_stage1_trigeminal_cluster_classification.R), [05_psychophysical_regression_analysis.R](R/05_psychophysical_regression_analysis.R), [06_stage2_modulator_identification_psychophysical.R](R/06_stage2_modulator_identification_psychophysical.R), [07_stage2_modulator_identification_clusters.R](R/07_stage2_modulator_identification_clusters.R)

**Critical design feature:** Train/validation split executed *before* imputation to prevent information leakage from held-out participants.

See manuscript Methods for full technical details.

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
- Imputation: `missRanger` (random forest-based with predictive mean matching)
- Dimensionality reduction: `FactoMineR`, `factoextra`, `EDOtrans`, `fastICA`, `Rtsne`, `uwot`
- Clustering: `cluster`, `NbClust`, `dendextend`
- Statistical analysis: `psych`, `MASS`, `boot`, `effsize`, `car`, `glmnet`, `randomForest`, `Boruta`, `nnet`

## Repository Structure

```
📂 R/
├── Core Infrastructure
│   ├── utils.R                                     [General utility functions: 25 helpers]
│   ├── globals.R                                   [Global config: libraries, colors, labels]
│   └── ProjectionsBiomed_MainFunctions_6_1core.R   [Helper library for projection methods]
│
├── Data Processing Pipeline (Required)
│   ├── read_data_and_basic_corrections.R             [Data import & validation]
│   |── build_analysis_dataset.R                      [Dataset construction: 242 variables]
│   ├── smoking_overview.R                            [Smoking patterns]
│   ├── facial_pain_overview.R                        [Facial pain]
│   ├── chronic_diseases_overview.R                   [Disease comorbidity]
│   ├── covid_overview.R                              [COVID-19 effects]
│   └── nasal_breathing_and_ENT_surgery_overview.R    [Nasal/ENT factors]
|
│
├── *** REPORTED ANALYSES (Manuscript) ***
│
│   ├── Data Preparation
│   │   ├── read_data_and_basic_corrections.R       [Data import & validation]
│   │   └── build_analysis_dataset.R                [Dataset construction: 241 variables; sample filtering]
│   │
│   ├── PART 1: EXPLORATORY ANALYSIS (Complete Imputed Dataset, n=1,001)
│   │   │
│   │   └── 01_stage1_trigeminal_data_space_exploration.R
│   │       ├─ Spearman rank correlations (210 pairs, 21 variables)
│   │       ├─ PCA of trigeminal variables (Kaiser-Guttman: eigenvalue > 1)
│   │       ├─ ABC analysis for variable importance ranking
│   │       └─ Exploratory regression: 220 candidate modulators → 3 psychophysical measures
│   │          (Hypothesis-generating; no validation; see Methods for multiple testing context)
│   │
│   ├── PART 2: MACHINE LEARNING FRAMEWORK (Train/Validation Split: 80%/20%, n=800/201)
│   │   │
│   │   ├── 02_stage1_trigeminal_clustering_phenotype_derivation.R
│   │   │   ├─ Training data only: 48 clustering solutions (6 projections × 8 algorithms)
│   │   │   ├─ Projection methods: EDOtrans, PCA, ICA, MDS, t-SNE, UMAP
│   │   │   ├─ Clustering algorithms: k-means, PAM, hierarchical (6 linkages)
│   │   │   ├─ Optimal clusters: NbClust (30 validity indices, majority vote)
│   │   │   ├─ Quality metric: Calinski-Harabasz index; null model comparison
│   │   │   └─ Validation: Nearest-centroid assignment; silhouette width
│   │   │
│   │   ├── 03_stage1_trigeminal_cluster_classification.R
│   │   │   ├─ Supervised classification: Cluster membership (categorical outcome)
│   │   │   ├─ Predictors: 20 remaining trigeminal variables
│   │   │   ├─ Models: Multinomial logistic, penalized variants, random forest
│   │   │   ├─ Training: Training set only
│   │   │   └─ Validation: Held-out validation set; balanced accuracy with 95% CI
│   │   │
│   │   ├── 04_stage1_trigeminal_cluster_interpretation.R
│   │   │   ├─ Cluster characterization: Kruskal-Wallis tests (effect size η²; 95% CI)
│   │   │   ├─ Psychometric interpretation: TriFunQ dimensions (sensitivity + avoidance)
│   │   │   └─ Cross-dataset concordance: Kendall's τ (training vs. validation)
│   │   │
│   │   ├── 05_psychophysical_regression_analysis.R
│   │   │   ├─ Continuous outcomes: 3 psychophysical measures (transformed)
│   │   │   ├─ Predictors: 20 remaining trigeminal variables
│   │   │   ├─ Models: OLS, ridge, lasso, elastic net (α=0.5), random forest
│   │   │   ├─ Training: Training set; λ tuned by 5-fold CV
│   │   │   ├─ Random forest: mtry & ntree optimized by grid search
│   │   │   └─ Validation: Nonlinear LS regression (slope & significance test)
│   │   │
│   │   ├── 06_stage2_modulator_identification_psychophysical.R
│   │   │   ├─ Same regression framework as Step 5
│   │   │   ├─ Predictors: 220 non-trigeminal candidate modulators
│   │   │   ├─ Outcomes: 3 psychophysical measures
│   │   │   ├─ Feature importance: Lasso selection + Boruta (≤1,000 iterations)
│   │   │   └─ Concordance: Wilcoxon/Spearman with Kendall's τ across datasets
│   │   │
│   │   └── 07_stage2_modulator_identification_clusters.R
│   │       ├─ Multinomial classification: Cluster membership (categorical outcome)
│   │       ├─ Predictors: 220 non-trigeminal candidate modulators
│   │       ├─ Models: Multinomial logistic, penalized variants, random forest
│   │       ├─ Feature importance: Boruta with reduced feature set
│   │       └─ Concordance: Wilcoxon/Spearman with Kendall's τ across datasets
│

│
├── *** EXPLORATORY ANALYSES ***
│   ├── descriptive_statistics_overview.R             [Cohort characteristics]
│
└── Data Files [Not version-controlled: raw and processed datasets stored locally]

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
- CO₂ sample filtering (breath-hold protocol only; n=453)
- Variable recoding: 241 variables
- **Train/validation split (80%/20%, n=800/201) executed BEFORE imputation** (critical for preventing information leakage)
- Imputation: missRanger (random forest with PMM k=5; 500 trees; applied independently to training and validation sets)
- Distribution assessment: Tukey's ladder of powers; transformations applied to AmmoLa and CO₂
- **Output: 21 trigeminal variables** (18 TriFunQ items + 3 transformed psychophysical measures)

---

### PART 1: Exploratory Analysis (Complete Imputed Dataset, n=1,001)

**Objective:** Characterize internal trigeminal data structure and screen for candidate modulators using hypothesis-generating methods without held-out validation.

#### Trigeminal Data Space Exploration

**File:** [01_stage1_trigeminal_data_space_exploration.R](R/01_stage1_trigeminal_data_space_exploration.R)

- **Correlations:** Spearman rank (210 unique pairs; 95% CI via bootstrapping, 1,000 iterations)
- **PCA:** Standardized 1,001 × 21 matrix; Kaiser-Guttman criterion (eigenvalue > 1)
- **Variable importance:** ABC analysis (Pareto-based) on absolute standardized loadings; variables ranked into A/B/C subsets
- **Exploratory screening:** Standard linear and 3 penalized approaches (ridge, lasso, elastic net; λ tuned by 5-fold CV) examining associations between 220 candidate modulators and each of 3 psychophysical measures
  - *Note:* These are descriptive association screens; interpretation requires manuscript Methods context on multiple testing framework

---

### PART 2: Machine Learning Framework (Train/Validation Split: n=800 training, n=201 validation)

**Objective:** Derive composite trigeminal phenotype via unsupervised clustering (training set only) and evaluate supervised predictive models on independent validation data.

#### Derivation of Composite Trigeminal Phenotype

**File:** [02_stage1_trigeminal_clustering_phenotype_derivation.R](R/02_stage1_trigeminal_clustering_phenotype_derivation.R)

**Training set only (800 × 21 matrix):**
- **48 clustering solutions** systematically compared:
  - 6 projection methods: EDOtrans, PCA, ICA, MDS, t-SNE, UMAP
  - 8 clustering algorithms: k-means, PAM, hierarchical (complete, single, average, weighted, centroid, McQuitty, median)
  - Dimensionality retained per method-specific criteria
- **Optimal cluster number:** NbClust with 30 validity indices; majority voting
- **Quality assessment:** Calinski-Harabasz index (primary criterion); null model comparison (3 permuted datasets)
- **Validation set:** Participants projected to fitted training space; nearest-centroid assignment; silhouette width quantification

#### Supervised Classification of Trigeminal Cluster Membership

**File:** [03_stage1_trigeminal_cluster_classification.R](R/03_stage1_trigeminal_cluster_classification.R)

- **Outcome:** Cluster membership (categorical)
- **Predictors:** 20 remaining trigeminal variables (excluding 3 psychophysical outcomes used in clustering)
- **Models:** Multinomial logistic, ridge, lasso, elastic net, random forest
- **Hyperparameter tuning:** λ by 5-fold CV (penalized); mtry & ntree by grid search (RF)
- **Training:** Training set only
- **Validation:** Balanced accuracy with 95% bootstrap CI on held-out validation set

#### Cluster Interpretation

**File:** [04_stage1_trigeminal_cluster_interpretation.R](R/04_stage1_trigeminal_cluster_interpretation.R)

- **Trigeminal variable comparisons:** Kruskal-Wallis tests; effect sizes (η²) with 95% bootstrap CI
- **Psychometric interpretation:** TriFunQ items a priori classified into sensitivity and avoidance dimensions
- **Cross-dataset validation:** Effect size concordance (training vs. validation) assessed by Kendall's τ

#### Supervised Prediction of Psychophysical Measures (from trigeminal variables)

**File:** [05_psychophysical_regression_analysis.R](R/05_psychophysical_regression_analysis.R)

- **Outcomes:** 3 transformed psychophysical measures (continuous)
- **Predictors:** 20 remaining trigeminal variables
- **Models:** OLS, ridge, lasso, elastic net (α=0.5), random forest
- **Training:** Training set; λ tuned by 5-fold CV; RF: mtry & ntree by grid search
- **Validation:** Nonlinear LS regression (predicted vs. observed) with slope significance test

#### Supervised Prediction of Psychophysical Measures (from candidate modulators)

**File:** [06_stage2_modulator_identification_psychophysical.R](R/06_stage2_modulator_identification_psychophysical.R)

- **Same regression framework as Step 5**
- **Predictors:** 220 non-trigeminal candidate modulators (demographic, medical, lifestyle, olfactory)
- **Outcomes:** 3 psychophysical measures
- **Feature importance:** Lasso-based variable selection; Boruta algorithm (≤1,000 iterations)
- **Cross-dataset validation:** For models generalizing successfully, selected variables characterized independently in training & validation sets (Wilcoxon/Spearman); effect size concordance via Kendall's τ

#### Supervised Classification of Cluster Membership (from candidate modulators)

**File:** [07_stage2_modulator_identification_clusters.R](R/07_stage2_modulator_identification_clusters.R)

- **Outcome:** Cluster membership (categorical)
- **Predictors:** 220 non-trigeminal candidate modulators
- **Models:** Multinomial logistic, ridge, lasso, elastic net, random forest (Boruta-reduced feature set)
- **Validation:** Balanced accuracy with 95% bootstrap CI on held-out data
- **Cross-dataset characterization:** Selected variables assessed independently in training & validation sets (effect size concordance via Kendall's τ)

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

### Sample Inclusion & Protocol
- CO₂ perception thresholds available for 453/1,001 participants
- Uncontrolled breathing protocol (n=117): 17.1% ceiling effects
- Breath-hold protocol (n=336): 5.4% ceiling effects
- **Data from uncontrolled protocol excluded; analysis uses controlled breath-hold protocol only**

### Variable Recoding & Construction
- Continuous demographics (age, weight, height): retained unchanged
- Categorical variables: one-hot encoded (gender) and manually standardized (diseases, surgeries, smoking, COVID-19)
- TriFunQ items (15) and onion-related variables (3): assigned to trigeminal output space
- Olfactory variables: continuous or recoded ordinal; odor identification untransformed
- Nasal irritation/airflow: included without transformation
- **Output: 241 variables × 1,001 participants; 1.75% missing (4,166 values)**

### Distribution Assessment & Transformations
- **Method:** Tukey's ladder of powers (λ = -2, -1, -0.5, 0, 0.5, 1, 2) with Box-Cox reference; D'Agostino's K² normality test
- **AmmoLa intensity** (left-skewed, ceiling effects): Reflected sign-preserving log transformation
  - Formula: 
```math
- \operatorname{sign}\bigl(\max(\mathrm{AmmoLa}) + 1 - \mathrm{AmmoLa}\bigr)
\cdot \log_{10}\bigl(\lvert \max(\mathrm{AmmoLa}) + 1 - \mathrm{AmmoLa}\rvert + 1\bigr)
```
- **CO₂ threshold** (right-skewed; lower = higher sensitivity): Sign-inverted log transformation  
  - Formula: 
```math
- \operatorname{sign}(\mathrm{CO}_2) \times \log_{10}\bigl(\lvert \mathrm{CO}_2\rvert + 1\bigr)
```
- **Lateralization** (approximately normal): Retained untransformed

### Dataset Split & Imputation
- **Split before imputation** (critical): Training (80%, n=800) and validation (20%, n=201) sets created first
- **Imputation method:** missRanger (random forest-based; https://cran.r-project.org/package=missRanger)
  - Predictive mean matching (k=5 nearest donors)
  - Configuration: 500 trees, minimum node size = 5
  - Type integrity enforced post-imputation (binary: [0,1]; ordinal: nearest integer; continuous: not rounded)
- **Imputation fidelity:** Kernel-density overlap and Anderson-Darling tests (p=1.0) confirmed no distributional differences for CO₂ and lateralization
- **Result:** Two independently imputed datasets (training and validation) used for all subsequent analyses

### Variable Transformations (Summary)
See [utils.R](R/utils.R) for transformation implementations.

---

## Usage

**Code Organization:**
All analysis scripts automatically load required utilities through the following dependency chain:
```
utils.R → globals.R → analysis scripts
```

**To Reproduce the Manuscript Analyses (in order):**

---

### Data Preparation (Required)

```r
source("R/read_data_and_basic_corrections.R")
source("R/build_analysis_dataset.R")
```
**Output:** Train/validation split (n=800/201), independently imputed 21-variable trigeminal matrices

---

### PART 1: Exploratory Analysis (Complete Imputed Dataset, n=1,001)

This part is hypothesis-generating and descriptive; uses merged training + validation data.

```r
source("R/01_stage1_trigeminal_data_space_exploration.R")
```
**Output:** Correlations, PCA, variable importance, exploratory modulator screening

---

### PART 2: Machine Learning Framework (Train/Validation Split)

All model development confined to training set; validation on held-out data.

#### Step 2a: Derive Composite Trigeminal Phenotype
```r
source("R/02_stage1_trigeminal_clustering_phenotype_derivation.R")
```
**Output:** Trigeminal cluster assignments (training set → validated on held-out), quality metrics

#### Step 2b: Supervised Cluster Classification
```r
source("R/03_stage1_trigeminal_cluster_classification.R")
source("R/04_stage1_trigeminal_cluster_interpretation.R")
```
**Output:** Classification models, cluster characterization, psychometric interpretation

#### Step 2c: Regression Models of Psychophysical Measures (trigeminal → psychophysical)
```r
source("R/05_psychophysical_regression_analysis.R")
```
**Output:** Regression models predicting AmmoLa, CO₂ threshold, Lateralization from trigeminal variables

#### Step 2d: Modulator Identification (external candidates → psychophysical & clusters)
```r
source("R/06_stage2_modulator_identification_psychophysical.R")
source("R/07_stage2_modulator_identification_clusters.R")
```
**Output:** Feature importance, modulator characterization, cross-dataset concordance

---

**Optional: Exploratory Population Context**

These supplementary descriptive scripts are not part of the core reported pipeline but provide population-level context:

```r
source("R/descriptive_statistics_overview.R")
source("R/smoking_overview.R")
source("R/facial_pain_overview.R")
source("R/chronic_diseases_overview.R")
source("R/covid_overview.R")
source("R/nasal_breathing_and_ENT_surgery_overview.R")
```

**Note:** All scripts automatically load utilities via `globals.R`. No need to manually source `utils.R` or `ProjectionsBiomed_MainFunctions_6_1core.R`.

---

## Citation

- Lötsch J, Weise S, Himmelspach A, Bormann A, and Hummel T
Phenotypes of trigeminal chemosensitivity in a general-population sample 2026 (in preparation)


## License

[Code and documentation](https://github.com/JornLotsch/TrigeminalChemosensitivity)

CC-BY 4.0

*This repository contains analysis code for scientific research. Use restricted to authorized collaborators and publication in peer-reviewed venues.*
