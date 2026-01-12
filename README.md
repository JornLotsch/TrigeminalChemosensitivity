# Trigeminal Chemosensitivity Analysis

A comprehensive computational study investigating modulatory effects of acquired variables (health factors, behavioral patterns, demographic characteristics) on trigeminal sensitivity phenotypes in a general population sample (n=1,001).

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

## Computational Setup

**Analysis Environment:**
- R version 4.5.2 (Linux)
- Python version 3.12.2
- IDEs: RStudio and PyCharm Professional (version 2025.2.3, JetBrains)
- **Operating System**: All code has been tested and validated on Linux (Ubuntu 24.04.3 LTS) only. Compatibility with macOS and Windows is not guaranteed.

**Primary R Packages:**
- Data manipulation: `dplyr`, `tidyr`, `stringr`
- Data visualization: `ggplot2`, `ggtext`
- Temporal operations: `lubridate`
- Statistical summaries: `psych`

**Specialized Packages:**
- Missing data imputation: `missForest`
- Gaussian Mixture Modeling: `opGMMassessment`
- Dimensionality reduction: `FactoMineR`, `fastICA`, `Rtsne`, `umap`
- Clustering: `cluster`, `NbClust`
- Statistical testing: `MASS`, `DescTools`

## Project Overview

This project investigates potential modulatory effects of health factors, behavioral patterns, and demographics (input space **X**) on trigeminal sensitivity (output space **Y**). Three complementary trigeminal sensitivity phenotypes were constructed from psychophysical measurements and self-reported sensations:

1. **Y_trig1**: HighLow_AmmoLa - univariate psychophysical classification
2. **Y_trig2**: marginal_trigeminal_sensitivity - multivariate continuous metric
3. **Y_trig3**: trigeminal_clusters - data-driven categorical phenotype

Modulatory factors examined include:
- **Health characteristics**: Chronic diseases, ENT surgeries, facial pain, nasal breathing problems
- **Behavioral factors**: Smoking history, pack-years, alcohol consumption
- **Temporal effects**: COVID-19 infection history and temporal patterns
- **Demographic factors**: Age, sex, weight, height
## Repository Structure

```
📂 R/
├── Data Preparation & Preprocessing
│   ├── read_data_and_basic_corrections.R      # Initial data import, validation, corrections
│   ├── build_analysis_dataset.R               # Construction of analysis-ready dataset
│   └── globals.R                              # Global variables and configuration
│
├── Exploratory & Descriptive Analysis
│   ├── descriptive_statistics_overview.R      # Cohort characteristics and summaries
│   ├── smoking_overview.R                     # Smoking patterns and temporal analysis
│   ├── facial_pain_overview.R                 # Facial pain prevalence and characteristics
│   ├── chronic_diseases_overview.R            # Chronic disease prevalence and patterns
│   ├── covid_overview.R                       # COVID-19 infection and temporal patterns
│   └── nasal_breathing_and_ENT_surgery_overview.R  # Nasal issues and ENT procedures
│
├── Advanced Statistical Analysis
│   ├── trigeminal_measures_distributions_correlations.R  # Phenotype construction & validation
│   ├── ammo_distribution_analysis.R           # Distribution analysis and transformations
│   └── ProjectionsBiomed_MainFunctions_6_1core.R  # Utilities for analysis
│
└── TrigeminalSensitivity.Rproj

📂 Python/
├── trigeminal_measures_exploreTukey.py        # Tukey transformation exploration
└── explore_tukey_lop.py                       # Additional transformation analysis

📂 Matlab/                                      # Reserved for future use
```
```
📂 R/
├── Data Preparation & Preprocessing
│   ├── read_data_and_basic_corrections.R      # Initial data import, validation, corrections
│   ├── build_analysis_dataset.R               # Construction of analysis-ready dataset
│   └── globals.R                              # Global variables and configuration
│
├── Exploratory & Descriptive Analysis
│   ├── descriptive_statistics_overview.R      # Cohort characteristics and summaries
│   ├── smoking_overview.R                     # Smoking patterns and temporal analysis
│   ├── facial_pain_overview.R                 # Facial pain prevalence and characteristics
│   ├── chronic_diseases_overview.R            # Chronic disease prevalence and patterns
│   ├── covid_overview.R                       # COVID-19 infection and temporal patterns
│   └── nasal_breathing_and_ENT_surgery_overview.R  # Nasal issues and ENT procedures
│
├── Advanced Statistical Analysis
│   ├── trigeminal_measures_distributions_correlations.R  # Phenotype construction & validation
│   ├── ammo_distribution_analysis.R           # Distribution analysis and transformations
│   └── ProjectionsBiomed_MainFunctions_6_1core.R  # Utilities for analysis
│
└── TrigeminalSensitivity.Rproj

📂 Python/
├── trigeminal_measures_exploreTukey.py        # Tukey transformation exploration
└── explore_tukey_lop.py                       # Additional transformation analysis

📂 Matlab/                                      # Reserved for future use
```

## Analysis Design

### Input and Output Spaces

The analysis framework defines:

- **Output Space (Y)**: Trigeminal sensitivity, represented through three distinct phenotypes constructed from psychophysical measurements (AmmoLa intensity, lateralization, CO₂ threshold) and self-reported trigeminal sensations (TriFunQ questionnaire items)

- **Input Space (X)**: Modulatory factors comprising 247 variables across 7 categories—demographics, chronic disease history, COVID-19 exposure, smoking and alcohol use, olfactory function, nasal irritation, and odor identification performance

---

## Construction of Trigeminal Sensitivity Phenotypes

### Trigeminal Sensitivity Measurements

The study employed three psychophysical measures of trigeminal sensitivity:

1. **Ammonia Lateralization Intensity (AmmoLa)**: Available for all 1,001 participants
2. **Lateralization Score**: Available in 458 participants (non-overlapping subset)
3. **CO₂ Detection Threshold**: Available in 453 participants (non-overlapping subset)
   - Subgroup without controlled breathing (n=177): 17.1% ceiling censoring
   - Subgroup with breath-hold protocol (n=336): 5.4% ceiling censoring (used for analysis)

These measures were combined with 15 TriFunQ questionnaire items (0.47% missing values) assessing self-reported trigeminal sensation.

### Univariate Phenotype (Y_trig1)

**[trigeminal_measures_distributions_correlations.R](R/trigeminal_measures_distributions_correlations.R)**

**Definition: HighLow_AmmoLa**

- Binary classification based on AmmoLa intensity (available for all n=1,001)
- High sensitivity: ≥90th percentile (n=100)
- Low-to-moderate sensitivity: remainder (n=901)

**Validation Analyses:**

1. **Correlation Analysis**
   - Spearman rank correlations between AmmoLa and other trigeminal measures
   - Pearson correlations with age
   - Kruskal-Wallis tests for biological sex associations
   - Result: Generally weak and non-significant correlations

2. **Agreement at Sensitivity Extremes**
   - Fisher's exact tests for 2×2 contingency tables
   - Result: Significant enrichment of highly sensitive individuals between AmmoLa and CO₂ measures (OR=1.99, p=0.014)
   - No significant association with lateralization scores (p=0.81)

3. **Group Comparisons**
   - Wilcoxon-Mann-Whitney U rank-sum tests
   - Effect sizes: Cohen's d with Hedges' correction
   - 95% bootstrap confidence intervals (1,000 iterations)
   - Result: Three TriFunQ items showed small but significant effects; AmmoLa showed expected large effect (d=-2.71)

### Composite Phenotypes (Y_trig2 & Y_trig3)

#### Data Preprocessing

**[ammo_distribution_analysis.R](R/ammo_distribution_analysis.R)**

**AmmoLa Transformation:**
- Problem: Left-skewness and ceiling effects in AmmoLa intensity ratings
- Solution: Reflected sign-log transformation
  $$\text{AmmoLa}_{\text{reflected,slog}} = -\text{sign}(\max(\text{AmmoLa})+1-\text{AmmoLa}) \cdot \log_{10}(|\max(\text{AmmoLa})+1-\text{AmmoLa}|+1)$$
- Rescaling: Transformed values rescaled to [0, 3] range for comparability with TriFunQ items

**Missing Data Handling:**
- 15 TriFunQ variables: 72 missing values (0.47%) imputed using random forests (`missForest`)
- Variables excluded: 3 onion-cutting items with >20% missing values

#### Marginal Trigeminal Sensitivity (Y_trig2)

**[trigeminal_measures_distributions_correlations.R](R/trigeminal_measures_distributions_correlations.R)**

**Definition:**
- Unidimensional phenotype: row means across 16 variables (15 TriFunQ + transformed AmmoLa)
- Distribution characteristics: Normal, unimodal (mean=1.42, SD=0.43)
- Anderson-Darling test: A=0.983, p=0.366 (non-significant, supporting normality)
- Gaussian Mixture Modeling: Single mode optimal (likelihood ratio tests)

**Principal Component Analysis:**
- PCA retained 4 components (eigenvalue >1): explained 43.6% of variance
- >90% variance requires 14 components (indicating complex data structure)
- PC1: Variables related to pungent odors
- PC2: AmmoLa intensity, minty/onion sensations, nasal breathing questions
- No clear discrete groups observed in PC1-PC2 projection

#### Cluster-Based Trigeminal Phenotype (Y_trig3)

**Projection Methods & Clustering Evaluation:**

Systematic evaluation across multiple dimensionality reduction techniques:
- **Projection methods** (6): No projection, PCA, ICA, MDS, t-SNE, UMAP
- **Clustering algorithms** (8): k-means, PAM, hierarchical (Ward, single, average, median, complete, centroid linkage)
- **Cluster determination**: NbClust package (30 validity indices, majority voting)

**Optimal Solution:**
- **Method**: UMAP projection + hierarchical centroid linkage
- **Result**: 2 clusters (Cluster 1: n=551; Cluster 2: n=450)
- **Quality metrics** (6): Silhouette index, Dunn index, Davies-Bouldin index, Calinski-Harabasz index, DBCV index, within-cluster inertia

**Cluster Characterization:**
- **Primary distinction**: Sensitivity to pungent odors (large effect sizes: d>1.1)
- **Item examples**: "I avoid burning or pungent smells" (d=1.25, p<10⁻²⁰); "When I smell something biting or pungent, I panic" (d=1.14, p<10⁻²⁰)
- **Secondary features**: Nasal breathing and discomfort (moderate-to-small effects)
- **Psychophysical measures**: AmmoLa intensity differed significantly (p=5.4×10⁻⁵, d=-0.38); lateralization and CO₂ thresholds showed minimal separation

---

## Processing of Modulatory Factors

### Feature Categories and Variable Construction

**[read_data_and_basic_corrections.R](R/read_data_and_basic_corrections.R) & [build_analysis_dataset.R](R/build_analysis_dataset.R)**

Data preprocessing pipeline:
1. Excel data import from Bormann Trigeminale Studie Daten.xlsx
2. Data validation, error correction, and typographical error fixing
3. German-to-English translation of categorical variables
4. Temporal variable extraction from free-text responses
5. Feature engineering and variable transformation
6. One-hot encoding of categorical variables

#### Category 1: Demographics (n=7)
- Age, weight, height (numeric)
- Gender: one-hot encoded (male, female, diverse, missing indicator)

#### Category 2: Disorders or Health Complaints (n=180)
- **Chronic diseases** (one-hot encoded): Hypertension, allergies, migraine, neurological disorders, chronic sinusitis, diabetes, asthma, hypothyroidism, atopic dermatitis, depression, epilepsy, MS, Hashimoto's thyroiditis, rheumatism, polyneuropathy, ADHD, allergic asthma, Crohn's disease, etc.
- **ENT surgeries** (one-hot encoded): Adenotomy, tonsillotomy, nasal polyps surgery, nasal septum surgery, tonsillectomy, inferior turbinate cauterization, otoplasty, etc.
- **Facial pain**: Presence (0=no, 1=past, 2=current), frequency categories, 4 quality descriptors (pulling, stabbing, pressing, burning)
- **Nasal breathing**: Medical consultation status, ongoing therapy status

#### Category 3: COVID-19 History (n=13)
- Infection status (binary) and infection count (1-5 episodes)
- Olfactory function: pre-infection and immediate post-infection ratings
- Time since infection: months for up to 3 separate episodes (longest and shortest intervals)
- Smell reduction: binary indicators for each infection period

#### Category 4: Smoking and Alcohol Use (n=9)
- **Smoking status**: Current, former, never-smoker (derived from multiple questions)
- **Temporal variables**: Extracted from free text → numerical formats (years/months)
- **Pack-years**: Calculated from daily consumption × duration
- **Cessation timing**: Years since smoking cessation (infinity for never-smokers)
- **Alcohol consumption**: Ordinal variable

#### Category 6: Rated Olfactory Function (n=12)
- Self-rated current smell ability (0-100% continuous)
- Problem onset: ordinal (gradual, sudden, uncertain)
- Temporal change: ordinal (no change, improved, worsened)
- Affected modalities: one-hot encoded (smell only, taste only, both, unknown)

#### Category 7: Nasal Irritation and Airflow (n=4)
- Sensitivity to stinging/burning stimuli (1-4 scale)
- Self-rated airflow: combined nostrils, right nostril, left nostril (0-100 scale)

#### Category 8: Psychophysical Measurements (n=1)
- Odor identification test: Score range 0-3 correct

### Missing Data Handling

**Exclusion Threshold:**
- Variables with >20% missing: excluded (7 variables removed)
  - Migraine-related variables (2): frequency per month, changes over 10 years
  - Onion-cutting variables (3): for trigeminal assessment
  - Psychophysical measures (2): Lateralization, CO₂ threshold

**Low-Level Missing Data:**
- After exclusion: 293 missing values across 240 variables (0.12% of total)
- Imputation method: Random forest-based `missForest` algorithm

**Phenotype-Specific Data Gaps:**
- Lateralization and CO₂ thresholds: subset-specific analyses only (non-overlapping subgroups)

## Statistical Methods

**Univariate Statistical Tests:**
- Spearman rank correlations (monotonic relationships)
- Pearson correlations (age associations)
- Kruskal-Wallis tests (biological sex associations, eta-squared effect sizes)
- Fisher's exact tests (2×2 contingency tables, sensitivity concordance)
- Wilcoxon-Mann-Whitney U rank-sum tests (distribution comparisons)
- Effect sizes: Cohen's d with Hedges' correction, 95% bootstrap confidence intervals (1,000 iterations)
- Anderson-Darling normality tests

**Dimensionality Reduction Methods:**
- **Principal Component Analysis (PCA)**: Kaiser-Guttman criterion (eigenvalue > 1) for component retention
- **Independent Component Analysis (ICA)**: Implementation via `fastICA`
- **Multidimensional Scaling (MDS)**: Via `MASS` package
- **t-Distributed Stochastic Neighbor Embedding (t-SNE)**: Via `Rtsne`
- **Uniform Manifold Approximation and Projection (UMAP)**: Via `umap`

**Clustering Algorithms:**
- **Partitioning methods**: k-means, Partitioning Around Medoids (PAM)
- **Hierarchical methods**: Ward linkage, single, average, median, complete, and centroid linkage
- **Cluster number determination**: NbClust package (30 validity indices, majority voting)
- **Cluster quality metrics**: Silhouette index, Dunn index, Davies-Bouldin index, Calinski-Harabasz index, DBCV index, within-cluster inertia

**Distribution Modeling:**
- **Gaussian Mixture Modeling**: 1-4 modes, evolutionary algorithm optimization via `opGMMassessment`
- **Model selection**: Likelihood ratio tests for comparing models
- **Transformation assessment**: Reflected sign-log, logarithmic variants

**Exploratory Analyses:**

**[descriptive_statistics_overview.R](R/descriptive_statistics_overview.R)**
- Cohort summary statistics and prevalence distributions

**[smoking_overview.R](R/smoking_overview.R)**
- Temporal analysis of smoking periods
- Pack-year calculations and cessation timing

**[facial_pain_overview.R](R/facial_pain_overview.R)**
- Pain prevalence and characteristic distributions

**[chronic_diseases_overview.R](R/chronic_diseases_overview.R)**
- Chronic disease prevalence and patterns

**[covid_overview.R](R/covid_overview.R)**
- COVID-19 exposure patterns and temporal characteristics

**[nasal_breathing_and_ENT_surgery_overview.R](R/nasal_breathing_and_ENT_surgery_overview.R)**
- Nasal and ENT procedure analysis

---

## Study Cohort Characteristics

---

## Usage

**Data Preparation** (required first):
```r
source("R/read_data_and_basic_corrections.R")
source("R/build_analysis_dataset.R")
```

**Exploratory Analyses** (can be run independently):
```r
source("R/descriptive_statistics_overview.R")
source("R/smoking_overview.R")
source("R/facial_pain_overview.R")
source("R/chronic_diseases_overview.R")
source("R/covid_overview.R")
source("R/nasal_breathing_and_ENT_surgery_overview.R")
```

**Advanced Statistical Analysis**:
```r
source("R/trigeminal_measures_distributions_correlations.R")
source("R/ammo_distribution_analysis.R")
```

---

- **Sample size**: 1,001 participants
- **Gender**: 598 women, 396 men, 5 diverse (2 unreported)
- **Age**: Mean 44.2 years (SD 16.8, range 18-86)
- **Weight**: Mean 76.6 kg (SD 17.1, range 43-172 kg)
- **Height**: Mean 172.2 cm (SD 9.6, range 146-208 cm)

**Health Characteristics:**
- Chronic diseases: 431 participants (43.0%)
- ENT surgeries: 242 participants (24.2%)
- Facial pain: 52 participants (5.2%)
- Nasal breathing problems: 154 participants (15.4%)
- COVID-19 exposure: 834 participants (83.3%)

**Behavioral Characteristics:**
- Current smokers: 172 participants (17.2%)
- Former smokers: 239 participants (23.9%)
- Never smokers: 590 participants (59.0%)
- Alcohol consumption: 822 participants (82.2%)

**Trigeminal Sensitivity Measures:**
- AmmoLa intensity: All 1,001 participants
- Lateralization test: 458 participants
- CO₂ threshold (breath-hold protocol): 336 participants
- TriFunQ questionnaire: 1,001 participants (0.47% missing)

---

## License

This project is part of ongoing research into trigeminal chemosensitivity and related medical factors. The analysis code and documentation are available at: https://github.com/JornLotsch/TrigeminalChemosensitivity
