# Trigeminal Chemosensitivity Analysis Companion

This repository contains the custom code used to perform the computational analyses reported in the accompanying manuscript on trigeminal chemosensitivity phenotypes and their candidate modulators.

It is intended to support transparency and review of the analytical workflow, and it accompanies the published study as a reproducible computational companion.

The custom code used to conduct the analyses reported in this paper is openly available at https://github.com/JornLotsch/TrigeminalChemosensitivity.

---

## Scope and purpose

The repository documents the computational workflow underlying the study in a transparent and reviewable form. It contains the scripts used for data preparation, exploratory analyses, phenotype derivation, and predictive modeling, together with generated outputs and summary tables.

The manuscript remains the source of the scientific interpretation and methodological detail.

---

## Data availability

The data supporting the findings reported in this study are available from the senior author upon reasonable request, subject to approval by the responsible ethics committee.

Accordingly, the underlying raw and processed study data are not distributed in this repository. The repository provides the analysis code and generated outputs rather than replacing access to the source data.

---

## Analysis structure

The analysis workflow is organized into three stages:

1. **Observed-data characterization**
   - Characterization of the three psychophysical trigeminal measures in the non-imputed dataset.
   - Assessment of associations among measures and their overlap at the high-sensitivity extremes.
   - Comparison of sex differences in observed data.

2. **Exploratory imputed-data analysis**
   - Characterization of the complete 21-variable trigeminal dataset after imputation.
   - Correlation analysis, principal component analysis, and variable-importance ranking.
   - Screening of candidate non-trigeminal modulators for associations with the psychophysical outcomes.
   - These analyses were explicitly exploratory and hypothesis-generating.

3. **Predictive modeling and validation**
   - Partitioning of the dataset into training and validation sets before imputation and before model development.
   - Independent imputation in the training and validation subsets.
   - Derivation of a composite trigeminal phenotype using unsupervised clustering and dimensionality-reduction approaches.
   - Interpretation of cluster structure using factor analysis and effect-size comparisons.
   - Training of supervised models to predict psychophysical outcomes and cluster membership.
   - Evaluation of generalization in the independent validation set only.

---

## Repository structure

```text
.
├── README.md
├── R/
│   ├── 01_stage1_2_trigeminal_data_space_exploration.R
│   ├── 02_stage3_trigeminal_clustering_phenotype_derivation.R
│   ├── 03_stage3_trigeminal_cluster_interpretation.R
│   ├── 04_stage3_psychophysical_regression_analysis.R
│   ├── 05_stage3_modulator_identification_psychophysical.R
│   ├── 06_stage3_modulator_identification_clusters.R
│   ├── 07_stage3_trigeminal_cluster_modulators_interpretation.R
│   ├── build_analysis_dataset.R
│   ├── read_data_and_basic_corrections.R
│   ├── globals.R
│   ├── ProjectionsBiomed_MainFunctions_6_1core.R
│   ├── utils.R
│   ├── descriptive_statistics_overview.R
│   ├── smoking_overview.R
│   ├── facial_pain_overview.R
│   ├── chronic_diseases_overview.R
│   ├── covid_overview.R
│   ├── nasal_breathing_and_ENT_surgery_overview.R
│   └── ...
├── Python/
│   ├── trigeminal_measures_PCA.py
│   ├── compute_PCA_altColors.py
│   ├── cABC_analysis.py
│   ├── trigeminal_measures_exploreTukey.py
│   ├── explore_tukey_lop.py
│   └── ...
├── results/
│   ├── summary tables
│   ├── generated figures
│   └── output CSV files
├── Matlab/
├── OLDstuff/
└── .gitignore
```

The repository contains both R and Python analyses, generated outputs, and supporting scripts. Some intermediate data products and local analysis objects are stored outside of version control and are therefore not included here.

---

## Computational workflow

The scripts are organized according to the staged analysis workflow described in the manuscript. In general, the order of processing is:

1. Data preparation and harmonization
2. Observed-data characterization
3. Exploratory trigeminal data-space analysis
4. Candidate-modulator screening
5. Training/validation split and independent imputation
6. Phenotype derivation and validation
7. Predictive modeling for psychophysical outcomes and cluster membership

This repository provides a transparent overview of the analysis pipeline; the full methodological description and interpretation remain in the manuscript.

---

## Environment

The project was run in a Linux environment using:

- R 4.5.x / 4.6.x
- Python 3.12.x

Main packages used in the core analysis workflow include `dplyr`, `tidyr`, `ggplot2`, `missRanger`, `FactoMineR`, `factoextra`, `cluster`, `NbClust`, `psych`, `glmnet`, `randomForest`, `Boruta`, `nnet`, `boot`, `car`, `ComplexHeatmap`, `uwot`, `Rtsne`, and `fastICA`. Additional helper packages are used in supporting scripts and exploratory analyses.

---

## Usage notes

This repository is intended for inspection of the computational pipeline and for transparent review of the analysis workflow. It is best used alongside the manuscript, which provides the scientific interpretation and detailed methodological context.

The code is organized to reflect the staged analysis workflow of the study.

---

## Citation

Lötsch J, Weise S, Himmelspach A, Bormann A, Hummel T.
Mapping population diversity in intranasal trigeminal chemosensitivity reveals multidimensional phenotypes.
Chemical Senses, 2026 (in revision).

---

## License

Code and documentation are shared under the CC BY 4.0 license.
