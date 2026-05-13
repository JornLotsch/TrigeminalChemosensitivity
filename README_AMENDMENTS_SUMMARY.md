# README.md Amendments Summary
**Date:** 13. Mai 2026  
**Purpose:** Align README with detailed analysis design documented in Methods section and improve research transparency

---

## Key Changes

### 1. **Analysis Design Section — Complete Rewrite**

**Before:**
- Oversimplified two-stage workflow description
- Unclear distinction between exploratory vs. predictive phases
- No mention of dataset scope for each phase

**After:**
- **Clear distinction of two complementary levels:**
  - Internal trigeminal data structure (21 variables)
  - External modulatory factors (220 non-trigeminal predictors)

- **Part 1: Exploratory Analysis** — explicitly labeled as hypothesis-generating
  - Complete imputed dataset (n=1,001)
  - No held-out validation
  
- **Part 2: Machine Learning Framework** — explicitly labeled as confirmatory/predictive
  - Train/validation split (80%/20%)
  - All development on training set only
  - Validation on held-out data
  
- **Critical design feature highlighted:** Train/validation split executed BEFORE imputation

---

### 2. **Environment Section — Corrected & Expanded**

**Before:**
```
- Imputation: `missForest`
```

**After:**
```
- Imputation: `missRanger` (random forest-based with predictive mean matching)
- Dimensionality reduction: `FactoMineR`, `factoextra`, `EDOtrans`, `fastICA`, `Rtsne`, `uwot`
- Statistical analysis: [...], `car`, `glmnet`, `randomForest`, `Boruta`, `nnet`
```

**Rationale:** Corrected to match actual implementation; added missing packages for reproducibility

---

### 3. **Repository Structure — Part 1/Part 2 Labeling**

**Reorganized to explicitly show:**
- PART 1: EXPLORATORY ANALYSIS (n=1,001)
- PART 2: MACHINE LEARNING FRAMEWORK (n=800/201 split)

**Each file now includes:**
- Specific analytical methods (e.g., "6 projections × 8 algorithms")
- Projection methods explicitly listed: EDOtrans, PCA, ICA, MDS, t-SNE, UMAP
- Clustering algorithms: k-means, PAM, hierarchical (6 linkages)
- Quality metrics: Calinski-Harabasz, silhouette width
- Validation approach: nearest-centroid, balanced accuracy with bootstrap CI

---

### 4. **Data Processing Section — Comprehensive New Content**

**Added:**

#### Sample Inclusion & Protocol
- CO₂ data availability: 453/1,001 participants
- Breath-hold protocol improvements: 17.1% → 5.4% ceiling effects
- **Explicit statement:** Data from uncontrolled protocol excluded

#### Variable Recoding & Construction
- Specific categories and one-hot encoding details
- Final output: 241 variables, 1.75% missing

#### Distribution Assessment & Transformations
- **Tukey's ladder of powers** methodology (λ values)
- **Box-Cox and D'Agostino's K² testing**
- **Mathematical formulas for transformations:**
  - AmmoLa: Reflected sign-preserving log transformation
  - CO₂: Sign-inverted log transformation
  - Lateralization: Untransformed (approximately normal)

#### Dataset Split & Imputation
- **Critical detail:** Split BEFORE imputation (prevents information leakage)
- **missRanger configuration:**
  - Predictive mean matching (k=5)
  - 500 trees, min node size = 5
- **Imputation fidelity validation:**
  - Kernel-density overlap check
  - Anderson-Darling tests (p=1.0)

---

### 5. **Reported Analyses Section — Reorganized & Expanded**

**Before:**
- Vague descriptions; "48 combinations" not explained
- No methodological specificity
- Missing quality metrics

**After:**

#### Data Preparation
- Explicit mention of CO₂ sample filtering
- Train/validation split timing (BEFORE imputation)
- Expected output: 21 trigeminal variables

#### Part 1: Exploratory Analysis
- Clear that it uses complete merged dataset (n=1,001)
- Methods: Spearman correlations, PCA, ABC analysis, exploratory regression
- Explicit note: These are hypothesis-generating; interpretation requires Methods context

#### Part 2: Machine Learning Framework
- **Each analysis step fully detailed:**
  
  **Clustering (Step 2a):**
  - 6 specific projection methods named
  - 8 specific clustering algorithms listed
  - NbClust with 30 validity indices
  - Calinski-Harabasz as primary criterion
  - Null model comparison

  **Supervised Classification (Step 2b):**
  - Outcome: cluster membership (categorical)
  - Predictors: 20 trigeminal variables (specified why 20, not 21)
  - Models: multinomial logistic, ridge, lasso, elastic net, RF
  - Hyperparameter tuning: 5-fold CV (penalized), grid search (RF)
  - Validation: balanced accuracy with bootstrap CI

  **Cluster Interpretation (Step 2c):**
  - Kruskal-Wallis with effect sizes (η²)
  - Psychometric dimensions: sensitivity + avoidance (a priori classification)
  - Cross-dataset concordance: Kendall's τ

  **Psychophysical Regression (Steps 2d-e):**
  - Distinct analyses for trigeminal vs. external predictors
  - Feature importance: Lasso + Boruta algorithms
  - Cross-dataset validation: Wilcoxon/Spearman + Kendall's τ concordance

---

### 6. **Usage Section — Restructured for Clarity**

**Before:**
- Mixed "Steps 1a, 1b, 1c" with "Steps 2a, 2b"
- Unclear which scripts to run in what order
- Optional descriptive scripts mixed with core pipeline

**After:**
- **Clear linear workflow:**
  1. Data Preparation (required)
  2. Part 1: Exploratory Analysis (complete dataset)
  3. Part 2: Machine Learning (train/validation with substeps a-d)
  4. Optional: Exploratory Population Context

- **Each step includes:**
  - Specific files to source
  - Expected outputs
  - Separation of required vs. optional

---

## Scientific Transparency Improvements

| Issue | Before | After |
|-------|--------|-------|
| Imputation method | Incorrect (missForest) | Corrected (missRanger) |
| Imputation parameters | Not mentioned | Fully specified (k=5, 500 trees, min node size=5) |
| CO₂ sample filtering | Not mentioned | Explicitly documented with protocol details |
| Clustering methods | "48 combinations" | 6 methods + 8 algorithms explicitly listed |
| Quality metrics | Not mentioned | Calinski-Harabasz, silhouette, null comparison |
| PCA variable importance | Not mentioned | ABC analysis method documented |
| Distribution testing | Not mentioned | Tukey's ladder + D'Agostino's K² detailed |
| Transformation formulas | Not mentioned | Mathematical formulas provided |
| Part 1 vs. Part 2 distinction | Vague | Crystal clear with n values |
| Multiple testing note | Missing | Flagged with reference to Methods context |
| Cross-dataset validation | Mentioned but vague | Kendall's τ concordance explicitly specified |

---

## Remaining Considerations

### Suggested Future Enhancements

1. **Add "Multiple Testing Framework" subsection** — Even brief note that Part 1 exploratory screening at α=0.05 expected ~11 false positives among 220 candidate variables; refer to manuscript Methods for correction approach (Bonferroni, FDR, etc.)

2. **Consider adding "Model Performance Criteria"** — Explicitly state generalization success criteria (e.g., "Successful generalization required significant model comparison and positive significant slope" for continuous outcomes; balanced accuracy thresholds for classification)

3. **Add note on external code repos** — Verify that links to external repos (EDOtrans, explore_tukey_lop.py) are current and accessible

4. **Consider "Limitations & Design Choices"** section — Document why certain design decisions were made (e.g., 80/20 split, why 5-fold CV, why specific validation metrics)

### Scientific Rigor Assessment

**Strengths:**
✓ Proper train/validation split BEFORE imputation (excellent practice)
✓ Systematic clustering comparison across 48 solutions
✓ Clear separation of exploratory vs. predictive phases
✓ Cross-dataset concordance validation
✓ Multiple regression frameworks (OLS, penalized, RF)
✓ Rigorous hyperparameter tuning

**Addressed Concerns:**
✓ Imputation method now correctly documented
✓ High missing data (CO₂, lateralization) acknowledged with validation checks
✓ Exploratory multiple testing now flagged

**Remaining Attention Points:**
⚠ Ensure Methods section specifies multiple testing correction for Part 1 screening
⚠ Confirm imputation assumptions appropriate for 66.5% missing CO₂ data
⚠ Verify external links remain accessible

