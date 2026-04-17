################################################################################
# Trigeminal Sensitivity Analysis - Clustering and Validation
#
# Author: Joern Lotsch
# Date: 2025-01-26
#
# Description:
#   Unsupervised clustering analysis of trigeminal sensitivity measures using
#   multiple projection methods (PCA, ICA, MDS, t-SNE, UMAP) and clustering
#   algorithms (k-means, hierarchical clustering). Includes validation on an
#   independent dataset and comparison of cluster-wise variable distributions.
#
# Input Files:
#   - analysis_dataset_training_imputed.csv: Training data with imputed values
#   - analysis_dataset_validation_imputed.csv: Validation data with imputed values
#   - trigeminale_daten_corrected_translated.csv: Raw complete dataset
#   - globals.R: Global variables and color palettes
#   - ProjectionsBiomed_MainFunctions_6_1core.R: Projection/clustering functions
#
# Output Files:
#   - Combined_projection_and_clustering_analysis_plot_*.svg: All clustering results
#   - dfClusterQuality.csv: Quality metrics for all clustering solutions
#   - p_GMM_row_means.svg/png: Distribution comparison (training vs validation)
#   - combined_trig_clustering_plot.svg/png: Best clustering diagnostics
#   - p_combined_clusters_train/valid.svg/png: Variable analysis by cluster
#   - heat_plot_clustered_variables.svg/png: Heatmap of clustered variables
#
# Key Analyses:
#   1. Data preprocessing and transformation
#   2. Distribution analysis (Pareto Density Estimation, Gaussian Mixture Models)
#   3. Projection-based clustering (6 methods × 8 algorithms)
#   4. Cluster quality evaluation and selection
#   5. Validation on independent dataset
#   6. Variable-wise cluster comparisons (effect sizes: rank-biserial or eta²)
#   7. Cross-dataset correlation of effect sizes
#
# Statistical Methods:
#   - Non-parametric tests: Wilcoxon (2 groups), Kruskal-Wallis (>2 groups)
#   - Effect sizes: Rank-biserial correlation (2 groups), Eta² (>2 groups)
#   - Bootstrap confidence intervals (1000 iterations)
#   - Kendall's tau for cross-dataset correlation
################################################################################

# ============================================================================ #
# 1. LOAD REQUIRED LIBRARIES
# ============================================================================ #

# Core data manipulation
library(dplyr)          # Data manipulation
library(tidyr)          # Data tidying
library(tibble)         # Modern data frames
library(purrr)          # Functional programming
library(reshape2)       # Data reshaping

# Visualization
library(ggplot2)        # Core plotting
library(ggpubr)         # Publication-ready plots
library(ggplotify)      # Convert base plots to ggplot
library(ggpmisc)        # Additional ggplot2 functionality
library(ggthemes)       # Extra themes
library(cowplot)        # Plot composition
library(patchwork)        # Plot composition
library(grid)           # Grid graphics
library(gridExtra)      # Multiple grid plots
library(ComplexHeatmap) # Advanced heatmaps
library(circlize)       # Color mapping for heatmaps
library(viridis)        # Color scales
library(colorspace)     # Color manipulation
library(scales)         # Scale functions

# Statistical analysis
library(boot)           # Bootstrap methods
library(effsize)        # Effect size calculations
library(Hmisc)          # Statistical tools
library(psych)          # Correlation analysis
library(MASS)           # Robust statistics

# Clustering and dimensionality reduction
library(cluster)        # Clustering algorithms
library(NbClust)        # Optimal cluster number
library(factoextra)     # Factor visualization
library(FactoMineR)     # Factor analysis

# Specialized analyses
library(DataVisualizations) # Pareto Density Estimation
library(opGMMassessment)    # Gaussian Mixture Models
library(missForest)         # Imputation methods
library(cvms)               # Confusion matrix visualization
library(vcd)                # Categorical data visualization

# Utilities
library(forcats)        # Factor handling
library(stringr)        # String manipulation
library(lubridate)      # Date/time handling
library(pbmcapply)      # Parallel apply functions

# Load custom functions and global variables
source("globals.R")
source("ProjectionsBiomed_MainFunctions_6_1core.R")
source("/home/joern/Aktuell/ABCplotGG.R")

# ============================================================================ #
# 2. GLOBAL OPTIONS AND ANALYSIS PARAMETERS
# ============================================================================ #

# Analysis control switches
remove_censored <- FALSE            # Whether to exclude censored values
scale_0_100 <- FALSE                # Whether to scale all measures to 0-100
analyze_only_untransformed <- FALSE # Skip transformation analysis
plot_only_untransformed <- TRUE     # Show only original data in plots

# Computational resources
nProc_possible <- parallel::detectCores() - 1
nProc_desired <- 48

# Dimensionality reduction methods to evaluate
projection_methods <- c("none", "PCA", "ICA", "MDS", "tSNE", "Umap")

# Clustering algorithms to evaluate
clustering_methods <- c("kmeans", "kmedoids",
                        "ward.D2", "single", "average",
                        "median", "complete", "centroid")

# Method for determining optimal number of clusters
cluster_number_methods <- "NbClust"

# Clustering quality metrics (full set)
unsupervised_metrics_full <- c("Silhouette_index", "Dunn_index",
                               "DaviesBouldin_index", "dbcv_index",
                               "CalinskiHarabasz_index", "inertia")

# Clustering quality metrics (reduced set)
unsupervised_metrics_reduced <- c("Silhouette_index", "Dunn_index",
                                  "CalinskiHarabasz_index")

# Validation and comparison metrics
valid_cluster_metrics <- c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                           "Rand_index", "DaviesBouldin_index", "dbcv_index",
                           "CalinskiHarabasz_index", "inertia",
                           "adjusted_mutual_information")

# Color palettes
original_colors <- c("cornsilk1", "cornsilk2", "cornsilk3", "cornsilk4")
dark_colors <- darken(original_colors, amount = 0.3)

# ============================================================================ #
# 3. HELPER FUNCTIONS
# ============================================================================ #

#' Sign-Preserving Logarithmic Transformation (Zero-Invariant)
#'
#' Applies logarithmic transformation while preserving the sign of the input.
#' Handles zero values gracefully using log1p.
#'
#' @param x Numeric vector to transform
#' @param base Logarithm base: 0 (natural log), 2, 10, or custom value
#' @return Transformed vector maintaining sign of original values
#' @examples
#' slog(c(-10, 0, 10), base = 10)
slog <- function(x, base = 10) {
  absX <- abs(x)
  s <- sign(x)

  if (base == 0) {
    return(s * log1p(absX))
  } else if (base == 2) {
    return(s * log2(absX + 1))
  } else if (base == 10) {
    return(s * log10(absX + 1))
  } else {
    return(s * log1p(absX) * log(base))
  }
}

#' Inverse Sign-Preserving Logarithmic Transformation
#'
#' Reverses the slog transformation to recover original values.
#'
#' @param y Transformed values from slog()
#' @param base Logarithm base used in forward transformation
#' @return Original scale values
inv_slog <- function(y, base = 10) {
  s <- sign(y)
  absY <- abs(y)

  if (base == 0) {
    val <- expm1(absY)
  } else if (base == 2) {
    val <- 2^absY - 1
  } else if (base == 10) {
    val <- 10^absY - 1
  } else {
    val <- expm1(absY / log(base))
  }
  return(s * val)
}

#' Reflected Logarithmic Transformation
#'
#' Applies slog transformation after reflecting data around maximum value.
#' Useful for right-skewed distributions.
#'
#' @param x Numeric vector to transform
#' @return Reflected and log-transformed vector
reflect_slog <- function(x) {
  slog(max(x, na.rm = TRUE) + 1 - x)
}

#' Reflected Logarithmic Transformation (Unflipped)
#'
#' Same as reflect_slog but maintains original direction by negating result.
#'
#' @param x Numeric vector to transform
#' @return Reflected, log-transformed, and negated vector
reflect_slog_unflipped <- function(x) {
  -slog(max(x, na.rm = TRUE) + 1 - x)
}

#' Inverse Reflected Logarithmic Transformation (Unflipped)
#'
#' Reverses the reflect_slog_unflipped transformation.
#'
#' @param y Transformed values
#' @param original_max Maximum value from original data
#' @param base Logarithm base used in transformation
#' @return Original scale values
inv_reflect_slog_unflipped <- function(y, original_max, base = 10) {
  M <- original_max
  x_original <- M + 1 - inv_slog(-y, base)
  return(x_original)
}

#' Create Pareto Density Estimation Plot with ggplot2
#'
#' Generates density plots using Pareto Density Estimation for each column
#' in the input data matrix.
#'
#' @param Data Numeric matrix or data frame
#' @return ggplot object with PDE curves
PDEplotGG <- function(Data) {
  Data <- as.matrix(Data)
  m <- matrix(NA, nrow = 0, ncol = 3)

  require(DataVisualizations)

  # Calculate PDE for each column
  for (i in seq_len(ncol(Data))) {
    PDE <- ParetoDensityEstimation(as.vector(na.omit(Data[, i])))
    m2 <- PDE$kernels
    m3 <- PDE$paretoDensity
    m1 <- rep(i, length(m2))
    m <- rbind(m, cbind(m1, m2, m3))
  }

  mdf <- data.frame(m)
  require(ggplot2)

  p <- ggplot(data = mdf, aes(x = m2, y = m3, colour = factor(m1))) +
    geom_line(aes(linewidth = 1)) +
    guides(linewidth = FALSE)

  return(p)
}

# Helper function to select colors for dendrogram branches
select_extremes <- function(group, k_clusters) {
  n_items <- length(group)
  if (n_items <= k_clusters) {
    return(group)
  } else {
    positions <- seq(1, n_items, length.out = k_clusters)
    return(group[round(positions)])
  }
}


# ============================================================================ #
# 4. LOAD AND PREPARE DATA
# ============================================================================ #

cat("\n=== Loading Data ===\n")

# Load training data (imputed)
trigeminale_training_data <- read.csv("analysis_dataset_training_imputed.csv",
                                      check.names = FALSE)
cat("Training data loaded:", nrow(trigeminale_training_data), "samples\n")

# Extract variables for clustering
variables_for_clustering_imputed <- trigeminale_training_data[,
  c(variables_by_categories$Nasal_chemosensory_perception,
    variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
rownames(variables_for_clustering_imputed) <- trigeminale_training_data$ID

# Load validation data (imputed)
trigeminale_validation_data <- read.csv("analysis_dataset_validation_imputed.csv",
                                        check.names = FALSE)
cat("Validation data loaded:", nrow(trigeminale_validation_data), "samples\n")

# Extract variables for clustering (validation)
variables_for_clustering_imputed_validation <- trigeminale_validation_data[,
  c(variables_by_categories$Nasal_chemosensory_perception,
    variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
rownames(variables_for_clustering_imputed_validation) <- trigeminale_validation_data$ID

# Load raw data (for reference)
trigeminale_all_data_raw <- read.csv("trigeminale_daten_corrected_translated.csv",
                                     check.names = FALSE)

# ============================================================================ #
# 5. DATA TRANSFORMATION AND PREPROCESSING
# ============================================================================ #

cat("\n=== Transforming Variables ===\n")

# Backup original variable names
original_names <- names(variables_for_clustering_imputed)

# Transform psychophysical variables (training)
# AmmoLa intensity: reflected slog (handles ceiling effects)
variables_for_clustering_imputed$`AmmoLa intensity` <-
  scaleRange_01(reflect_slog_unflipped(variables_for_clustering_imputed$`AmmoLa intensity`)) * 3

# CO2 threshold: negative slog (lower = more sensitive)
variables_for_clustering_imputed$`CO2 threshold` <-
  scaleRange_01(-slog(variables_for_clustering_imputed$`CO2 threshold`)) * 3

# Lateralization: negative slog (lower = more sensitive)
variables_for_clustering_imputed$`Lateralization (x/20)` <-
  scaleRange_01(-slog(variables_for_clustering_imputed$`Lateralization (x/20)`)) * 3

# Transform psychophysical variables (validation)
variables_for_clustering_imputed_validation$`AmmoLa intensity` <-
  scaleRange_01(reflect_slog_unflipped(variables_for_clustering_imputed_validation$`AmmoLa intensity`)) * 3

variables_for_clustering_imputed_validation$`CO2 threshold` <-
  scaleRange_01(-slog(variables_for_clustering_imputed_validation$`CO2 threshold`)) * 3

variables_for_clustering_imputed_validation$`Lateralization (x/20)` <-
  scaleRange_01(-slog(variables_for_clustering_imputed_validation$`Lateralization (x/20)`)) * 3

cat("Variables transformed to [0, 3] scale\n")

# Assign generic variable names for internal processing
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed <-
  variables_for_clustering_imputed
names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed) <-
  paste0("Var", seq_len(ncol(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed)))

# Prepare dataset for analysis (adds Target and Label columns)
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared <-
  prepare_dataset(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed)

# Inspect dataset structure
cat("\nDataset structure:\n")
print(names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared))
cat("\nTarget distribution:\n")
print(table(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Target))

# ============================================================================ #
# 6. DISTRIBUTION ANALYSIS
# ============================================================================ #

cat("\n=== Analyzing Distributions ===\n")

# Prepare matrices for visualization (rename transformed variables)
heat_matrix <- variables_for_clustering_imputed %>%
  rename(AmmoLa_intensity_transformed = `AmmoLa intensity`,
         CO2_threshold_transformed = `CO2 threshold`)

heat_matrix_validation <- variables_for_clustering_imputed_validation %>%
  rename(AmmoLa_intensity_transformed = `AmmoLa intensity`,
         CO2_threshold_transformed = `CO2 threshold`)

# Calculate row means (average trigeminal sensitivity per subject)
row_means <- rowMeans(heat_matrix)
row_means_validation <- rowMeans(heat_matrix_validation)

# Pareto Density Estimation plots
p_PDE_chemosensory_perception_for_clustering_imputed <- PDEplotGG(row_means)
print(p_PDE_chemosensory_perception_for_clustering_imputed)

p_PDE_chemosensory_perception_for_clustering_imputed_validation <- PDEplotGG(row_means_validation)
print(p_PDE_chemosensory_perception_for_clustering_imputed_validation)

# Gaussian Mixture Model fitting (training)
cat("Fitting Gaussian Mixture Models (training)...\n")
GMM_chemosensory_perception_for_clustering_imputed <-
  opGMMassessment::opGMMassessment(row_means, MaxModes = 4, FitAlg = "DO",
                                   MaxCores = 4, PlotIt = TRUE, Seed = 42)
cat("GMM means (training):", GMM_chemosensory_perception_for_clustering_imputed$Means, "\n")
cat("GMM SDs (training):", GMM_chemosensory_perception_for_clustering_imputed$SDs, "\n")

# Gaussian Mixture Model fitting (validation)
cat("Fitting Gaussian Mixture Models (validation)...\n")
GMM_chemosensory_perception_for_clustering_imputed_validation <-
  opGMMassessment::opGMMassessment(row_means_validation, MaxModes = 4, FitAlg = "DO",
                                   MaxCores = 4, PlotIt = TRUE, Seed = 42)
cat("GMM means (validation):", GMM_chemosensory_perception_for_clustering_imputed_validation$Means, "\n")
cat("GMM SDs (validation):", GMM_chemosensory_perception_for_clustering_imputed_validation$SDs, "\n")

# Normality tests
cat("\nAnderson-Darling test for normality (training):\n")
print(DescTools::AndersonDarlingTest(row_means, "pnorm",
                                     mean = mean(row_means), sd = sd(row_means)))

cat("\nAnderson-Darling test for normality (validation):\n")
print(DescTools::AndersonDarlingTest(row_means_validation, "pnorm",
                                     mean = mean(row_means_validation),
                                     sd = sd(row_means_validation)))

# Compare training vs validation distributions
cat("\nComparing training vs validation distributions...\n")
ad_result <- twosamples::ad_test(row_means, row_means_validation)
ad_p <- ad_result["P-Value"]
cat("Anderson-Darling test p-value:", format(ad_p, scientific = FALSE, digits = 3), "\n")

# Visualization of distribution comparison
df_row_means <- rbind.data.frame(
  cbind.data.frame(Data = "Training", values = row_means),
  cbind.data.frame(Data = "Validation", values = row_means_validation)
)

p_GMM_row_means <- ggplot(df_row_means, aes(x = values, color = Data, fill = Data)) +
  geom_density(alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  scale_color_manual(values = c("cornsilk2", "cornsilk4")) +
  labs(title = "Row means of trigeminal measures matrix",
       subtitle = "Training versus validation",
       x = "Means", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.2, .8),
        legend.direction = "vertical") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Anderson-Darling test\np = ",
                         format(ad_p, scientific = FALSE, digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           color = "#444444", family = "Libre Franklin",
           fontface = "plain")

print(p_GMM_row_means)

# Save distribution comparison
ggsave("p_GMM_row_means.svg", p_GMM_row_means, width = 8, height = 8, dpi = 300)
ggsave("p_GMM_row_means.png", p_GMM_row_means, width = 8, height = 8, dpi = 300)

# ============================================================================ #
# 7. PROJECTION AND CLUSTERING ANALYSIS
# ============================================================================ #

cat("\n=== Performing Projection and Clustering Analysis ===\n")
cat("This may take several minutes...\n")

# Set dataset names for analysis
DatasetNames <- c("variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared")

# Set seed for reproducibility
set.seed(42)

# Run comprehensive projection and clustering analysis
# Tests all combinations of projection methods and clustering algorithms
projectionsAndPlots_TriFunQ <- perform_analysis(
  datasets = DatasetNames,
  projection_methods = projection_methods,
  clustering_methods = clustering_methods,
  cluster_number_methods = cluster_number_methods,
  label_points = FALSE,
  highlight_misclassified = FALSE,
  selected_cluster_metrics = unsupervised_metrics_full,
  highlight_best_clustering = TRUE,
  palette_target = cb_palette,
  palette_cluster = cb_palette,
  seed = 42,
  cells_colored_for = "Cluster",
  points_colored_for = "Cluster",
  nProc = max(1, min(nProc_desired, nProc_possible)),
  max_clusters = 5
)

cat("Analysis complete!\n")

# Optionally save results for later reloading
# saveRDS(projectionsAndPlots_TriFunQ, file = "projectionsAndPlots_TriFunQ_training.RData")
# projectionsAndPlots_TriFunQ <- readRDS("projectionsAndPlots_TriFunQ_training.RData")

# ============================================================================ #
# 8. COMBINE AND SAVE PROJECTION PLOTS
# ============================================================================ #

cat("\n=== Creating Combined Projection Plots ===\n")

# Combine all projection plots into single figure

# remove plots with projection = "none" (makes no sense to plot variables 1 and 2)

no_none_projection_plots <- as.numeric(strsplit(paste(setdiff(1:length(names(projectionsAndPlots_TriFunQ$projections_plots$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared))
                                    , grep("none", names(projectionsAndPlots_TriFunQ$projections_plots$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared))
), collapse = ","), ",")[[1]])

combined_plots <- combine_all_plots(
  datasets = DatasetNames,
  projection_plots = projectionsAndPlots_TriFunQ$projections_plots$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[no_none_projection_plots],
  projection_methods = projection_methods[projection_methods != "none"],
  clustering_methods = clustering_methods,
  cluster_number_methods = cluster_number_methods
)

print(combined_plots)
# Apply title updates
original_plots_updated <- update_plot_titles(original_plots, title_size = 12, wrap_chars = 15)

# Recombine with updated titles
combined_plots_fixed <- combine_all_plots(
  datasets = DatasetNames,
  projection_plots = original_plots_updated$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[no_none_projection_plots],
  projection_methods = projection_methods[projection_methods != "none"],
  clustering_methods = clustering_methods,
  cluster_number_methods = cluster_number_methods
)


# Add title above (cowplot style)
combined_plots_fixed_title <- cowplot::plot_grid(
  ggplot() +
    labs(title = "Comprehensive combinations of projection and clustering methods.",
         subtitle = "Rows: Projection Methods, Columns: Clustering Methods;\nLabelling: Projection method - clustering method - cluster number detection method (always NbClust)") +
    theme_void() +
    theme(plot.title = element_text(size = 25, face = "plain", margin = margin(b = 4)),
          plot.subtitle = element_text(size = 20, face = "plain", margin = margin(b = 1))),
  combined_plots_fixed$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared,
  ncol = 1,
  rel_heights = c(0.05, 1),
  align = "v"
)

print(combined_plots_fixed_title)

# Save combined plot
ggsave(
  filename = paste0("Combined_projection_and_clustering_analysis_plot_",
                    cluster_number_methods, "_clusters_", DatasetNames[1], ".svg"),
  plot = combined_plots_fixed_title,
  width = 4 * (length(clustering_methods) + length(cluster_number_methods)),
  height = 4.5 * length(projection_methods[projection_methods != "none"]),
  limitsize = FALSE
)
ggsave(
  filename = paste0("Combined_projection_and_clustering_analysis_plot_",
                    cluster_number_methods, "_clusters_", DatasetNames[1], ".png"),
  plot = combined_plots_fixed_title,
  width = 4 * (length(clustering_methods) + length(cluster_number_methods)),
  height = 4.5 * length(projection_methods[projection_methods != "none"]),
  limitsize = FALSE,
  dpi = 300
)

cat("Combined projection plot saved\n")

# ============================================================================ #
# 9. CLUSTER QUALITY EVALUATION
# ============================================================================ #

cat("\n=== Evaluating Cluster Quality ===\n")

# Extract cluster quality results
dfClusterQuality <- projectionsAndPlots_TriFunQ$cluster_quality_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed

# Remove metrics not used in current analysis
excluded_metrics <- setdiff(valid_cluster_metrics, "CalinskiHarabasz_index")
excluded_metrics <- setdiff(valid_cluster_metrics, unsupervised_metrics_full)

dfClusterQuality <- dfClusterQuality[, !names(dfClusterQuality) %in%
                                       c(excluded_metrics,
                                         paste0(excluded_metrics, "_rank"),
                                         "combined_rank_metrics_with_orig_classes",
                                         "combined_rank")]

# Calculate combined rank score (product of individual ranks)
dfClusterQuality <- dfClusterQuality %>%
  rownames_to_column(var = "row_id") %>%
  rowwise() %>%
  mutate(
    combined_rank_metrics_without_orig_classes = {
      rank_cols <- grep("_rank$", names(cur_data()), value = TRUE)
      if (length(rank_cols) == 0) {
        return(NA_real_)
      }
      prod(as.numeric(cur_data()[1, rank_cols, drop = FALSE]), na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  column_to_rownames(var = "row_id")

# Generate unique method identifiers
rownames(dfClusterQuality) <- apply(dfClusterQuality[, 1:3], 1, paste0, collapse = "_")
dfClusterQuality$Method <- row.names(dfClusterQuality)

# Sort by best overall performance (highest rank product)
dfClusterQuality <- dfClusterQuality[order(-dfClusterQuality$combined_rank_metrics_without_orig_classes), ]

cat("Best clustering method:", rownames(dfClusterQuality)[1], "\n")
cat("Combined rank score:", dfClusterQuality$combined_rank_metrics_without_orig_classes[1], "\n")

# Visualize cluster quality rankings as heatmap
heatmap_data <- melt(toPercent(dfClusterQuality[, grep("_rank", names(dfClusterQuality))]))
colnames(heatmap_data) <- c("row", "column", "value")
heatmap_data$row <- forcats::fct_rev(as.factor(heatmap_data$row))

clustering_heatmap <- ggplot(heatmap_data, aes(x = column, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("cornsilk", "cornsilk4")) +
  theme_plot() +
  theme(
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6)
  ) +
  labs(title = "Ranked cluster quality scores", fill = "Scaled\nrank",
       y = "Projection_clustering_clusternumber", x = NULL)

print(clustering_heatmap)

# Bar plot of combined rank scores
ggplot(dfClusterQuality, aes(y = combined_rank_metrics_without_orig_classes,
                              x = reorder(Method, -combined_rank_metrics_without_orig_classes))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

# Save cluster quality results
dfClusterQuality[] <- lapply(dfClusterQuality, function(x) {
  if (is.numeric(x)) round(x, 3) else x
})
write.csv(dfClusterQuality, "dfClusterQuality.csv", row.names = TRUE)
cat("Cluster quality metrics saved to dfClusterQuality.csv\n")

# ============================================================================ #
# 10. EXTRACT BEST CLUSTERING SOLUTION
# ============================================================================ #

cat("\n=== Extracting Best Clustering Solution ===\n")

# Extract best methods from quality evaluation
best_projection_method <- dfClusterQuality$projection_method[1]
best_clustering_method <- dfClusterQuality$clustering_method[1]
best_cluster_number_method <- dfClusterQuality$cluster_number_method[1]

cat("Best projection method:", best_projection_method, "\n")
cat("Best clustering method:", best_clustering_method, "\n")
cat("Cluster number method:", best_cluster_number_method, "\n")

# Extract cluster assignments
TriFunQ_clusters <-
  projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[
    best_projection_method]][[best_clustering_method]][[best_cluster_number_method]]$Cluster

k_clusters <- length(unique(TriFunQ_clusters))
cat("Number of clusters:", k_clusters, "\n")
cat("Cluster sizes:\n")
print(table(TriFunQ_clusters))

# Extract instance labels
analysis_training_instances <- projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[
  best_projection_method]][[best_clustering_method]][[best_cluster_number_method]]$Label

cat("Number of instances:", length(TriFunQ_clusters), "\n")

# ============================================================================ #
# 11. AGREEMENT WITH AmmoLA GROUPING (FISHER'S EXACT TEST)
# ============================================================================ #

cat("\n=== Testing Agreement with AmmoLA Grouping ===\n")

# Prepare data for agreement analysis
trigeminal_measures_data_grouped_breath_hold_ID <- trigeminal_measures_data_grouped_breath_hold
trigeminal_measures_data_grouped_breath_hold_ID$ID <- trigeminale_data$ID

# Create AmmoLa intensity grouping (match to training instances)
AmmoLa_intensity <- trigeminal_measures_data_grouped_breath_hold_ID$AmmoLa_intensity[
  trigeminal_measures_data_grouped_breath_hold_ID$ID %in% analysis_training_instances
] + 1

# Create contingency table
table_AmmoLa_vs_TriFunQ_clusters <- table(AmmoLa_intensity, TriFunQ_clusters)
cat("\nContingency table:\n")
print(table_AmmoLa_vs_TriFunQ_clusters)

# Perform Fisher's exact test
ftest_AmmoLa_vs_TriFunQ_clusters <- fisher.test(table_AmmoLa_vs_TriFunQ_clusters)
cat("\nFisher's exact test:\n")
print(ftest_AmmoLa_vs_TriFunQ_clusters)

# Format test results for visualization
expected_counts <- chisq.test(table_AmmoLa_vs_TriFunQ_clusters)$expected
min_expected <- min(expected_counts)
stat_text_TriFunQ_clusters <- paste0(
  "Fisher's exact test\n",
  "p = ", signif(ftest_AmmoLa_vs_TriFunQ_clusters$p.value, 3), "\n",
  "Min expected = ", round(min_expected, 1)
)

# Create heatmap-style agreement table
p_table_AmmoLa_vs_TriFunQ_clusters <- ggplot(
  as.data.frame(table_AmmoLa_vs_TriFunQ_clusters),
  aes(x = TriFunQ_clusters, y = AmmoLa_intensity, fill = Freq)
) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "cornsilk1", high = "cornsilk4") +
  theme_minimal() +
  labs(
    title = "AmmoLa grouping vs TriFunQ clusters",
    x = "TriFunQ clusters",
    y = "AmmoLa_intensity",
    fill = "Count"
  ) +
  annotate(
    "rect", xmin = 1.1, xmax = 1.9, ymin = 1.2, ymax = 1.8,
    alpha = 0.5, fill = "white", color = NA
  ) +
  annotate(
    "text", x = 1.5, y = 1.5,
    label = stat_text_TriFunQ_clusters, size = 3, color = "black"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  coord_fixed(ratio = 1)

print(p_table_AmmoLa_vs_TriFunQ_clusters)

# ============================================================================ #
# 12. TRAINING CLUSTER DIAGNOSTICS
# ============================================================================ #

cat("\n=== Generating Training Cluster Diagnostics ===\n")

# A. Silhouette analysis
# Calculate dissimilarity matrix from projected data
diss_matrix_train <- dist(projectionsAndPlots_TriFunQ$projection_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[best_projection_method]]$Projected)

# Calculate silhouette widths
sil_train <- silhouette(projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[best_projection_method]][[best_clustering_method]][[best_cluster_number_method]]$Cluster, diss_matrix_train)

# Prepare silhouette data for plotting
sil_df_train <- as.data.frame(sil_train)
sil_df_train$case <- as.numeric(rownames(sil_df_train))
sil_df_train <- sil_df_train %>%
  arrange(cluster, -sil_width) %>%
  mutate(id = row_number())

# Create silhouette plot
p_TriFunQ_clusters_silhouette <- ggplot(sil_df_train, aes(x = id, y = sil_width, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
  labs(title = "Silhouette Plot - TriFunQ Training Clusters",
       x = "Case", y = "Silhouette Width", fill = "Cluster") +
  theme_plot() +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", strip.background = element_rect(fill = "cornsilk2")) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid")

print(p_TriFunQ_clusters_silhouette)

# B. Dendrogram (for hierarchical clustering methods only)
if (best_clustering_method %in% c("ward.D2", "single", "average", "median", "complete", "centroid")) {

  cat("Creating dendrogram for hierarchical clustering...\n")

  # Create dendrogram
  hc <- hclust(diss_matrix_train, method = best_clustering_method)
  require(magrittr)
  require(dendextend)

  dend <- hc %>%
    as.dendrogram %>%
    set("branches_k_color", value = select_extremes(c(dark_colors, "grey31"), k_clusters), k = k_clusters) %>%
    set("branches_lwd", 0.7) %>%
    set("labels_cex", 0.6) %>%
    set("labels_colors", value = select_extremes(c(dark_colors, "grey31"), k_clusters), k = k_clusters) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 0.5)

  ggd1 <- as.ggdend(dend)
  p_TriFunQ_clusters_dend <- ggplot(ggd1, horiz = FALSE) +
    labs(title = "Dendrogram - TriFunQ training", x = "Case", y = "Dissimilarity") +
    theme_plot() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  print(p_TriFunQ_clusters_dend)
}

# C. Projection scatter plot with clusters
df_projected_train_tri <- cbind.data.frame(
  Target = projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[best_projection_method]][[best_clustering_method]][[best_cluster_number_method]]$Cluster,
  projectionsAndPlots_TriFunQ$projection_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[[best_projection_method]]$Projected
)

p_TriFunQ_clusters_scatter <- ggplot(df_projected_train_tri,
                                     aes(x = Dim1, y = Dim2, shape = as.factor(Target),
                                         fill = as.factor(Target), color = as.factor(Target))) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(level = 0.95, aes(fill = as.factor(Target)), alpha = 0.3, geom = "polygon") +
  theme_plot() +
  labs(title = "Projection scatter plot with confidence ellipses",
       x = "Dim 1", y = "Dim 2", shape = "Cluster") +
  scale_color_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  theme(legend.position.inside = TRUE, legend.position = c(.1, .2)) +
  guides(color = "none", fill = "none")

print(p_TriFunQ_clusters_scatter)

# ============================================================================ #
# 13. VALIDATION CLUSTERING
# ============================================================================ #

cat("\n=== Applying Clustering to Validation Data ===\n")

# Prepare validation dataset
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_validation <-
  variables_for_clustering_imputed_validation
names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_validation) <-
  paste0("Var", seq_len(ncol(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_validation)))

# Extract training data (without Target and Label columns)
train_projection_data <- variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared[,
  !names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared) %in% c("Target", "Label")]

# Perform projection on training data
train_result_projection <- performProjection(X = train_projection_data,
                                             projection_method = best_projection_method)
train_model <- train_result_projection$projection_object

# Apply trained projection to validation data
valid_projection_data <- variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_validation
valid_result_projection <- performProjection(X = train_projection_data,
                                             projection_method = best_projection_method,
                                             projection_object = train_model,
                                             newdata = valid_projection_data)

# Apply clustering to validation data
if (best_clustering_method == "kmeans") {
  # K-means: fit to training, assign validation to nearest centers
  cat("Using k-means clustering...\n")
  set.seed(42)
  clusters_train <- kmeans(valid_result_projection$projection_object$embedding,
                          centers = length(unique(TriFunQ_clusters)), nstart = 100)

  # Function to assign new points to nearest cluster
  assign_to_nearest_cluster <- function(new_points, cluster_centers) {
    apply(new_points, 1, function(point) {
      which.min(apply(cluster_centers, 1, function(center) {
        sum((point - center)^2)
      }))
    })
  }

  clusters_val <- assign_to_nearest_cluster(valid_result_projection$Projected,
                                            clusters_train$centers)
  clusters_train <- clusters_train$cluster

} else {
  # Hierarchical clustering: fit to training, assign validation to nearest centers
  cat("Using hierarchical clustering...\n")
  set.seed(42)
  dist_matrix_train <- dist(train_result_projection$Projected, method = "euclidean")
  hclust_train <- hclust(dist_matrix_train, method = best_clustering_method)
  clusters_train <- cutree(hclust_train, k = length(unique(TriFunQ_clusters)))

  # Extract cluster centers
  cluster_centers_hclust <- tapply(train_result_projection$Projected,
                                   list(clusters_train),
                                   colMeans)
  cluster_centers_hclust <- do.call(rbind, cluster_centers_hclust)

  # Assign validation to nearest centers
  clusters_val <- assign_to_nearest_cluster(valid_result_projection$Projected,
                                           cluster_centers_hclust)
}

# Display cluster sizes
cat("Training projection dimensions:", dim(valid_result_projection$projection_object$embedding), "\n")
cat("Validation projection dimensions:", dim(valid_result_projection$Projected), "\n")
cat("Training cluster sizes:\n")
print(table(clusters_train))
cat("Validation cluster sizes:\n")
print(table(clusters_val))

# ============================================================================ #
# 14. VALIDATION CLUSTER DIAGNOSTICS
# ============================================================================ #

cat("\n=== Generating Validation Cluster Diagnostics ===\n")

# A. Silhouette analysis
diss_matrix_val <- dist(valid_result_projection$Projected)
sil_val <- silhouette(clusters_val, diss_matrix_val)
sil_df_val <- as.data.frame(sil_val)
sil_df_val$case <- as.numeric(rownames(sil_df_val))
sil_df_val <- sil_df_val %>%
  arrange(cluster, -sil_width) %>%
  mutate(id = row_number())

p_validation_clusters_silhouette <- ggplot(sil_df_val,
                                           aes(x = id, y = sil_width, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
  labs(title = "Silhouette Plot - Validation Clusters",
       x = "Case", y = "Silhouette Width", fill = "Cluster") +
  theme_plot() +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", strip.background = element_rect(fill = "cornsilk2")) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid")

print(p_validation_clusters_silhouette)

# B. Dendrogram (for hierarchical methods)
if (best_clustering_method %in% c("ward.D2", "single", "average", "median", "complete", "centroid")) {

  select_extremes <- function(group, k_clusters) {
    n_items <- length(group)
    if (n_items <= k_clusters) {
      return(group)
    } else {
      positions <- seq(1, n_items, length.out = k_clusters)
      return(group[round(positions)])
    }
  }

  hc <- hclust(diss_matrix_val, method = best_clustering_method)
  require(magrittr)
  require(dendextend)

  dend <- hc %>%
    as.dendrogram %>%
    set("branches_k_color", value = select_extremes(c(dark_colors, "grey31"), k_clusters), k = k_clusters) %>%
    set("branches_lwd", 0.7) %>%
    set("labels_cex", 0.6) %>%
    set("labels_colors", value = select_extremes(c(dark_colors, "grey31"), k_clusters), k = k_clusters) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 0.5)

  ggd1 <- as.ggdend(dend)
  p_TriFunQ_clusters_dend_valid <- ggplot(ggd1, horiz = FALSE) +
    labs(title = "Dendrogram - TriFunQ validation", x = "Case", y = "Dissimilarity") +
    theme_plot() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  print(p_TriFunQ_clusters_dend_valid)
}

# C. Projection scatter plots
# Training scatter
df_projected_train <- data.frame(
  Target = clusters_train,
  valid_result_projection$projection_object$embedding
)

p_training_clusters_scatter <- ggplot(df_projected_train,
                                      aes(x = X1, y = X2, shape = as.factor(Target),
                                          fill = as.factor(Target), color = as.factor(Target))) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(level = 0.95, aes(fill = as.factor(Target)), alpha = 0.3, geom = "polygon") +
  theme_plot() +
  labs(title = "Training Clusters in Projection Space",
       x = "Dim 1", y = "Dim 2", shape = "Cluster") +
  scale_color_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  theme(legend.position.inside = TRUE, legend.position = c(.1, .2)) +
  guides(color = "none", fill = "none")

print(p_training_clusters_scatter)

# Validation scatter
df_projected_val <- data.frame(
  Target = clusters_val,
  valid_result_projection$Projected
)

p_validation_clusters_scatter <- ggplot(df_projected_val,
                                        aes(x = Dim1, y = Dim2, shape = as.factor(Target),
                                            fill = as.factor(Target), color = as.factor(Target))) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(level = 0.95, aes(fill = as.factor(Target)), alpha = 0.3, geom = "polygon") +
  theme_plot() +
  labs(title = "Validation Clusters in Projection Space",
       x = "Dim 1", y = "Dim 2", shape = "Cluster") +
  scale_color_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  theme(legend.position.inside = TRUE, legend.position = c(.1, .2)) +
  guides(color = "none", fill = "none")

print(p_validation_clusters_scatter)

# D. Silhouette width comparison
sil_summary <- data.frame(
  Solution = c("Training TriFunQ", "Validation"),
  Avg_Silhouette = c(mean(sil_train[, "sil_width"]), mean(sil_val[, "sil_width"]))
)

p_silhouette_comparison <- ggplot(sil_summary, aes(x = Solution, y = Avg_Silhouette, fill = Solution)) +
  geom_col(alpha = 0.8, color = "black", width = 0.5) +
  geom_text(aes(label = round(Avg_Silhouette, 3)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = select_extremes(c(dark_colors, "grey31"), k_clusters)) +
  labs(title = "Average Silhouette Width Comparison", y = "Average Silhouette Width") +
  theme_plot() +
  theme(legend.position = "none") +
  ylim(0, NA)

print(p_silhouette_comparison)

# E. Combined diagnostic plot (training + validation)
cat("Creating combined diagnostic plot...\n")

# Top row: training diagnostics
top_r <- (p_TriFunQ_clusters_dend + theme(plot.margin = margin(5, 5, 5, 5))) |
  (p_TriFunQ_clusters_silhouette + theme(plot.margin = margin(5, 5, 5, 5))) |
  (p_TriFunQ_clusters_scatter + theme(plot.margin = margin(5, 5, 5, 5)))

# Bottom row: validation diagnostics
bottom_r <- (p_TriFunQ_clusters_dend_valid + theme(plot.margin = margin(5, 5, 5, 5))) |
  (p_validation_clusters_silhouette + theme(plot.margin = margin(5, 5, 5, 5))) |
  (p_validation_clusters_scatter + theme(plot.margin = margin(5, 5, 5, 5)))

# Combine with labels
combined_trig_clustering_plot <- (top_r / bottom_r) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Predictive modeling: Best clustering solution",
    subtitle = "Upper row: Training | Lower row: Validation",
    tag_levels = LETTERS,
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

print(combined_trig_clustering_plot)

# Save combined diagnostic plot
ggsave("combined_trig_clustering_plot.svg", combined_trig_clustering_plot,
       width = 18, height = 12, dpi = 300, bg = "white")
ggsave("combined_trig_clustering_plot.png", combined_trig_clustering_plot,
       width = 18, height = 12, dpi = 300, bg = "white")

cat("Combined diagnostic plot saved\n")

# ============================================================================ #
# 15. EFFECT SIZE FUNCTIONS (RANK-BISERIAL AND ETA-SQUARED)
# ============================================================================ #

cat("\n=== Defining Effect Size Functions ===\n")

#' Calculate Rank-Biserial Correlation
#'
#' Non-parametric effect size for two-group comparisons, companion to
#' Mann-Whitney U / Wilcoxon rank-sum test.
#'
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @return Rank-biserial correlation coefficient (range: -1 to 1)
#' @details
#' Formula: r = 1 - (2U)/(n1*n2) where U is the Mann-Whitney U statistic.
#' Interpretation: |r| = 0.1 (small), 0.3 (medium), 0.5 (large)
rank_biserial_val <- function(x, y) {
  wtest <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wtest$statistic)
  n1 <- length(x)
  n2 <- length(y)
  r <- 1 - (2 * U) / (n1 * n2)
  return(r)
}

#' Bootstrap Confidence Intervals for Rank-Biserial Correlation
#'
#' Generates 95% confidence intervals using percentile bootstrap method.
#'
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @param R Number of bootstrap iterations (default: 1000)
#' @return Vector with lower and upper 95% CI bounds
boot_rank_biserial <- function(x, y, R = 1000) {
  data <- data.frame(group = rep(c("x", "y"), times = c(length(x), length(y))),
                     value = c(x, y))

  stat_fun <- function(data, indices) {
    d_x <- data$value[indices][data$group[indices] == "x"]
    d_y <- data$value[indices][data$group[indices] == "y"]
    if (length(d_x) < 2 || length(d_y) < 2) return(NA)
    return(rank_biserial_val(d_x, d_y))
  }

  boot_out <- boot(data, stat_fun, R = R)
  ci <- boot.ci(boot_out, type = "perc")

  if (is.null(ci)) {
    return(c(NA, NA))
  } else {
    return(ci$percent[4:5])
  }
}

#' Calculate Eta-Squared for Kruskal-Wallis Test
#'
#' Non-parametric effect size for multi-group comparisons (>2 groups).
#'
#' @param data Data frame with 'value' and 'Cluster' columns
#' @return Eta-squared value (range: 0 to 1)
#' @details
#' Formula: η² = (H - k + 1) / (n - k) where H is Kruskal-Wallis H statistic.
#' Interpretation: η² = 0.01 (small), 0.06 (medium), 0.14 (large)
eta_squared <- function(data) {
  kw <- kruskal.test(value ~ Cluster, data = data)
  H <- kw$statistic
  k <- length(unique(data$Cluster))
  n <- nrow(data)
  eta2 <- (H - k + 1) / (n - k)
  return(eta2)
}

#' Bootstrap Confidence Intervals for Eta-Squared
#'
#' Generates 95% confidence intervals using percentile bootstrap method.
#'
#' @param data Data frame with 'value' and 'Cluster' columns
#' @param R Number of bootstrap iterations (default: 200)
#' @return Vector with lower and upper 95% CI bounds
boot_eta_squared <- function(data, R = 200) {
  if (nrow(data) < 20) return(c(NA, NA))

  n <- nrow(data)
  boot_values <- numeric(R)

  for (i in 1:R) {
    indices <- sample(n, replace = TRUE)
    boot_data <- data[indices, ]
    if (nrow(boot_data) >= 10) {
      boot_values[i] <- eta_squared(boot_data)
    }
  }

  boot_values <- boot_values[!is.na(boot_values)]
  if (length(boot_values) < 50) return(c(NA, NA))

  ci_lower <- quantile(boot_values, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_values, 0.975, na.rm = TRUE)
  return(c(ci_lower, ci_upper))
}

cat("Effect size functions defined\n")

# ============================================================================ #
# 16. PREPARE DATA FOR VARIABLE ANALYSIS BY CLUSTER
# ============================================================================ #

cat("\n=== Preparing Data for Variable Analysis ===\n")

# Training data: combine clusters with original variables
trigeminal_clustered_training_data <- cbind.data.frame(
  Cluster = as.factor(TriFunQ_clusters),
  trigeminale_training_data[, names(trigeminale_training_data) %in% c(
    variables_by_categories$Nasal_chemosensory_perception,
    variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
)

trigeminal_clustered_training_data$ID <- trigeminale_training_data$ID

# Transform CO2 threshold (exclude early dataset, negate for direction)
trigeminal_clustered_training_data$`CO2 threshold`[trigeminal_clustered_training_data$ID <= 549] <- NA
trigeminal_clustered_training_data$`CO2 threshold` <-
  -slog(trigeminal_clustered_training_data$`CO2 threshold`)

# Transform AmmoLa intensity
trigeminal_clustered_training_data$`AmmoLa intensity` <-
  reflect_slog_unflipped(trigeminal_clustered_training_data$`AmmoLa intensity`)

# Rename for clarity
trigeminal_clustered_training_data <- trigeminal_clustered_training_data %>%
  rename(AmmoLa_intensity_transformed = `AmmoLa intensity`,
         CO2_threshold_transformed = `CO2 threshold`)

# Convert to long format for plotting
trigeminal_clustered_training_data_long <- reshape2::melt(
  trigeminal_clustered_training_data[, !names(trigeminal_clustered_training_data) %in% c("ID")],
  id.vars = "Cluster"
)

# Validation data: combine clusters with original variables
trigeminal_clustered_validation_data <- cbind.data.frame(
  Cluster = as.factor(clusters_val),
  trigeminale_validation_data[, names(trigeminale_validation_data) %in% c(
    variables_by_categories$Nasal_chemosensory_perception,
    variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
)

trigeminal_clustered_validation_data$ID <- trigeminale_validation_data$ID

# Transform variables (same as training)
trigeminal_clustered_validation_data$`CO2 threshold`[trigeminal_clustered_validation_data$ID <= 549] <- NA
trigeminal_clustered_validation_data$`CO2 threshold` <-
  -slog(trigeminal_clustered_validation_data$`CO2 threshold`)
trigeminal_clustered_validation_data$`AmmoLa intensity` <-
  reflect_slog_unflipped(trigeminal_clustered_validation_data$`AmmoLa intensity`)

trigeminal_clustered_validation_data <- trigeminal_clustered_validation_data %>%
  rename(AmmoLa_intensity_transformed = `AmmoLa intensity`,
         CO2_threshold_transformed = `CO2 threshold`)

# Convert to long format
trigeminal_clustered_validation_data_long <- reshape2::melt(
  trigeminal_clustered_validation_data[, !names(trigeminal_clustered_validation_data) %in% c("ID")],
  id.vars = "Cluster"
)

# Prepare dataset list
datasets <- list(
  training = trigeminal_clustered_training_data_long,
  valid = trigeminal_clustered_validation_data_long
)

cat("Data prepared for", length(unique(trigeminal_clustered_training_data_long$variable)),
    "variables across", k_clusters, "clusters\n")

# Common plotting settings
dodge_width <- 0.8
dark_colors <- c("#DDD6B2", "#D2CC9F", "#C8C28C", "#BEB879")

# To be continued...
# ============================================================================ #
# ============================================================================ #
# 17. VISUALIZATION FUNCTIONS FOR VARIABLE ANALYSIS
# ============================================================================ #

cat("\n=== Defining Visualization Functions ===\n")

#' Create Violin Plot with Cluster Comparisons
#'
#' Generates violin/box/jitter plots for each variable, comparing values
#' across clusters. Includes statistical comparison (Kruskal-Wallis p-value)
#' and median trend lines.
#'
#' @param data Long-format data frame with 'variable', 'value', 'Cluster' columns
#' @param dataset_name Name of dataset for plot title
#' @param k_clusters Number of clusters (for positioning calculations)
#' @return ggplot object
create_violin_plot <- function(data, dataset_name, k_clusters) {
  dodge_width <- 0.8

  # Calculate medians for trend lines
  median_data <- data %>%
    group_by(variable, Cluster) %>%
    dplyr::summarise(median_value = mean(value, na.rm = TRUE), .groups = "drop")

  # Generate dodged x-positions for median lines (0.8, 1.2, 1.6, etc.)
  positions <- seq(0.8, 0.8 + (k_clusters - 1) * 0.4, by = 0.4)
  median_data$x_position <- rep(positions, times = n_distinct(data$variable))

  ggplot(data, aes(x = variable, y = value, color = Cluster, fill = Cluster)) +
    geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = dodge_width)) +
    geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = dodge_width), outlier.shape = NA) +
    geom_jitter(alpha = 1, size = 0.3,
                position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = dodge_width)) +
    ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                               label.y.npc = "top", vjust = -0.2, size = 3) +
    geom_line(data = median_data, aes(x = x_position, y = median_value, group = variable),
              color = "black", size = 0.8, linetype = "dashed", inherit.aes = FALSE) +
    facet_wrap(variable ~ ., nrow = 3, scales = "free", labeller = label_wrap_gen(width = 20)) +
    scale_color_manual(values = dark_colors) +
    scale_fill_manual(values = dark_colors) +
    labs(title = paste("Raw trigeminal data per variable and cluster -", dataset_name),
         fill = "Cluster", x = NULL) +
    theme_plot() +
    theme(
      legend.direction = "horizontal",
      legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black", size = 6, face = "plain"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(color = "none")
}

#' Create Rank-Biserial Effect Size Plot
#'
#' For 2-cluster comparisons: calculates and plots rank-biserial correlation
#' with bootstrap confidence intervals for each variable.
#'
#' @param data Long-format data frame
#' @param dataset_name Name of dataset for plot title
#' @return ggplot object
create_rank_biserial_plot <- function(data, dataset_name) {
  variables <- unique(data$variable)

  # Calculate rank-biserial for each variable
  results_rb <- map_df(variables, function(var) {
    data_var <- data %>% filter(variable == var, !is.na(value))

    x <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[1]]
    y <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[2]]

    if (length(x) < 2 || length(y) < 2) {
      return(tibble(variable = var, r = NA, ci_lower = NA, ci_upper = NA))
    }

    r <- rank_biserial_val(x, y)
    ci <- boot_rank_biserial(x, y)

    tibble(variable = var, r = r, ci_lower = ci[1], ci_upper = ci[2])
  })

  # Create bar plot with error bars
  ggplot(results_rb, aes(x = reorder(variable, r), y = r)) +
    geom_col(fill = "cornsilk3") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    coord_flip() +
    labs(title = paste("Rank-biserial correlation (95% CI) -", dataset_name),
         y = "Rank-biserial r", x = "Variables") +
    theme_plot() +
    theme(axis.text.y = element_text(size = 8)) +
    geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = "dashed",
               color = c("salmon", "grey55", "salmon")) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35))
}

#' Create Eta-Squared Effect Size Plot
#'
#' For >2 cluster comparisons: calculates and plots eta-squared (η²) from
#' Kruskal-Wallis test with bootstrap confidence intervals for each variable.
#'
#' @param data Long-format data frame
#' @param dataset_name Name of dataset for plot title
#' @return ggplot object
create_eta_plot <- function(data, dataset_name) {
  variables <- unique(data$variable)

  # Calculate eta-squared for each variable
  results_eta <- map_df(variables, function(var) {
    data_var <- data %>% filter(variable == var, !is.na(value))
    if (nrow(data_var) < 10 || length(unique(data_var$Cluster)) < 2) {
      return(tibble(variable = var, eta2 = NA, ci_lower = NA, ci_upper = NA))
    }
    eta_val <- eta_squared(data_var)
    ci <- boot_eta_squared(data_var, R = 100)
    tibble(variable = var, eta2 = eta_val, ci_lower = ci[1], ci_upper = ci[2])
  })

  # Create bar plot with error bars
  ggplot(results_eta, aes(x = reorder(variable, eta2), y = eta2)) +
    geom_col(fill = "cornsilk3") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    coord_flip() +
    labs(title = paste("Kruskal-Wallis η² effect size (95% CI) -", dataset_name),
         y = "η² (eta-squared)", x = "Variables") +
    theme_plot() +
    theme(axis.text.y = element_text(size = 8)) +
    geom_hline(yintercept = c(0, 0.06, 0.14), linetype = "dashed",
               color = c("grey55", "orange", "salmon")) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35))
}

cat("Visualization functions defined\n")

# ============================================================================ #
# 18. GENERATE PLOTS AND EFFECT SIZE TABLES
# ============================================================================ #

cat("\n=== Generating Variable Analysis Plots ===\n")
cat("Effect size method:", ifelse(k_clusters == 2, "Rank-biserial (2 groups)", "Eta-squared (>2 groups)"), "\n")

# Generate all plots and tables for both datasets
plots_and_tables <- map2(datasets, names(datasets), function(data, name) {

  # Create violin plot
  violin_plot <- create_violin_plot(data, stringr::str_to_title(name), k_clusters = k_clusters)

  # Determine effect size based on number of clusters
  n_clusters <- length(unique(data$Cluster))

  if (n_clusters == 2) {
    # Two clusters: use rank-biserial correlation
    effect_plot <- create_rank_biserial_plot(data, stringr::str_to_title(name))

    variables <- unique(data$variable)
    effect_results <- map_df(variables, function(var) {
      data_var <- data %>% filter(variable == var, !is.na(value))

      x <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[1]]
      y <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[2]]

      if (length(x) < 2 || length(y) < 2) {
        return(tibble(variable = var, r = NA, ci_lower = NA, ci_upper = NA))
      }

      r <- rank_biserial_val(x, y)
      ci <- boot_rank_biserial(x, y)

      tibble(variable = var, r = r, ci_lower = ci[1], ci_upper = ci[2])
    })

  } else {
    # More than two clusters: use eta-squared
    effect_plot <- create_eta_plot(data, stringr::str_to_title(name))

    variables <- unique(data$variable)
    effect_results <- map_df(variables, function(var) {
      data_var <- data %>% filter(variable == var, !is.na(value))
      if (nrow(data_var) < 10 || length(unique(data_var$Cluster)) < 2) {
        return(tibble(variable = var, eta2 = NA, ci_lower = NA, ci_upper = NA))
      }
      eta_val <- eta_squared(data_var)
      ci <- boot_eta_squared(data_var)
      tibble(variable = var, eta2 = eta_val, ci_lower = ci[1], ci_upper = ci[2])
    })
  }

  list(
    violin = violin_plot,
    effect_plot = effect_plot,
    effect_table = effect_results
  )
})

# Extract results
p_trigeminal_clustered_training_data <- plots_and_tables$training$violin
p_trigeminal_clustered_valid_data <- plots_and_tables$valid$violin
p_effect_training <- plots_and_tables$training$effect_plot
p_effect_valid <- plots_and_tables$valid$effect_plot
results_effect_training <- plots_and_tables$training$effect_table
results_effect_valid <- plots_and_tables$valid$effect_table

# Display plots
print(p_trigeminal_clustered_training_data)
print(p_trigeminal_clustered_valid_data)
print(p_effect_training)
print(p_effect_valid)

# ============================================================================ #
# 19. COMBINED PLOTS FOR TRAINING DATA
# ============================================================================ #

cat("\n=== Creating Combined Training Plot ===\n")

# Combine effect size plot (left) + violin plot (right)
p_left_clusters_train <- p_effect_training + theme(plot.margin = margin(5, 5, 5, 5))
p_right_clusters_train <- p_trigeminal_clustered_training_data + theme(plot.margin = margin(5, 5, 5, 5))

p_combined_clusters_train <- (p_left_clusters_train | p_right_clusters_train) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Cluster differences in trigeminal measures: Training data",
    tag_levels = list(c("A", "B")),
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

print(p_combined_clusters_train)

# Save combined training plot
ggsave("p_combined_clusters_train.svg", p_combined_clusters_train,
       width = 17, height = 15, dpi = 300, limitsize = FALSE)
ggsave("p_combined_clusters_train.png", p_combined_clusters_train,
       width = 17, height = 15, dpi = 300, limitsize = FALSE)

cat("Combined training plot saved\n")

# ============================================================================ #
# 20. STATISTICAL TESTS AND CROSS-DATASET CORRELATION
# ============================================================================ #

cat("\n=== Calculating P-Values and Cross-Dataset Correlation ===\n")

# Print effect size tables
cat("\nTraining effect sizes:\n")
print(as.data.frame(results_effect_training))
cat("\nValidation effect sizes:\n")
print(as.data.frame(results_effect_valid))

# Calculate p-values (Wilcoxon test for 2 groups, Kruskal-Wallis for >2)
cat("\nCalculating p-values...\n")
p_vals_train <- apply(trigeminal_clustered_training_data[, !names(trigeminal_clustered_training_data) %in% c("Cluster", "ID")], 2,
                      function(x) wilcox.test(x ~ trigeminal_clustered_training_data$Cluster)$p.value)
cat("Training p-values:\n")
print(p_vals_train)

p_vals_valid <- apply(trigeminal_clustered_validation_data[, !names(trigeminal_clustered_validation_data) %in% c("Cluster", "ID")], 2,
                      function(x) wilcox.test(x ~ trigeminal_clustered_validation_data$Cluster)$p.value)
cat("Validation p-values:\n")
print(p_vals_valid)

# Sort both tables by variable (ensure same order)
results_effect_training <- results_effect_training %>% arrange(variable)
results_effect_valid <- results_effect_valid %>% arrange(variable)

# Kendall's tau correlation between training and validation effect sizes
cat("\nCorrelating effect sizes across datasets (Kendall's tau)...\n")
if (k_clusters == 2) {
  # Rank-biserial correlation
  plot(results_effect_valid$r ~ results_effect_training$r)
  Tau_clusters <- cor.test(x = results_effect_training$r,
                           y = results_effect_valid$r, method = "kendall")
  df_effsize <- cbind.data.frame(Training = results_effect_training$r,
                          Validation = results_effect_valid$r)
} else {
  # Eta-squared
  plot(results_effect_valid$eta2 ~ results_effect_training$eta2)
  Tau_clusters <- cor.test(x = results_effect_training$eta2,
                           y = results_effect_valid$eta2, method = "kendall")
  df_effsize <- data.frame(Training = results_effect_training$eta2,
                          Validation = results_effect_valid$eta2)
}
df_effsize$mean_effect <- rowMeans(df_effsize)
df_effsize$variable <- str_wrap(results_effect_training$variable, width = 22)

print(Tau_clusters)

# Create correlation plot
p_cor_Effsizes_tau <-
  ggplot(df_effsize, aes(x = Training, y = Validation, color = as.factor(mean_effect))) +
  geom_point(size = 4, show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(x = Training, y = Validation, label=variable), force = TRUE, force_pull = TRUE, size = 1.5, show.legend = FALSE) +
  stat_cor(method = "kendall", label.x = .1, label.y = .6) +
  guides(color = "none") +
  theme_plot()

# ============================================================================ #
# 21. COMBINED PLOTS FOR VALIDATION DATA
# ============================================================================ #

cat("\n=== Creating Combined Validation Plot ===\n")

# Combine effect size + correlation (left) + violin plot (right)
p_left_clusters_valid <- (p_effect_valid + theme(plot.margin = margin(5, 5, 5, 5))) /
  (p_cor_Effsizes_tau + theme(plot.margin = margin(5, 5, 5, 5))) +
  plot_layout(heights = c(2, 1))
p_right_clusters_valid <- p_trigeminal_clustered_valid_data + theme(plot.margin = margin(5, 5, 5, 5))

p_combined_clusters_valid <- (p_left_clusters_valid | p_right_clusters_valid) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Cluster differences in trigeminal measures: Validation data",
    tag_levels = list(c("A", "B")),
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

print(p_combined_clusters_valid)

# Save combined validation plot
ggsave("p_combined_clusters_valid.svg", p_combined_clusters_valid,
       width = 17, height = 15, dpi = 300, limitsize = FALSE)
ggsave("p_combined_clusters_valid.png", p_combined_clusters_valid,
       width = 17, height = 15, dpi = 300, limitsize = FALSE)

cat("Combined validation plot saved\n")

# ============================================================================ #
# 22. CLUSTERED VARIABLES HEATMAP
# ============================================================================ #

cat("\n=== Creating Clustered Variables Heatmap ===\n")

# Prepare data for heatmap
row_means_all <- rowMeans(rbind(heat_matrix, heat_matrix_validation))
row_order <- order(row_means_all)
clusters_all <- c(TriFunQ_clusters, clusters_val)

# Combined matrix (rows ordered by mean)
heat_matrix_main <- rbind(heat_matrix, heat_matrix_validation)[row_order, ]

# Dataset annotation (training vs validation)
train_valid_anno <- rep(c("training", "validation"),
                       c(nrow(heat_matrix), nrow(heat_matrix_validation)))[row_order]

# Define color palette
pal2 <- colorRampPalette(c("cornsilk", "cornsilk3", "cornsilk4", "grey33"))

# Right annotation: row means barplot (colored by dataset)
row_means_ha <- rowAnnotation(
  bar = anno_barplot(row_means_all[row_order],
                     gp = gpar(fill = ifelse(train_valid_anno == "training", "cornsilk1", "cornsilk4"),
                               col = NA),
                     border = TRUE),
  show_legend = FALSE,
  annotation_label = "Row means",
  annotation_name_rot = 90,
  show_annotation_name = TRUE,
  width = unit(2, "cm")
)

# Left annotation: dataset label
train_valid_ha <- rowAnnotation(
  Dataset = train_valid_anno,
  col = list(Dataset = c("training" = "cornsilk1", "validation" = "cornsilk4")),
  show_annotation_name = TRUE,
  annotation_name_rot = 90,
  width = unit(1.5, "cm")
)

# Create heatmap function
create_heatmap_trig_clustered_data <- function() {
  grid::grid.grabExpr({
    ht <- ComplexHeatmap::Heatmap(
      as.matrix(heat_matrix_main),
      left_annotation = train_valid_ha,
      right_annotation = row_means_ha,
      col = pal2(112),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 90,
      heatmap_legend_param = list(
        title = "Scaled value [0,3]",
        direction = "horizontal",
        title_position = "lefttop"
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      column_names_max_height = grid::unit(16, "cm"),
      row_split = clusters_all[row_order]
    )
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  })
}

# Generate heatmap
gp_clustered_variables <- create_heatmap_trig_clustered_data()

# Create title plot
p_title_heat <- ggplot() +
  labs(title = "Clustered variables") +
  theme_void() +
  theme(
    plot.title = element_text(face = "plain", size = 12, color = "#222222",
                             hjust = 0, margin = margin(b = 5)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 0, 20)
  )

# Combine title + heatmap
heat_plot_clustered_variables <- cowplot::plot_grid(
  p_title_heat,
  cowplot::as_grob(gp_clustered_variables),
  ncol = 1,
  rel_heights = c(0.05, 1),
  align = "v"
)

print(heat_plot_clustered_variables)

# Save heatmap
ggsave("heat_plot_clustered_variables.svg", heat_plot_clustered_variables,
       width = 12, height = 15)
ggsave("heat_plot_clustered_variables.png", heat_plot_clustered_variables,
       width = 12, height = 15, dpi = 300)

cat("Clustered variables heatmap saved\n")

# ============================================================================ #
# 23. END OF ANALYSIS
# ============================================================================ #

cat("\n" , rep("=", 78), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 78), "\n\n", sep = "")

cat("Summary:\n")
cat("  - Best clustering method:", best_projection_method, "+", best_clustering_method, "\n")
cat("  - Number of clusters:", k_clusters, "\n")
cat("  - Training samples:", length(TriFunQ_clusters), "\n")
cat("  - Validation samples:", length(clusters_val), "\n")
cat("  - Effect size method:", ifelse(k_clusters == 2, "Rank-biserial", "Eta-squared"), "\n")
cat("  - Cross-dataset correlation (Kendall's tau):", round(Tau_clusters$estimate, 3), "\n")
cat("  - Correlation p-value:", format.pval(Tau_clusters$p.value, digits = 3), "\n\n")

cat("Output files generated:\n")
cat("  - Combined_projection_and_clustering_analysis_plot_*.svg\n")
cat("  - dfClusterQuality.csv\n")
cat("  - p_GMM_row_means.svg/png\n")
cat("  - combined_trig_clustering_plot.svg/png\n")
cat("  - p_combined_clusters_train.svg/png\n")
cat("  - p_combined_clusters_valid.svg/png\n")
cat("  - heat_plot_clustered_variables.svg/png\n\n")

cat("Script completed successfully!\n")
cat(date(), "\n\n")

################################################################################
# END OF SCRIPT
################################################################################
