################################################################################
# Script: Trigeminal Cluster Analysis and Visualization
# Author: Jorn Lotsch
# Date: 2025-11-04
# Purpose:
#   Process and interpret clustered trigeminal sensitivity data.
#   Includes loading precomputed clustering results, summarizing cluster-wise
#   variable means, selecting relevant variables via effect size (eta²) and
#   feature importance (Boruta), and generating heatmaps for cluster profiling.
#
# Main components:
#   - Clustered trigeminal dataset (training subset)
#   - Heatmap input matrix for cluster characterization
#   - Effect size estimates (eta²)
#   - Boruta feature importance scores
################################################################################

# ==============================================================================#
# Libraries and helper functions
# ==============================================================================#

# Load global settings and utility functions
source("globals.R")
source("utils.R")


# ==============================================================================#
# Load data
# ==============================================================================#

cat("\\n=== Loading Data ===\\n")

# Load clustered trigeminal dataset
all_trig_clustered_data_and_clusters <- read.csv(
  "trigeminal_clustered_data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Ensure correct data types
all_trig_clustered_data_and_clusters$ID <- as.character(all_trig_clustered_data_and_clusters$ID)
all_trig_clustered_data_and_clusters$Set <- as.character(all_trig_clustered_data_and_clusters$Set)
all_trig_clustered_data_and_clusters$Cluster <- as.factor(all_trig_clustered_data_and_clusters$Cluster)

# Subset to training data only
trig_clusters_training_data <- subset(all_trig_clustered_data_and_clusters, Set == "Training")

# Load heatmap input matrix
heat_matrix_clustered <- read.csv(
  "heat_matrix_trig_for_clustering.csv",
  row.names = 1,
  check.names = FALSE
)

# Load effect size and feature importance results
eta_square_training <- read.csv(
  "effsizes_trig_cluster_variables.csv",
  row.names = 1
)

Boruta_training <- read.csv(
  "Boruta_importance_trigeminal_sensitivity_clusters.csv",
  row.names = 1
)


# ==============================================================================#
# Cluster interpretation
# ==============================================================================#

# Combine cluster labels with feature matrix
clusters_and_clustered_data <- cbind.data.frame(
  Cluster = trig_clusters_training_data$Cluster,
  heat_matrix_clustered
)

# Compute mean values per cluster and reshape for heatmap
mat_clusterdata <- clusters_and_clustered_data %>%
  group_by(Cluster) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-Cluster, names_to = "Variable", values_to = "Mean") %>%
  pivot_wider(names_from = Cluster, values_from = Mean) %>%
  column_to_rownames("Variable") %>%
  as.matrix()

# Define color scale for heatmap
col_fun <- colorRamp2(
  c(
    min(mat_clusterdata, na.rm = TRUE),
    mean(mat_clusterdata, na.rm = TRUE),
    max(mat_clusterdata, na.rm = TRUE)
  ),
  c(actual_palette[4], "ghostwhite", actual_palette[13])
)

# Perform cABC analysis on effect sizes
cABC_trig_var_effectsize <- cABCanalysis::cABC_analysis(
  results_effect_training$eta2,
  PlotIt = TRUE
)

# Select variables based on cABC results
vars_to_keep <- results_effect_training$variable[cABC_trig_var_effectsize$Aind]

# Subset matrix to selected variables
mat_clusterdata_1 <- mat_clusterdata[
  rownames(mat_clusterdata) %in% vars_to_keep,
]

# Generate heatmap of cluster means
ComplexHeatmap::Heatmap(
  mat_clusterdata,
  name = "Mean",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Clusters",
  row_title = "Variables"
)

# Add effect size and Boruta importance to matrix
mat_clusterdata_and_eff_sizes <- as.data.frame(mat_clusterdata)

# Match eta² values
idx <- match(
  rownames(mat_clusterdata_and_eff_sizes),
  eta_square_training$variable_orig
)
mat_clusterdata_and_eff_sizes$Eta_square <- eta_square_training$Training[idx]

# Match Boruta importance values
idx <- match(
  rownames(mat_clusterdata_and_eff_sizes),
  Boruta_training$Feature
)
mat_clusterdata_and_eff_sizes$Boruta <- Boruta_training$median[idx]


# ==============================================================================#
# Confusion matrix analysis
# ==============================================================================#

# ------------------------------------------------------------------------------
# Load feature selection results
# ------------------------------------------------------------------------------
clusters_trig_modulators_sig_predictors_interpretation <- read.csv(
  "clusters_trig_modulators_coef_table_full.csv",
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# Filter and aggregate selected variables across clusters
# ------------------------------------------------------------------------------
clusters_trig_modulators_sig_predictors_interpretation_1 <-
  clusters_trig_modulators_sig_predictors_interpretation %>%
  # Keep clusters 1–4 only
  filter(class %in% 1:4) %>%
  # Aggregate selection flags per variable
  group_by(variable) %>%
  summarise(
    glm_selected = any(glm_selected, na.rm = TRUE),
    ridge_selected = any(ridge_selected, na.rm = TRUE),
    lasso_selected = any(lasso_selected, na.rm = TRUE),
    elastic_selected = any(elastic_selected, na.rm = TRUE),

    # Mark Boruta as confirmed if selected in any cluster
    Boruta = if (any(Boruta == "Confirmed", na.rm = TRUE)) {
      "Confirmed"
    } else {
      "Rejected"
    },
    .groups = "drop"
  ) %>%
  # Keep variables selected by at least one method
  filter(
    glm_selected | ridge_selected | lasso_selected |
      elastic_selected | Boruta == "Confirmed"
  )

# ------------------------------------------------------------------------------
# Count number of selected features per method
# ------------------------------------------------------------------------------
n_features_glm <- sum(clusters_trig_modulators_sig_predictors_interpretation_1$glm_selected, na.rm = TRUE)
n_features_ridge <- sum(clusters_trig_modulators_sig_predictors_interpretation_1$ridge_selected, na.rm = TRUE)
n_features_lasso <- sum(clusters_trig_modulators_sig_predictors_interpretation_1$lasso_selected, na.rm = TRUE)
n_features_elastic <- sum(clusters_trig_modulators_sig_predictors_interpretation_1$elastic_selected, na.rm = TRUE)
n_features_RF <- sum(clusters_trig_modulators_sig_predictors_interpretation_1$Boruta == "Confirmed", na.rm = TRUE)

# Extract Boruta-confirmed features
Boruta_features <- clusters_trig_modulators_sig_predictors_interpretation_1$variable[
  clusters_trig_modulators_sig_predictors_interpretation_1$Boruta == "Confirmed"
]

cat("\\nGenerating confusion matrix analysis...\\n")

# ------------------------------------------------------------------------------
# Generate confusion matrix plots for selected categorical variables
# ------------------------------------------------------------------------------
cv_plots <- lapply(c("Gender_m", "turbinate.surgery"), function(actual_feature) {
  # Prepare dataset with feature and cluster labels
  df_cm_cvms_clus <- data.frame(
    Gender = trig_train_final[actual_feature],
    Cluster = trig_train_final$Cluster,
    check.names = FALSE
  )

  # Remove missing values
  df_cm_cvms_clus <- df_cm_cvms_clus[complete.cases(df_cm_cvms_clus), ]

  # Convert binary feature to factor
  df_cm_cvms_clus[[actual_feature]] <- factor(
    df_cm_cvms_clus[[actual_feature]],
    levels = c(0, 1),
    labels = na.omit(unique(df_cm_cvms_clus[[actual_feature]]))
  )

  # Convert cluster labels to factor
  df_cm_cvms_clus$Cluster <- factor(
    df_cm_cvms_clus$Cluster,
    levels = c(1, 2, 3, 4),
    labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
  )

  # Create contingency table
  tbl_clus <- table(
    df_cm_cvms_clus[[actual_feature]],
    df_cm_cvms_clus$Cluster
  )

  # Perform Fisher's exact test
  fisher_result_clus <- fisher.test(tbl_clus)
  fisher_p_clus <- signif(fisher_result_clus$p.value, 3)

  # Convert table to data frame for plotting
  plot_df <- as.data.frame(tbl_clus)
  colnames(plot_df) <- c("Gender", "Cluster", "Freq")

  # Create heatmap-style confusion matrix plot
  p_cm_clus_vs_var <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Cluster, y = Gender, fill = Freq)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 5) +
    ggplot2::scale_fill_gradient(low = "lightcyan", high = "steelblue") +
    ggplot2::labs(
      title = paste0(actual_feature, " vs Cluster"),
      subtitle = paste0("Fisher's exact test p-value: ", fisher_p_clus),
      x = "Cluster",
      y = actual_feature,
      fill = "Count"
    ) +
    theme_plot()

  # Display plot
  print(p_cm_clus_vs_var)

  # Return results
  return(list(
    p = p_cm_clus_vs_var,
    plot_df = plot_df,
    fisher_result_clus = fisher_result_clus
  ))
})

# Assign names to plot list
names(cv_plots) <- c("Gender_m", "turbinate surgery")


# ------------------------------------------------------------------------------
# Analysis of continuous (gradual) variables across clusters
# ------------------------------------------------------------------------------
Boruta_gradual_features <- make.names(c(
  "Current smell ability",
  "Sensitivity of the nose to stinging/burning stimuli",
  "Smell ability before COVID-19"
))

gradual_var_exploration <- lapply(Boruta_gradual_features, function(actual_feature) {
  # Prepare dataset
  df_cm_cvms_clus <- data.frame(
    Gender = trig_train_final[actual_feature],
    Cluster = trig_train_final$Cluster,
    check.names = FALSE
  )

  # Remove missing values
  df_cm_cvms_clus <- df_cm_cvms_clus[complete.cases(df_cm_cvms_clus), ]

  library(dplyr)

  # Perform Kruskal-Wallis test
  kruskal_res <- kruskal.test(
    df_cm_cvms_clus[[actual_feature]] ~ as.factor(df_cm_cvms_clus$Cluster)
  )

  # Compute descriptive statistics per cluster
  desc_stats <- df_cm_cvms_clus %>%
    group_by(Cluster) %>%
    summarise(
      n = sum(!is.na(.data[[actual_feature]])),
      mean = mean(.data[[actual_feature]], na.rm = TRUE),
      sd = sd(.data[[actual_feature]], na.rm = TRUE),
      median = median(.data[[actual_feature]], na.rm = TRUE),
      IQR = IQR(.data[[actual_feature]], na.rm = TRUE),
      min = min(.data[[actual_feature]], na.rm = TRUE),
      max = max(.data[[actual_feature]], na.rm = TRUE),
      .groups = "drop"
    )

  # Return results
  return(list(
    kruskal_res = kruskal_res,
    plot_df = plot_df,
    desc_stats = desc_stats
  ))
})

# Assign names to results
names(gradual_var_exploration) <- Boruta_gradual_features
