################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Joern Lotsch
# Date: 2025-11-04
# Description: Comprehensive analysis of trigeminal sensitivity study data
#              including preprocessing, transformation, distribution analysis,
#              correlation analysis, and publication-ready visualizations.
#
# Data includes three main measures:
#   - AmmoLa intensity (R28): Ammonia lateralization intensity
#   - Lateralization score: Correct lateralization count (x/20)
#   - CO2 threshold: Carbon dioxide detection threshold
#
# Key analyses:
#   1. Data quality assessment (censoring, missing values)
#   2. Distribution visualization (PDE, histograms, heatmaps)
#   3. Agreement analysis for most sensitive subjects (Fisher's exact tests)
#   4. Correlation analysis with age and sex effects
#   5. Combined visualization for publication
################################################################################

# ========================================================================== #
# 1. LOAD REQUIRED LIBRARIES AND GLOBALS
# ========================================================================== #

library(circlize) # Color mapping for heatmaps
library(ComplexHeatmap) # Advanced heatmap visualization
library(dplyr) # Data manipulation
library(forcats) # Factor handling
library(ggpmisc) # Additional ggplot2 functionality
library(ggplot2) # Data visualization
library(ggthemes) # Extra themes for ggplot2
library(grid) # Grid graphics
library(lubridate) # Date/time handling
library(MASS) # Robust statistical methods
library(purrr) # Functional programming tools
library(reshape2) # Data reshaping
library(scales) # Scale functions for visualization
library(stringr) # String manipulation
library(tidyr) # Data tidying
library(viridis) # Color scales for plots
library(vcd) # Categorical data visualization
library(cvms) # Confusion matrix visualization
library(psych) # Correlation analysis
library(cowplot) # Plot composition and alignment
library(Hmisc) # Statistical tools
library(DataVisualizations) # Pareto Density Estimation
library(NbClust) # Cluster number detection
library(missForest) # Imputation
library(opGMMassessment) # Gaussian mixture analysis
library(FactoMineR) # Factor analysis
library(factoextra) # Factor visualization
library(cluster) # Clustering
library(forcats) # for fct_rev()
library(gridExtra)
library(ggthemes)
library(ComplexHeatmap)
library(ggplotify)

source("globals.R")

# ========================================================================== #
# 2. GLOBAL OPTIONS AND ANALYSIS SWITCHES
# ========================================================================== #

# Control analysis behavior
remove_censored <- FALSE # Whether to exclude censored values
scale_0_100 <- FALSE # Whether to scale all measures to 0-100
analyze_only_untransformed <- FALSE # Skip transformation analysis
plot_only_untransformed <- TRUE # Show only original data in plots

# ========================================================================== #
# 3. HELPER FUNCTIONS
# ========================================================================== #

#' Sign-preserving logarithmic transformation (zero-invariant)
#'
#' @param x Numeric vector to transform
#' @param base Logarithm base (0 = natural log, 2, 10, or custom)
#' @return Transformed vector maintaining sign of original values
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

#' Inverse slog transformation
inv_slog <- function(y, base = 10) {
  s <- sign(y)
  absY <- abs(y)

  if (base == 0) {
    # inverse of slog with natural log: x = exp(absY) - 1
    val <- expm1(absY)
  } else if (base == 2) {
    # inverse of slog with log2: x = 2^(absY) - 1
    val <- 2 ^ absY - 1
  } else if (base == 10) {
    # inverse of slog with log10: x = 10^(absY) - 1
    val <- 10 ^ absY - 1
  } else {
    # inverse with custom base: solve log1p(absX)*log(base) = absY => log1p(absX) = absY / log(base)
    val <- expm1(absY / log(base))
  }
  return(s * val)
}


#' Reflected logarithmic transformation
reflect_slog <- function(x) {
  slog(max(x, na.rm = TRUE) + 1 - x)
}

reflect_slog_unflipped <- function(x) {
  -slog(max(x, na.rm = TRUE) + 1 - x)
}

# Inverse of reflect_slog_unflipped transform
inv_reflect_slog_unflipped <- function(y, original_max, base = 10) {
  M <- original_max
  x_original <- M + 1 - inv_slog(-y, base)
  return(x_original)
}

#' Create Pareto Density Estimation plot with ggplot2
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

# ========================================================================== #
# 4. READ DATA
# ========================================================================== #

trigeminale_data <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)

# Extract and rename relevant trigeminal measures
trigeminal_measures <- c("AmmoLa intensity", "Lateralization (x/20)", "CO2 threshold")
trigeminal_measures_data <- trigeminale_data[, trigeminal_measures]
names(trigeminal_measures_data) <- c("AmmoLa_intensity", "Lateralization", "CO2_threshold")

# Report sample sizes for each measure
cat("\nSample sizes per measure:\n")
print(apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x))))

# ========================================================================== #
# 5. INITIAL DATA VISUALIZATION - HEATMAP OVERVIEW
# ========================================================================== #

# Scale all measures to 0-100 for visualization comparability
trigeminal_measures_data_scaled <- trigeminal_measures_data %>%
  mutate(
    Lateralization = scale01minmax(Lateralization, minX = 0, maxX = 20) * 100,
    CO2_threshold = scale01minmax(CO2_threshold, minX = 100, maxX = 2000) * 100,
    Segment = if_else(row_number() <= 549, "first_part", "second_part")
  )

# Reshape data for heatmap visualization
heatmap_data <- trigeminal_measures_data_scaled %>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols = c(AmmoLa_intensity, Lateralization, CO2_threshold),
    names_to = "Measure",
    values_to = "Value"
  )

# Custom ggplot2 theme for publication-quality plots
theme_plot <- function() {
  theme_minimal(base_family = "Libre Franklin") +
    theme(
      plot.title = element_text(
        face = "plain", size = 12, color = "#222222",
        hjust = 0, margin = margin(b = 10)
      ),
      axis.title = element_text(face = "plain", size = 10, color = "#444444"),
      axis.text = element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = element_text(
        size = 8, color = "#888888",
        hjust = 0, margin = margin(t = 10)
      ),
      panel.grid.major.y = element_line(
        color = "#dddddd", linetype = "dashed", size = 0.3
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks.length = unit(5, "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = margin(20, 20, 20, 20),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "#222222")
    )
}

# Create comprehensive heatmap showing all measures across observations
cat("\nGenerating overview heatmap...\n")
p_trigeminal_measures_done <- ggplot(
  heatmap_data,
  aes(x = Row, y = Measure, fill = Value)
) +
  geom_tile() +
  scale_fill_gradient(low = "cornsilk", high = "cornsilk4", na.value = "grey90", name = "Value [%]") +
  theme_plot() +
  labs(
    x = "Observation (Subject #)",
    y = "Measure",
    title = "Trigeminal measures",
    fill = "Value [%]"
  ) +
  geom_vline(xintercept = 549.5, linetype = "dashed", color = "black", linewidth = 1) +
  annotate(
    "rect", xmin = 0, xmax = 548.5, ymin = 1.8, ymax = 2.2,
    fill = "white", alpha = 0.7
  ) +
  annotate(
    "text", x = 275, y = 2,
    label = "Breath not held during CO2 threshold measurement",
    size = 5, color = "black"
  ) +
  annotate(
    "rect", xmin = 550.5, xmax = 1001, ymin = 1.8, ymax = 2.2,
    fill = "white", alpha = 0.7
  ) +
  annotate(
    "text", x = 775, y = 2,
    label = "Breath held during CO2 threshold measurement",
    size = 5, color = "black"
  )

print(p_trigeminal_measures_done)
ggsave("p_trigeminal_measures_done.svg", p_trigeminal_measures_done,
       width = 16, height = 4)

# ========================================================================== #
# 6. DATA QUALITY ASSESSMENT - CENSORING ANALYSIS
# ========================================================================== #

cat("\nAnalyzing data censoring patterns...\n")

# Maximum possible values for each measure (ceiling effects)
max_vals <- c(100, 20, 2000)

# Calculate censoring percentages for full dataset
counts <- apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x)))
censored <- sapply(seq_along(max_vals), function(i) sum(trigeminal_measures_data[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored <- (censored / counts) * 100

cat("\nOverall censoring:\n")
print(percentage_censored)

# Censoring in first subset (rows 1-549, breath not controlled)
subset1 <- trigeminal_measures_data
subset1$CO2_threshold[550:nrow(subset1)] <- NA
counts1 <- apply(subset1, 2, function(x) sum(!is.na(x)))
censored1 <- sapply(seq_along(max_vals), function(i) sum(subset1[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored1 <- (censored1 / counts1) * 100

cat("\nCensoring in subset 1 (breath not controlled):\n")
print(percentage_censored1)

# Censoring in second subset (rows 549+, breath hold protocol)
subset2 <- trigeminal_measures_data
subset2$CO2_threshold[1:549] <- NA
counts2 <- apply(subset2, 2, function(x) sum(!is.na(x)))
censored2 <- sapply(seq_along(max_vals), function(i) sum(subset2[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored2 <- (censored2 / counts2) * 100

cat("\nCensoring in subset 2 (breath controlled):\n")
print(percentage_censored2)

# ========================================================================== #
# 7. GROUP AGREEMENT ANALYSIS - TOP 10% SENSITIVITY
# ========================================================================== #

cat("\nAnalyzing agreement among most sensitive subjects...\n")

# Create binary grouping: most sensitive 10% vs others
trigeminal_measures_data_grouped <- trigeminal_measures_data_scaled

# For AmmoLa and Lateralization: top 10% = most sensitive
trigeminal_measures_data_grouped[, 1:2] <- apply(
  trigeminal_measures_data_grouped[, 1:2],
  2,
  function(x) ifelse(x >= 90, 1, 0)
)

# For CO2: lower threshold = more sensitive (inverted scale)
trigeminal_measures_data_grouped$CO2_threshold <- ifelse(
  trigeminal_measures_data_grouped$CO2_threshold <= 10, 1, 0
)

# Summary of sensitivity thresholding
cat("\nDistribution of most sensitive subjects:\n")
print(apply(trigeminal_measures_data_grouped[, 1:3], 2, table))

# --- Agreement analysis for ALL subjects ---
trigeminal_measures_data_grouped_all <- trigeminal_measures_data_grouped

# Fisher's exact test: AmmoLa vs Lateralization
cat("\nFisher's test: AmmoLa vs Lateralization (all subjects)\n")
fisher_test_result_lat <- fisher.test(
  trigeminal_measures_data_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_grouped_all$Lateralization
)
print(fisher_test_result_lat)
print(table(
  trigeminal_measures_data_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_grouped_all$Lateralization
))

# Fisher's exact test: AmmoLa vs CO2
cat("\nFisher's test: AmmoLa vs CO2 (all subjects)\n")
fisher_test_result_CO2 <- fisher.test(
  trigeminal_measures_data_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_grouped_all$CO2_threshold
)
print(fisher_test_result_CO2)
print(table(
  trigeminal_measures_data_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_grouped_all$CO2_threshold
))

# --- Agreement analysis ONLY for breath hold subjects ---
trigeminal_measures_data_grouped_breath_hold <- trigeminal_measures_data_grouped_all

# Exclude CO2 data from breath-not-controlled group
trigeminal_measures_data_grouped_breath_hold$CO2_threshold[1:549] <- NA

# Fisher's tests for breath hold group only
cat("\nFisher's test: AmmoLa vs Lateralization (breath hold only)\n")
fisher_test_lat_breath_hold <- fisher.test(table(
  trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity,
  trigeminal_measures_data_grouped_breath_hold$Lateralization
))
print(fisher_test_lat_breath_hold)

cat("\nFisher's test: AmmoLa vs CO2 (breath hold only)\n")
fisher_test_CO2_breath_hold <- fisher.test(table(
  trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity,
  trigeminal_measures_data_grouped_breath_hold$CO2_threshold
))
print(fisher_test_CO2_breath_hold)

# Store tables and test results for visualization
table_AmmoLa_vs_CO2_most_sensitive <- table(
  trigeminal_measures_data_grouped_breath_hold[, c(1, 3)]
)
ftest_AmmoLa_vs_CO2_most_sensitive <- fisher.test(table_AmmoLa_vs_CO2_most_sensitive)
cat("\nAmmoLa vs CO2 agreement (breath hold):\n")
print(ftest_AmmoLa_vs_CO2_most_sensitive)

table_AmmoLa_vs_lateralization_most_sensitive <- table(
  trigeminal_measures_data_grouped_breath_hold[, c(1, 2)]
)
ftest_AmmoLa_vs_lateralization_most_sensitive <- fisher.test(
  table_AmmoLa_vs_lateralization_most_sensitive
)
cat("\nAmmoLa vs Lateralization agreement (breath hold):\n")
print(ftest_AmmoLa_vs_lateralization_most_sensitive)

# ========================================================================== #
# 8. CONFUSION MATRIX ANALYSIS
# ========================================================================== #

cat("\nGenerating confusion matrix analysis...\n")

# Prepare data for confusion matrix (complete cases only)
df_cm_cvms <- cbind.data.frame(
  AmmoLa = trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity,
  CO2 = trigeminal_measures_data_grouped_breath_hold$CO2_threshold
)
df_cm_cvms <- df_cm_cvms[complete.cases(df_cm_cvms),]

# Calculate confusion matrix with comprehensive statistics
cm <- confusion_matrix(targets = df_cm_cvms$CO2, predictions = df_cm_cvms$AmmoLa)
cm_stats <- cm[1, c(
  "Balanced Accuracy", "F1", "Sensitivity", "Specificity",
  "Pos Pred Value", "Neg Pred Value", "Kappa", "MCC", "Detection Rate"
)]

# Format statistics for plot annotation
vec <- unlist(cm_stats)
formatted <- paste(names(vec), signif(vec, 3), sep = ": ")
midpoint <- ceiling(length(formatted) / 2)
stat_string <- paste(
  paste(formatted[1:midpoint], collapse = " | "),
  paste(formatted[(midpoint + 1):length(formatted)], collapse = " | "),
  sep = "\n"
)

# Fisher's exact test for the confusion matrix
tbl <- table(df_cm_cvms$AmmoLa, df_cm_cvms$CO2)
fisher_result <- fisher.test(tbl)
fisher_p <- signif(fisher_result$p.value, 3)

# Plot detailed confusion matrix
p_cm_AmmoLa_vs_CO2 <- plot_confusion_matrix(
  cm$`Confusion Matrix`[[1]],
  add_sums = FALSE,
  intensity_by = "counts"
) +
  ggplot2::labs(
    subtitle = stat_string,
    caption = paste("Fisher's exact test p-value:", fisher_p)
  ) +
  ggplot2::xlab("Predictions: AmmoLa") +
  ggplot2::ylab("Target: CO2") +
  theme_plot()

print(p_cm_AmmoLa_vs_CO2)
# ggsave("p_cm_AmmoLa_vs_CO2.svg", p_cm_AmmoLa_vs_CO2, width = 8, height = 8)

# --- Simplified 2x2 agreement tables for publication ---
# Format Fisher's test results for CO2 comparison
stat_text_Co2 <- paste0(
  "Fisher's exact test\n",
  "p = ", signif(ftest_AmmoLa_vs_CO2_most_sensitive$p.value, 3), "\n",
  "OR = ", signif(ftest_AmmoLa_vs_CO2_most_sensitive$estimate, 3)
)

# Create heatmap-style agreement table for AmmoLa vs CO2
p_table_AmmoLa_vs_CO2_most_sensitive <- ggplot(
  as.data.frame(table_AmmoLa_vs_CO2_most_sensitive),
  aes(x = CO2_threshold, y = AmmoLa_intensity, fill = Freq)
) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "cornsilk1", high = "cornsilk4") +
  theme_minimal() +
  labs(
    title = "Agreement of most sensitive subjects across tests",
    x = "CO2_threshold",
    y = "AmmoLa_intensity",
    fill = "Count"
  ) +
  annotate(
    "rect", xmin = 1.1, xmax = 1.9, ymin = 1.2, ymax = 1.8,
    alpha = 0.5, fill = "white", color = NA
  ) +
  annotate(
    "text", x = 1.5, y = 1.5,
    label = stat_text_Co2, size = 3, color = "black"
  ) +
  theme_plot() +
  coord_fixed(ratio = 1)

print(p_table_AmmoLa_vs_CO2_most_sensitive)

# Format Fisher's test results for Lateralization comparison
stat_text_lateralization <- paste0(
  "Fisher's exact test\n",
  "p = ", signif(ftest_AmmoLa_vs_lateralization_most_sensitive$p.value, 3), "\n",
  "OR = ", signif(ftest_AmmoLa_vs_lateralization_most_sensitive$estimate, 3)
)

# Create heatmap-style agreement table for AmmoLa vs Lateralization
p_table_AmmoLa_vs_lateralization_most_sensitive <- ggplot(
  as.data.frame(table_AmmoLa_vs_lateralization_most_sensitive),
  aes(x = Lateralization, y = AmmoLa_intensity, fill = Freq)
) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "cornsilk1", high = "cornsilk4") +
  theme_minimal() +
  labs(
    title = "Agreement of most sensitive subjects across tests",
    x = "Lateralization",
    y = "AmmoLa_intensity",
    fill = "Count"
  ) +
  annotate(
    "rect", xmin = 1.1, xmax = 1.9, ymin = 1.2, ymax = 1.8,
    alpha = 0.5, fill = "white", color = NA
  ) +
  annotate(
    "text", x = 1.5, y = 1.5,
    label = stat_text_lateralization, size = 3, color = "black"
  ) +
  theme_plot() +
  coord_fixed(ratio = 1)

print(p_table_AmmoLa_vs_lateralization_most_sensitive)

# ========================================================================== #
# 9. VARIABLE TRANSFORMATION AND ALTERNATIVE GMM BASED GROUPING ANALYSIS
# ========================================================================== #

print("Minimum value in AmmoLa_intensity")
print(min(trigeminal_measures_data$AmmoLa_intensity))

if (!analyze_only_untransformed) {
  cat("\nApplying variable transformations...\n")

  # Define transformation helper functions
  square <- function(x) x ^ 2

  # Apply transformations to address skewness
  trigeminal_measures_data <- trigeminal_measures_data %>%
    mutate(
      AmmoLa_intensity_reflect_slog = reflect_slog_unflipped(AmmoLa_intensity),
      Lateralization_square = square(Lateralization),
      CO2_threshold_log = slog(CO2_threshold) # replace regular log by slog for consistency
    )
}


# Assess univariate grouping in AmmoLA
max_modes <- 4

AmmoLa_GMM <- opGMMassessment::opGMMassessment(
  trigeminal_measures_data$AmmoLa_intensity_reflect_slog,
  MaxModes = max_modes,
  FitAlg = "DO",
  MaxCores = max_modes,
  PlotIt = TRUE,
  Seed = 42
)

# Create plots
# PDE plot for transformed AmmoLa intensity plus GMM boundary
pPDE_AmmoLA_GMM <- AmmoLa_GMM$Plot +
  theme_plot() +
  labs(title = "GMM of transformed AmmoLa intensity estimates", x = "Value", y = "PDE") +
  guides(color = "none")

# PDE plot for untransformed AmmoLa intensity plus GMM boundary
pPDE_AmmoLA_plus_GMM_boundary <- PDEplotGG(trigeminal_measures_data$AmmoLa_intensity) +
  theme_plot() +
  labs(title = "Distribution of AmmoLa intensity estimates", x = "Value", y = "PDE") +
  guides(color = "none") +
  geom_vline(xintercept = inv_reflect_slog_unflipped(y = tail(AmmoLa_GMM$Boundaries), original_max = max(trigeminal_measures_data$AmmoLa_intensity))) +
  geom_vline(xintercept = 90)

Lateralization_GMM <- opGMMassessment::opGMMassessment(
  trigeminal_measures_data$Lateralization_square,
  MaxModes = max_modes,
  FitAlg = "DO",
  MaxCores = max_modes,
  PlotIt = TRUE,
  Seed = 42
)

CO2_GMM <- opGMMassessment::opGMMassessment(
  trigeminal_measures_data$CO2_threshold_log[550:nrow(trigeminal_measures_data)],
  MaxModes = max_modes,
  FitAlg = "DO",
  MaxCores = max_modes,
  PlotIt = TRUE,
  Seed = 412
)

# Check if groups of AmmoLA have also higher values in Lateralization and CO2
# Upper 10% sensitivity
wilcox.test(trigeminal_measures_data$AmmoLa_intensity ~ trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity)
wilcox.test(trigeminal_measures_data$Lateralization ~ trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity)
wilcox.test(trigeminal_measures_data$CO2_threshold[550:nrow(trigeminal_measures_data)] ~ trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity[550:nrow(trigeminal_measures_data)])

cat("\nAnalyzing agreement among most sensitive subjects...\n")

trigeminal_measures_data_grouped_GMM <- trigeminal_measures_data[, 4:6]
trigeminal_measures_data_grouped_GMM$CO2_threshold_log[1:549] <- NA

# Create binary grouping: Upper mode = most sensitive vs others
trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog <- ifelse(
  trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog > tail(AmmoLa_GMM$Boundaries), 1, 0
)
trigeminal_measures_data_grouped_GMM$Lateralization_square <- ifelse(
  trigeminal_measures_data_grouped_GMM$Lateralization_square > tail(Lateralization_GMM$Boundaries), 1, 0
)
trigeminal_measures_data_grouped_GMM$CO2_threshold_log <- ifelse(
  trigeminal_measures_data_grouped_GMM$CO2_threshold <= tail(CO2_GMM$Boundaries), 1, 0
)

# Summary of sensitivity thresholding
cat("\nDistribution of most sensitive subjects:\n")
print(apply(trigeminal_measures_data_grouped_GMM[, 1:3], 2, table))

# Fisher's tests for breath hold group only
cat("\nFisher's test: AmmoLa vs Lateralization (breath hold only)\n")
print(fisher.test(table(
  trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog,
  trigeminal_measures_data_grouped_GMM$Lateralization_square
)))

cat("\nFisher's test: AmmoLa vs CO2 (breath hold only)\n")
print(fisher.test(table(
  trigeminal_measures_data_grouped_GMM$AmmoLa_intensity,
  trigeminal_measures_data_grouped_GMM$CO2_threshold
)))

# Upper sensitivity GMM mode
wilcox.test(trigeminal_measures_data$AmmoLa_intensity ~ trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog)
wilcox.test(trigeminal_measures_data$Lateralization ~ trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog)
wilcox.test(trigeminal_measures_data$CO2_threshold[550:nrow(trigeminal_measures_data)] ~ trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog[550:nrow(trigeminal_measures_data)])

# Prepare variables for clustering analysis
trigeminal_measures_data_CO2_breath_hold <- trigeminal_measures_data
trigeminal_measures_data_CO2_breath_hold$CO2_threshold[1:549] <- NA
trigeminal_measures_data_CO2_breath_hold$CO2_threshold_log[1:549] <- NA

vals <- unique(trigeminal_measures_data_CO2_breath_hold$AmmoLa_intensity)
vals_cut <- vals[!(vals %in% c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))]

# Perform brute force grouping analysis
brute_force_grouping <- lapply(vals_cut, function(AmmoLa) {
  group_AmmoLa <- ifelse(trigeminal_measures_data_CO2_breath_hold$AmmoLa_intensity < AmmoLa, 0, 1)

  # For Lateralization
  valid_lateralization <- !is.na(trigeminal_measures_data_CO2_breath_hold$Lateralization) & !is.na(group_AmmoLa)
  groups_lat <- unique(group_AmmoLa[valid_lateralization])
  if (length(groups_lat) == 2) {
    wilcoxon_lateralization <- wilcox.test(
      trigeminal_measures_data_CO2_breath_hold$Lateralization[valid_lateralization] ~
        group_AmmoLa[valid_lateralization]
    )
    stat_lat <- wilcoxon_lateralization$statistic
    p_lat <- wilcoxon_lateralization$p.value
  } else {
    stat_lat <- NA
    p_lat <- NA
  }

  # For CO2_threshold
  valid_CO2 <- !is.na(trigeminal_measures_data_CO2_breath_hold$CO2_threshold) & !is.na(group_AmmoLa)
  groups_CO2 <- unique(group_AmmoLa[valid_CO2])
  if (length(groups_CO2) == 2) {
    wilcoxon_CO2 <- wilcox.test(
      trigeminal_measures_data_CO2_breath_hold$CO2_threshold[valid_CO2] ~
        group_AmmoLa[valid_CO2]
    )
    stat_CO2 <- wilcoxon_CO2$statistic
    p_CO2 <- wilcoxon_CO2$p.value
  } else {
    stat_CO2 <- NA
    p_CO2 <- NA
  }

  list(
    wilcoxon_lateralization = stat_lat,
    wilcoxon_CO2 = stat_CO2,
    wilcoxon_lateralization_p = p_lat,
    wilcoxon_CO2_p = p_CO2
  )
})

results_list <- lapply(vals_cut, function(AmmoLa) {
  group_AmmoLa <- ifelse(trigeminal_measures_data_CO2_breath_hold$AmmoLa_intensity < AmmoLa, 0, 1)

  # For Lateralization
  valid_lateralization <- !is.na(trigeminal_measures_data_CO2_breath_hold$Lateralization) & !is.na(group_AmmoLa)
  groups_lat <- unique(group_AmmoLa[valid_lateralization])
  if (length(groups_lat) == 2) {
    wilcoxon_lateralization <- wilcox.test(
      trigeminal_measures_data_CO2_breath_hold$Lateralization[valid_lateralization] ~
        group_AmmoLa[valid_lateralization]
    )
    stat_lat <- wilcoxon_lateralization$statistic
    p_lat <- wilcoxon_lateralization$p.value
  } else {
    stat_lat <- NA
    p_lat <- NA
  }

  # For CO2_threshold
  valid_CO2 <- !is.na(trigeminal_measures_data_CO2_breath_hold$CO2_threshold) & !is.na(group_AmmoLa)
  groups_CO2 <- unique(group_AmmoLa[valid_CO2])
  if (length(groups_CO2) == 2) {
    wilcoxon_CO2 <- wilcox.test(
      trigeminal_measures_data_CO2_breath_hold$CO2_threshold[valid_CO2] ~
        group_AmmoLa[valid_CO2]
    )
    stat_CO2 <- wilcoxon_CO2$statistic
    p_CO2 <- wilcoxon_CO2$p.value
  } else {
    stat_CO2 <- NA
    p_CO2 <- NA
  }

  data.frame(
    AmmoLa = AmmoLa,
    wilcoxon_lateralization = stat_lat,
    wilcoxon_CO2 = stat_CO2,
    wilcoxon_lateralization_p = p_lat,
    wilcoxon_CO2_p = p_CO2
  )
})

# Combine all results into one data frame
brute_force_grouping_df <- do.call(rbind, results_list)
df <- brute_force_grouping_df

# Normalize columns for comparability
df_norm <- data.frame(
  wilcoxon_lateralization = df$wilcoxon_lateralization / max(df$wilcoxon_lateralization, na.rm = TRUE),
  wilcoxon_CO2 = df$wilcoxon_CO2 / max(df$wilcoxon_CO2, na.rm = TRUE),
  wilcoxon_lateralization_p = 1 - (df$wilcoxon_lateralization_p / max(df$wilcoxon_lateralization_p, na.rm = TRUE)),
  wilcoxon_CO2_p = 1 - (df$wilcoxon_CO2_p / max(df$wilcoxon_CO2_p, na.rm = TRUE))
)

# Compute Euclidean distance from the ideal point (1,1,1,1)
distances <- sqrt(rowSums((df_norm - 1) ^ 2))

# Find row(s) with minimum distance (closest to ideal)
best_row_index <- which.min(distances)
best_row <- df[best_row_index,]

best_row

# Normalize CO2 statistic and p-value
df_norm_CO2 <- data.frame(
  wilcoxon_CO2 = df$wilcoxon_CO2 / max(df$wilcoxon_CO2, na.rm = TRUE),
  wilcoxon_CO2_p = 1 - (df$wilcoxon_CO2_p / max(df$wilcoxon_CO2_p, na.rm = TRUE))
)

# Compute Euclidean distance to ideal point (1,1)
distances_CO2 <- sqrt(rowSums((df_norm_CO2 - 1) ^ 2))

# Get best row index and row
best_CO2_index <- which.min(distances_CO2)
best_CO2_row <- df[best_CO2_index,]

best_CO2_row

# ========================================================================== #
# 10. AGE AND SEX EFFECTS ANALYSIS
# ========================================================================== #

cat("\nAnalyzing age correlations...\n")

# Add age to dataset and calculate correlations
trigeminal_measures_data_age <- trigeminal_measures_data
trigeminal_measures_data_age$age <- as.numeric(trigeminale_data$Age)
corr_mat_age <- Hmisc::rcorr(
  as.matrix(trigeminal_measures_data_age),
  type = "pearson"
)

cat("\nAnalyzing sex differences...\n")

# Add sex to dataset
trigeminal_measures_data_sex <- trigeminal_measures_data
trigeminal_measures_data_sex$sex <- as.factor(trigeminale_data$Gender)

# Kruskal-Wallis tests for sex differences
sex_diff_trig <- apply(
  trigeminal_measures_data_sex[, 1:(ncol(trigeminal_measures_data_sex) - 1)],
  2,
  function(x) kruskal.test(x ~ as.factor(trigeminal_measures_data_sex$sex))
)

# Calculate effect sizes (eta-squared) for sex differences
variables <- colnames(trigeminal_measures_data_sex)[
  colnames(trigeminal_measures_data_sex) != "sex"
]

sex_stats <- lapply(variables, function(var) {
  x <- trigeminal_measures_data_sex[[var]]
  s <- trigeminal_measures_data_sex$sex
  idx <- complete.cases(x, s)
  x <- x[idx]
  s <- s[idx]

  # ANOVA for eta-squared calculation
  aovmod <- aov(x ~ s)
  ss_total <- sum((x - mean(x)) ^ 2)
  ss_between <- sum(tapply(x, s, function(g) length(g) * (mean(g) - mean(x)) ^ 2))
  eta2 <- ss_between / ss_total

  # Non-parametric test
  kruskal_p <- tryCatch(
    kruskal.test(x ~ s)$p.value,
    error = function(e) NA
  )

  n <- length(x)
  data.frame(variable = var, eta2 = eta2, kruskal_p = kruskal_p, n = n)
})

sex_stats_df <- bind_rows(sex_stats) %>%
  arrange(desc(eta2))

cat("\nSex effect sizes:\n")
print(sex_stats_df)

# ========================================================================== #
# 11. DISTRIBUTION PLOTS - HISTOGRAMS AND DENSITY ESTIMATES
# ========================================================================== #

cat("\nGenerating distribution plots...\n")

# Reshape data for plotting
trigeminal_measures_data_long <- reshape2::melt(trigeminal_measures_data)

# Filter to untransformed variables if requested
if (plot_only_untransformed) {
  trigeminal_measures_data_long <- trigeminal_measures_data_long[
    -grep("log|square", trigeminal_measures_data_long$variable),
  ]
}

# --- Non-CO2 measures ---
trigeminal_measures_data_long_1 <- trigeminal_measures_data_long %>%
  filter(!str_detect(variable, "CO2"))

p_distribution_non_CO2 <- ggplot(
  trigeminal_measures_data_long_1,
  aes(x = value)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = "cornsilk3", color = "white", alpha = 0.7
  ) +
  geom_density(color = "cornsilk4", linewidth = 1) +
  facet_wrap(~variable, scales = "free") +
  labs(
    title = "Distribution of non-CO2 related trigeminal measures",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_plot()

print(p_distribution_non_CO2)
# ggsave("p_distribution_tigeminal_nonCO2.svg",
#        p_distribution_tigeminal_nonCO2,
#        width = 10, height = ifelse(plot_only_untransformed, 5, 10))

# --- CO2 measures split by breath control protocol ---
trigeminal_measures_data_CO2 <- trigeminal_measures_data[
  , str_detect(names(trigeminal_measures_data), "CO2")
]

trigeminal_measures_data_CO2 <- trigeminal_measures_data_CO2 %>%
  mutate(
    row_id = row_number(),
    subset_group = ifelse(row_id <= 549, "Normal breath", "Breath hold")
  )

trigeminal_measures_data_CO2_long <- trigeminal_measures_data_CO2 %>%
  pivot_longer(
    cols = -c(row_id, subset_group),
    names_to = "variable",
    values_to = "value"
  )

if (plot_only_untransformed) {
  trigeminal_measures_data_CO2_long <- trigeminal_measures_data_CO2_long[
    -grep("log|square", trigeminal_measures_data_CO2_long$variable),
  ]
}

p_distribution_CO2 <- ggplot(
  trigeminal_measures_data_CO2_long,
  aes(x = value)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    fill = "cornsilk3", color = "white", alpha = 0.7
  ) +
  geom_density(color = "cornsilk4", size = 1) +
  facet_wrap(subset_group * variable ~ ., scales = "free", ncol = 2) +
  labs(
    title = "Distribution of CO2-related trigeminal measures, split by procedure",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_plot()

print(p_distribution_CO2)
# ggsave("p_distribution_tigeminal_CO2.svg",
#        p_distribution_tigeminal_CO2,
#        width = 10, height = ifelse(plot_only_untransformed, 5, 10))

# ========================================================================== #
# 12. PARETO DENSITY ESTIMATION (PDE) PLOTS
# ========================================================================== #

cat("\nGenerating Pareto Density Estimation plots...\n")

# PDE plot for AmmoLa intensity
pPDE_AmmoLA <- PDEplotGG(trigeminal_measures_data$AmmoLa_intensity) +
  theme_plot() +
  labs(title = "Distribution of AmmoLa intensity estimates", x = "Value", y = "PDE") +
  guides(color = "none") +
  scale_color_manual(values = "cornsilk4")

# PDE plot for Lateralization
pPDE_Lateralization <- PDEplotGG(trigeminal_measures_data$Lateralization) +
  theme_plot() +
  labs(title = "Distribution of the number of correct lateralizations", x = "Value", y = "PDE") +
  guides(color = "none") +
  scale_color_manual(values = "cornsilk4")

# PDE plot for CO2 thresholds (split by breath protocol)
CO2_threshold_breathing <- trigeminal_measures_data$CO2_threshold[1:549]
CO2_threshold_breath_hold <- trigeminal_measures_data$CO2_threshold[550:length(trigeminal_measures_data$CO2_threshold)]

# Pad vectors to equal length with NA
max_len <- max(length(CO2_threshold_breathing), length(CO2_threshold_breath_hold))
length(CO2_threshold_breathing) <- max_len
length(CO2_threshold_breath_hold) <- max_len

# Combine into data frame for PDE
CO2_df_for_PDE <- data.frame(
  CO2_threshold_breathing = CO2_threshold_breathing,
  CO2_threshold_breath_hold = CO2_threshold_breath_hold
)

pPDE_CO2 <- PDEplotGG(CO2_df_for_PDE) +
  theme_plot() +
  labs(title = "Distribution of CO2 thresholds", color = "Breath", x = "Value", y = "PDE") +
  scale_color_manual(
    values = c("grey83", "cornsilk4"),
    labels = c("uncontrolled", "hold")
  ) +
  theme(legend.position.inside = TRUE, legend.position = c(.2, .85))

# ========================================================================== #
# 13. CORRELATION MATRIX PREPARATION
# ========================================================================== #

cat("\nPreparing correlation matrix...\n")

# Handle transformed variables if enabled
if (!analyze_only_untransformed) {
  trigeminal_measures_data_transformed <- trigeminal_measures_data[, c("AmmoLa_intensity_reflect_slog", "Lateralization_square", "CO2_threshold_log")]

  # Invert CO2 so higher = more sensitive
  trigeminal_measures_data_transformed$CO2_threshold_log <- -trigeminal_measures_data_transformed$CO2_threshold_log

  # Remove breath-not-controlled CO2 data
  trigeminal_measures_data_transformed <- trigeminal_measures_data_transformed %>%
    mutate(CO2_threshold_log = ifelse(row_number() <= 549, NA, CO2_threshold_log))

  # Define variable pairs for correlation analysis
  var_pairs <- list(
    c("AmmoLa_intensity_reflect_slog", "Lateralization_square"),
    c("AmmoLa_intensity_reflect_slog", "CO2_threshold_log")
  )

  plot_list <- map(var_pairs, function(vars) {
    df <- trigeminal_measures_data_transformed[, vars]
    names(df) <- c("x", "y")
    df$pair <- paste(vars, collapse = " vs ")
    df
  })

  plot_df <- bind_rows(plot_list)
}

# Prepare comprehensive correlation dataset
all_measures_correlations <- trigeminal_measures_data

# Invert CO2 scales (lower threshold = higher sensitivity)
all_measures_correlations$CO2_threshold <- -all_measures_correlations$CO2_threshold
all_measures_correlations$CO2_threshold_log <- -all_measures_correlations$CO2_threshold_log

# Split CO2 by breath control protocol
all_measures_correlations$CO2_threshold_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold[1:549]

all_measures_correlations$CO2_threshold_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold[550:nrow(all_measures_correlations)]

all_measures_correlations$CO2_threshold_log_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_log_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold_log[1:549]

all_measures_correlations$CO2_threshold_log_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_log_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold_log[550:nrow(all_measures_correlations)]

# Filter variables based on analysis settings
if (analyze_only_untransformed | plot_only_untransformed) {
  all_measures_correlations <- all_measures_correlations[, - c(grep("log|square", names(all_measures_correlations)))]
}

# Remove redundant columns
all_measures_correlations <- all_measures_correlations[, !names(all_measures_correlations) %in% c("CO2_threshold_breath_not_hold", "CO2_threshold")]

# ========================================================================== #
# 14. CORRELATION MATRIX CALCULATION AND HEATMAP
# ========================================================================== #

cat("\nCalculating correlation matrix...\n")

# Select only numeric columns
numeric_data <- all_measures_correlations %>%
  dplyr::select(where(is.numeric))

# Calculate Spearman correlations with p-values
cor_method <- "spearman"
cor_results <- corr.test(numeric_data, use = "pairwise", method = cor_method)
print(cor_results)

corr_mat <- cor_results$r # Correlation coefficients
p_mat <- cor_results$p # P-values

#' Convert p-values to significance stars
signif_code <- function(p) {
  if (is.na(p)) ""
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

# Define color palette for correlation heatmap (NYT-inspired)
breaks <- seq(-1, 1, by = 0.1)
nyt_colors <- c("ghostwhite", "#f5f5dc", "#ede8d0",
                "cornsilk", "cornsilk2", "cornsilk3", "cornsilk4")
color_vec <- colorRampPalette(c(rev(nyt_colors), nyt_colors))(length(breaks))
col_fun <- colorRamp2(breaks, color_vec)

#' Determine text color based on background brightness
text_color_fun <- function(fill_color) {
  rgb_val <- col2rgb(fill_color) / 255
  brightness <- 0.299 * rgb_val[1,] + 0.587 * rgb_val[2,] + 0.114 * rgb_val[3,]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}

#' Create correlation heatmap with ComplexHeatmap
create_heatmap_trig <- function() {
  ht <- Heatmap(
    corr_mat,
    name = "Correlation",
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE,
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    cluster_columns = FALSE,
    clustering_method_columns = "ward.D2",
    show_column_dend = FALSE,
    column_names_rot = 45,
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(4, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- sprintf("%.2f", corr_mat[i, j])
      pval <- p_mat[i, j]
      star <- signif_code(pval)
      lbl <- paste0(val, star)
      col_txt <- text_color_fun(fill)
      grid.text(lbl, x, y, gp = gpar(fontsize = 10, col = col_txt))
    },
    rect_gp = gpar(col = NA),
    border = FALSE,
    row_names_side = "left",
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter",
      title_gp = gpar(fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    show_heatmap_legend = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )

  grid.newpage()
  draw(ht, heatmap_legend_side = "bottom", newpage = FALSE)
}

# Capture heatmap as graphical object
gp <- grid.grabExpr(create_heatmap_trig())

# Create separate title plot for consistent alignment
p_title <- ggplot() +
  labs(title = paste0("Correlation matrix (", cor_method, ")")) +
  theme_minimal(base_family = "Libre Franklin") +
  theme(
    plot.title = element_text(
      face = "plain", size = 12, color = "#222222",
      hjust = 0, margin = margin(b = 10)
    ),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 0, 20)
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "plain", size = 12, color = "#222222",
      hjust = 0, margin = margin(b = 5)
    ),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine title and heatmap
corr_plot <- plot_grid(
  p_title, gp,
  ncol = 1,
  rel_heights = c(0.05, 1),
  align = "v"
) + coord_fixed(ratio = 1)

print(corr_plot)
# ggsave("trigeminal_correlation_heatmap.svg", corr_plot, width = 8, height = 9)


# For comparison, do the violin plot for the AmmoLa_High_low tareg

# Clusters
trigeminal_ammola_high_low_data <- cbind.data.frame(Cluster = as.factor(trigeminal_measures_data_grouped$AmmoLa_intensity),
                                                    trigeminale_data[, names(trigeminale_data) %in%
                                                                       c(variables_by_categories$Nasal_chemosensory_perception,
                                                                         variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])])

head(trigeminal_ammola_high_low_data)
trigeminal_ammola_high_low_data$`CO2 threshold`[1:549] <- NA
trigeminal_ammola_high_low_data$`CO2 threshold` <- -trigeminal_ammola_high_low_data$`CO2 threshold`

trigeminal_ammola_high_low_data_long <- reshape2::melt(trigeminal_ammola_high_low_data, id.vars = "Cluster")
head(trigeminal_ammola_high_low_data_long)

round(apply(trigeminal_ammola_high_low_data[, -1], 2, function(x) wilcox.test(x ~ trigeminal_ammola_high_low_data$Cluster, na.rm = TRUE)$p.value), 2)


library(dplyr)
library(ggplot2)

# Dodge width for positioning (same as in jitterdodge and boxplot)
dodge_width <- 0.8

# Calculate the medians for each (Cluster, variable) combination
median_data <- trigeminal_ammola_high_low_data_long %>%
  group_by(variable, Cluster) %>%
  dplyr::summarise(median_value = mean(value, na.rm = T), .groups = "drop")
median_data$x_position <- rep(c(0.8, 1.2), nrow(median_data) / 2)

# Plot data and add a line for medians
set.seed(42)
p_trigeminal_ammola_high_low_data <- ggplot(trigeminal_ammola_high_low_data_long, aes(x = variable, y = value, color = Cluster, fill = Cluster)) +
  geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = dodge_width)) +
  geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = dodge_width), outlier.shape = NA) +
  geom_jitter(alpha = 1, size = 0.3,
              position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = dodge_width)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                             label.y.npc = "top",
                             vjust = -0.2, size = 3) +
# Add the median lines within dodge positions
geom_line(data = median_data, aes(x = x_position, y = median_value, group = variable),
            color = "black", size = 0.8, linetype = "solid", inherit.aes = FALSE) +
  facet_wrap(variable ~ ., nrow = 1, scales = "free", labeller = label_wrap_gen(width = 20)) +
  guides(color = "none") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk4", "cornsilk4")) +
  labs(title = "Raw trigeminal data per variable and AmmoLa sensitvity class", fill = "Cluster", x = NULL) +
  theme_plot() +
  theme(
    legend.direction = "horizontal",
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(colour = "black", size = 6, face = "plain"),
    axis.text.x = element_blank(), # Remove x axis tick labels
    axis.ticks.x = element_blank() # Remove x axis ticks (optional)
  ) +
  guides(color = "none")

p_trigeminal_ammola_high_low_data


library(effsize)
library(dplyr)
library(tidyr)
library(purrr)
library(boot)

# Define helper function to compute Cohen's d for two numeric vectors
cohen_d_val <- function(x, y) {
  cohen.d(x, y, hedges.correction = TRUE)$estimate
}

# Bootstrap function for Cohen's d
boot_cohen_d <- function(x, y, R = 1000) {
  data <- data.frame(group = rep(c("x", "y"), times = c(length(x), length(y))),
                     value = c(x, y))

  stat_fun <- function(data, indices) {
    d_x <- data$value[indices][data$group[indices] == "x"]
    d_y <- data$value[indices][data$group[indices] == "y"]
    if (length(d_x) < 2 || length(d_y) < 2) return(NA)
    return(cohen_d_val(d_x, d_y))
  }

  boot_out <- boot(data, stat_fun, R = R)

  ci <- boot.ci(boot_out, type = "perc")
  if (is.null(ci)) {
    return(c(NA, NA))
  } else {
    return(ci$percent[4:5])
  }
}

# List variables
variables <- unique(trigeminal_ammola_high_low_data_long$variable)

# Calculate Cohen's d and CIs per variable using map
results_cohen_AmmoLa <- map_df(variables, function(var) {
  data_var <- trigeminal_ammola_high_low_data_long %>%
    filter(variable == var) %>%
    filter(!is.na(value))

  x <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[1]]
  y <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[2]]

  # Skip variables with insufficient data
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(variable = var, d = NA, ci_lower = NA, ci_upper = NA))
  }

  d <- cohen_d_val(x, y)
  ci <- boot_cohen_d(x, y)

  tibble(variable = var, d = d, ci_lower = ci[1], ci_upper = ci[2])
})

# Plot Cohen's d with CI error bars
p_cohens_d_AmmoLa_gigh_low <- ggplot(results_cohen_AmmoLa, aes(x = reorder(variable, d), y = d)) +
  geom_col(fill = "cornsilk3") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  coord_flip() +
  labs(title = "Cohen's d (95% CI)",
       y = "Cohen's d",
       x = "Variables") +
  theme_plot() +
  geom_hline(yintercept = -0.8, linetype = "dashed", color = "salmon")

p_cohens_d_AmmoLa_gigh_low

# ========================================================================== #
# 15. COMBINED PUBLICATION FIGURE
# ========================================================================== #

cat("\nCreating combined publication figure...\n")

# Combine all main plots into publication-ready figure
combined_trigeminal_analysis_plot <- cowplot::plot_grid(
  cowplot::plot_grid(
# Top row: Distribution plots (A, B, C)
  cowplot::plot_grid(
    pPDE_CO2,
    pPDE_AmmoLA + geom_vline(xintercept = 90, linetype = "dashed", color = "salmon"),
    pPDE_Lateralization,
    labels = LETTERS[1:3],
    nrow = 1,
    align = "h",
    axis = "tb",
    label_y = 0.97
  ),
# Bottom row: Agreement tables (D, E) and correlation matrix (F)
  cowplot::plot_grid(
    cowplot::plot_grid(
      corr_plot,
      labels = LETTERS[4],
      nrow = 1,
      align = "h",
      axis = "tb",
      label_y = 0.97
    ),
    cowplot::plot_grid(
      p_table_AmmoLa_vs_CO2_most_sensitive,
      p_table_AmmoLa_vs_lateralization_most_sensitive,
      labels = LETTERS[5:6],
      nrow = 1,
      align = "h",
      axis = "tb",
      label_y = 0.97
    ),
    nrow = 1,
    rel_widths = c(1, 2),
    align = "h",
    axis = "tb"
  ),
cowplot::plot_grid(
  p_trigeminal_ammola_high_low_data,
  nrow = 1,
  align = "h",
  axis = "tb",
  label_y = 0.97,
  labels = LETTERS[7]
),
  nrow = 3,
  rel_heights = c(2, 2, 3),
  align = "v",
  axis = "lr"
),
cowplot::plot_grid(
  p_cohens_d_AmmoLa_gigh_low,
  labels = LETTERS[8]),
ncol = 2,
rel_widths = c(5, 3),
align = "h",
axis = "tb"
)

# Display combined figure
print(combined_trigeminal_analysis_plot)

# Save combined figure
ggsave(
  "combined_trigeminal_analysis_plot.svg",
  combined_trigeminal_analysis_plot,
  width = 40,
  height = 18
)


# ========================================================================== #
# 16. CLUSTERING
# ========================================================================== #

# Extract relevant variables for clustering analysis
variables_nasal_chemosensory_perception <- trigeminale_data[c("ID", variables_by_categories$Nasal_chemosensory_perception)]

# Check for missing values in the dataset
cat("Checking which variables or cases have missing values above acceptable thresholds...\n")
n_missings_per_variable <- colSums(is.na(variables_nasal_chemosensory_perception))

# Define threshold for acceptable missing values
f_missings <- 0.2
f_missings_per_variable <- n_missings_per_variable / nrow(variables_nasal_chemosensory_perception)
to_drop <- names(variables_nasal_chemosensory_perception)[which(f_missings_per_variable > f_missings)]

# Filter out variables with excessive missing values
variables_for_clustering <- variables_nasal_chemosensory_perception[, !names(variables_nasal_chemosensory_perception) %in% c("ID", to_drop)]
n_missings_per_variable_for_clustering <- colSums(is.na(variables_for_clustering))
cat("n_missings_per_variable_for_clustering\n")
print(n_missings_per_variable_for_clustering)
cat("Sum of missings: ")
print(sum(n_missings_per_variable_for_clustering))
print(sum(n_missings_per_variable_for_clustering) / length(unlist(as.vector(variables_for_clustering))))


# Impute missing values using Random Forest
set.seed(42)
Impu <- try(suppressWarnings(missForest::missForest(variables_for_clustering, maxiter = 10000)), TRUE)

# Check for successful imputation
if (!inherits(Impu, "try-error")) {
  ImputedData <- round(Impu$ximp)
}

# Report any remaining missing values post-imputation
n_missings_per_variable_imputed <- colSums(is.na(ImputedData))
print(n_missings_per_variable_imputed)

# Prepare the imputed dataset for clustering analysis
variables_for_clustering_imputed <- ImputedData
rownames(variables_for_clustering_imputed) <- variables_nasal_chemosensory_perception$ID

# # ========================================================================== #
# # CLUSTER NUMBER DETECTION
# # ========================================================================== #
#
# # Detect optimal number of clusters using NbClust
# set.seed(42)
# ClusterIndices <- NbClust::NbClust(data = variables_for_clustering_imputed,
#                                    min.nc = 2,
#                                    max.nc = 5,
#                                    method = "centroid",
#                                    index = "all")
#
# # Display the best number of clusters
# cat("Best number of clusters detected:\n")
# table(ClusterIndices$Best.nc[1,])
# nCluster <- as.numeric(names(which.max(table(ClusterIndices$Best.nc[1,]))))
# par(mfrow = c(1, 1))  # Reset plotting parameters
#
# # ========================================================================== #
# # CLUSTERING ANALYSIS
# # ========================================================================== #
#
# # Perform hierarchical clustering using Ward's method
# set.seed(42)
# clusterObject <- stats::hclust(dist(variables_for_clustering_imputed), method = "ward.D2")
#
# # Visualize the dendrogram for hierarchical clustering
# dendro <- as.dendrogram(clusterObject)
# plot(dendro)
#
# # Cut the dendrogram to form clusters
# clusters <- cutree(clusterObject, k = nCluster)
# cat("Cluster assignments:\n")
# print(table(clusters))
#
# # ========================================================================== #
# # PRINCIPAL COMPONENT ANALYSIS (PCA)
# # ========================================================================== #
#
# # Perform PCA on the imputed clustering data
# res.pca <- FactoMineR::PCA(variables_for_clustering_imputed, scale.unit = TRUE, graph = TRUE, ncp = 8)
#
# # Display eigenvalues and variable contributions
# cat("PCA Eigenvalues:\n")
# print(res.pca$eig)
# cat("PCA Variable Contributions:\n")
# print(res.pca$var)
#
# # Perform HCPC for hierarchical clustering on PCA results
# res.hcpc <- HCPC(res.pca, graph = FALSE, nb.clust = 2)
#
# # ========================================================================== #
# # PCA SCREE PLOT
# # ========================================================================== #
# pScree <- fviz_screeplot(res.pca, addlabels = TRUE,
#                          barfill = "cornsilk1", barcolor = "cornsilk3", baralpha = 0.7) +
#   theme_light() +
#   labs(title = "PCA Scree Plot") +
#   ylim(0, 1.1 * max(res.pca$eig[, 2]))
#
#
# # ========================================================================== #
# # PCA VARIABLE CONTRIBUTION PLOTS
# # ========================================================================== #
#
# # Top 20 variable contributions to first principal component
# pContrib <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20,
#                          fill = "cornsilk1", color = "cornsilk3",
#                          ncp = sum(res.pca$eig[, 1] > 1)) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = paste0("PCA Variable Contribution to Relevant PCs"), x = NULL)
#
# # PCA variable contribution to PC1
# pContrib1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20,
#                           fill = "cornsilk1", color = "cornsilk3", ncp = 1) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "PCA Variable Contribution to PC1", x = NULL)
#
#
# # PCA overal variable contributions
# pPCA_biplot <- fviz_pca_var(res.pca, col.var="contrib",
#                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                             repel = TRUE )
#
#
# # ========================================================================== #
# # PCA FACTOR PLOT
# # ========================================================================== #
# pFactorplot <- fviz_ellipses(res.pca, habillage = factor(res.hcpc$data.clust$clust)) +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "PCA Factor Plot")
#
# # ========================================================================== #
# # DENDROGRAM AND CLUSTER VISUALIZATION
# # ========================================================================== #
#
# # Visualize dendrogram with clusters
# fviz_dend(res.hcpc,
#           cex = 0.7,                     # Label size
#           palette = "jco",               # Color palette
#           rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#           rect_border = "jco",           # Rectangle color
#           labels_track_height = 0.8      # Augment the room for labels
# )
#
# # Visualize clusters on the PCA plot
# fviz_cluster(res.hcpc,
#              repel = TRUE,            # Avoid label overlapping
#              show.clust.cent = TRUE, # Show cluster centers
#              palette = "jco",         # Color palette
#              ggtheme = theme_minimal(),
#              main = "Factor Map"
# )
#
# # ========================================================================== #
# # U-MATRIX FOR CLUSTER VISUALIZATION
# # ========================================================================== #
# Umatrix::iEsomTrain(Data = as.matrix(variables_for_clustering_imputed),
#                     Cls = clusters, Toroid = TRUE)

# ========================================================================== #
# END OF STANDARD CLUSTERING SECTION
# ========================================================================== #


# ========================================================================== #
# Second Unsupervised Projection and Clustering Analysis
# ========================================================================== #

# ========================================================================== #
# 1. LOAD REQUIRED FUNCTIONS AND SOURCES
# ========================================================================== #

# Load main projection and plotting functions
# source("/home/joern/.Datenplatte/Joerns Dateien/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/ProjectionsBiomed_MainFunctions_6_1core.R")
source("ProjectionsBiomed_MainFunctions_6_1core.R")

# Load ABC plotting utility
source("/home/joern/Aktuell/ABCplotGG.R")

# ========================================================================== #
# 2. SET ANALYSIS PARAMETERS
# ========================================================================== #

# Define computational resources
nProc_possible <- parallel::detectCores() - 1
nProc_desired <- 48

# Define dimensionality reduction (projection) methods
projection_methods <- c("none", "PCA", "ICA", "MDS", "tSNE", "Umap")

# Define clustering algorithms
clustering_methods <- c("kmeans", "kmedoids",
                        "ward.D2", "single", "average",
                        "median", "complete", "centroid")

# Define method for determining number of clusters
cluster_number_methods <- "NbClust"

# Define clustering evaluation indices
unsupervised_metrics_full <- c("Silhouette_index", "Dunn_index",
                               "DaviesBouldin_index", "dbcv_index",
                               "CalinskiHarabasz_index", "inertia")

unsupervised_metrics_reduced <- c("Silhouette_index", "Dunn_index",
                                  "CalinskiHarabasz_index")

# Define validation and comparison metrics
valid_cluster_metrics <- c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                           "Rand_index", "DaviesBouldin_index", "dbcv_index",
                           "CalinskiHarabasz_index", "inertia",
                           "adjusted_mutual_information")

# ========================================================================== #
# 3. LOAD AND PREPARE DATA
# ========================================================================== #

# Backup original variable names
original_names <- names(variables_nasal_chemosensory_perception_for_clustering_imputed)

# Assign generic variable names for internal use
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed <-
  variables_nasal_chemosensory_perception_for_clustering_imputed
names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed) <-
  paste0("Var", seq_len(ncol(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed)))

#Add AmmoLa
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed$AmmoLa <-
  scaleRange_01(trigeminal_measures_data$AmmoLa_intensity_reflect_slog) * 3

# Prepare dataset for analysis
variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared <-
  prepare_dataset(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed)

# Inspect dataset structure
names(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared)
table(variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Target)

# Distribution analysis
# Inspection of the pdf
heat_matrix <- cbind.data.frame(variables_nasal_chemosensory_perception_for_clustering_imputed,
                                AmmoLa = variables_nasal_chemosensory_perception_for_clustering_imputed_renamed$AmmoLa)
row_means <- rowMeans(heat_matrix)
p_PDE_chemosensory_perception_for_clustering_imputed <- PDEplotGG(row_means)
print(p_PDE_chemosensory_perception_for_clustering_imputed)

GMM_chemosensory_perception_for_clustering_imputed <-
  opGMMassessment::opGMMassessment(row_means, MaxModes = 4, FitAlg = "DO", MaxCores = 4, PlotIt = TRUE, Seed = 42)

DescTools::AndersonDarlingTest(row_means, "pnorm", mean = mean(row_means), sd = sd(row_means))

# Clustering
DatasetNames <- c("variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared")

# Ensure data object is available globally
# assign("variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared",
#        variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared,
#        envir = .GlobalEnv)

# ========================================================================== #
# 4. PERFORM PROJECTION AND CLUSTERING ANALYSIS
# ========================================================================== #

set.seed(42)

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
  max_clusters = 5,
)

# Optionally save analysis results
# saveRDS(projectionsAndPlots_TriFunQ, file = "projectionsAndPlots_TriFunQ.RData")
# projectionsAndPlots_TriFunQ <- readRDS("projectionsAndPlots_TriFunQ.RData")

# ========================================================================== #
# 5. COMBINE AND SAVE PROJECTION PLOTS
# ========================================================================== #

combined_plots <- combine_all_plots(
  datasets = DatasetNames,
  projection_plots = projectionsAndPlots_TriFunQ$projections_plots,
  projection_methods = projection_methods,
  clustering_methods = clustering_methods,
  cluster_number_methods = cluster_number_methods
)

print(combined_plots)

ggsave(
  filename = paste0("Combined_projection_and_clustering_analysis_plot_", cluster_number_methods, "_clusters_",
                    DatasetNames[1], ".svg"),
  plot = combined_plots$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed,
  width = 3.5 * (length(clustering_methods) + length(cluster_number_methods)),
  height = 3.5 * length(projection_methods),
  limitsize = FALSE
)

# ========================================================================== #
# 6. CLUSTER QUALITY EVALUATION
# ========================================================================== #

dfClusterQuality <- projectionsAndPlots_TriFunQ$cluster_quality_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed

# Exclude metrics not used in current visualization
excluded_metrics <- setdiff(valid_cluster_metrics, unsupervised_metrics_full)
dfClusterQuality <- dfClusterQuality[, !names(dfClusterQuality) %in%
                                       c(excluded_metrics,
                                         paste0(excluded_metrics, "_rank"),
                                         "combined_rank_metrics_with_orig_classes",
                                         "combined_rank")]

# Generate unique method identifiers
rownames(dfClusterQuality) <- apply(dfClusterQuality[, 1:3], 1, paste0, collapse = "_")
dfClusterQuality$Method <- row.names(dfClusterQuality)

# Sort by best overall performance
dfClusterQuality <- dfClusterQuality[order(-dfClusterQuality$combined_rank_metrics_without_orig_classes),]

# Heatmap overview of clustering ranks
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
  labs(title = "Ranked cluster quality scores", fill = "Scaled\nrank", y = "Projection_clustering_clusternumber", x = NULL)

print(clustering_heatmap)


# Bar plot of Calinski-Harabasz index across methods
ggplot(dfClusterQuality, aes(y = CalinskiHarabasz_index, x = reorder(Method, - CalinskiHarabasz_index))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

write.csv(round(dfClusterQuality, 3), "dfClusterQuality.csv")

# ========================================================================== #
# 7. Assign CLUSTER LABEL TO ONE SAMPLE WHICH WAS OMITTED DUE TO PROJETCION PROBLEMS
# ========================================================================== #
TriFunQ_clusters <-
  projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Umap$centroid$NbClust$Cluster
length(TriFunQ_clusters)
k_clusters <- length(unique(TriFunQ_clusters))


# omitted_sample <- variables_nasal_chemosensory_perception_for_clustering_imputed_renamed[479, , drop = FALSE]
#
# # Assuming cluster_labels is your vector of cluster labels with NA inserted at position 479
# TriFunQ_clusters <- append(TriFunQ_clusters, NA, after = 478)
# length(TriFunQ_clusters)
#
# # Calculate distances from omitted sample to all clustered samples
# clustered_samples <- variables_nasal_chemosensory_perception_for_clustering_imputed_renamed[-479, , drop = FALSE]
# valid_labels <- TriFunQ_clusters[!is.na(TriFunQ_clusters)]
# dists <- apply(clustered_samples, 1, function(x) sqrt(sum((x - omitted_sample)^2)))
#
# # Find the index of the nearest neighbor
# nn_index <- which.min(dists)
#
# # Assign the cluster label of the nearest neighbor to omitted sample
# TriFunQ_clusters[479] <- valid_labels[nn_index]
# length(TriFunQ_clusters)

# ========================================================================== #
# 7. CHECK AGREEMENT WITH AmmoLA Grouping
# ========================================================================== #

# --- Simplified 2x2 agreement tables for publication ---

# Fisher's exact test for the confusion matrix
AmmoLa_intensity <- trigeminal_measures_data_grouped_breath_hold$AmmoLa_intensity + 1
table_AmmoLa_vs_TriFunQ_clusters <- table(AmmoLa_intensity, TriFunQ_clusters)
ftest_AmmoLa_vs_TriFunQ_clusters <- fisher.test(table_AmmoLa_vs_TriFunQ_clusters)
cat("\nAmmoLa vs TriFunQ clusters:\n")
print(ftest_AmmoLa_vs_TriFunQ_clusters)


# Format Fisher's test results for CO2 comparison
stat_text_TriFunQ_clusters <- paste0(
  "Fisher's exact test\n",
  "p = ", signif(ftest_AmmoLa_vs_TriFunQ_clusters$p.value, 3), "\n",
  "OR = ", signif(ftest_AmmoLa_vs_TriFunQ_clusters$estimate, 3)
)

# Create heatmap-style agreement table for AmmoLa vs CO2
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
  theme_plot() +
  coord_fixed(ratio = 1)

print(p_table_AmmoLa_vs_TriFunQ_clusters)

# ========================================================================== #
# 7. ADD CLUSTER DIAGNOSTICS PLOTS
# ========================================================================== #

# Silhouette plot
# Compute silhouette
diss_matrix <- dist(projectionsAndPlots_TriFunQ$projection_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Umap$Projected)
sil <- silhouette(projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Umap$centroid$NbClust$Cluster, diss_matrix)
sil_df <- as.data.frame(sil)

# Add a case/order column and sort by cluster, then silhouette width
sil_df$case <- as.numeric(rownames(sil_df))
sil_df <- sil_df %>%
  arrange(cluster, - sil_width) %>%
  mutate(id = row_number())

# Plot with ggplot2
p_TriFunQ_clusters_silhouette <- ggplot(sil_df, aes(x = id, y = sil_width, fill = as.factor(cluster))) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
  labs(title = "Silhouette plot by case", x = "Case", y = "Silhouette width", fill = "Cluster") +
  theme_plot() + scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position.inside = TRUE, legend.position = c(.2, .2), legend.direction = "horizontal")

p_TriFunQ_clusters_silhouette

# Dedrogram
# Source - https://stackoverflow.com/a
# Posted by andschar, modified by community. See post 'Timeline' for change history
# Retrieved 2025-11-14, License - CC BY-SA 4.0

hc <- hclust(diss_matrix, method = "ward.D2")

require(magrittr)
require(ggplot2)
require(dendextend)

dend <- hc %>% as.dendrogram %>%
  set("branches_k_color", value = c("cornsilk3", "cornsilk4"), k = k_clusters) %>% set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>% set("labels_colors", k = 2) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>%
  set("labels_colors", value = c("cornsilk3", "cornsilk4"), k = k_clusters)

ggd1 <- as.ggdend(dend)

p_TriFunQ_clusters_dend <- ggplot(ggd1, horiz = FALSE) +
  labs(title = "Dendrogram", x = "Case", y = "Dissimilarity") +
  theme_plot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position.inside = TRUE, legend.position = c(.2, .2))

p_TriFunQ_clusters_dend

# Projectcion and cluster plot

df_projected_and_clustered <- cbind.data.frame(
  Target = projectionsAndPlots_TriFunQ$clustering_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Umap$centroid$NbClust$Cluster,
  projectionsAndPlots_TriFunQ$projection_results$variables_nasal_chemosensory_perception_for_clustering_imputed_renamed_prepared$Umap$Projected
)

# Plot with points for each class and 95% confidence ellipses
p_TriFunQ_clusters_cluster <- ggplot(df_projected_and_clustered, aes(x = Dim1, y = Dim2, shape = as.factor(Target), fill = as.factor(Target), color = as.factor(Target))) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(level = 0.95, aes(fill = as.factor(Target)), alpha = 0.3, geom = "polygon") +
  theme_plot() +
  labs(title = "Projection scatter plot with confidence ellipses",
       x = "Dim 1", y = "Dim 2", shape = "Cluster") +
  scale_color_manual(values = c("grey33", "grey11")) +
  scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  theme(legend.position.inside = TRUE, legend.position = c(.1, .2)) +
  guides(color = "none", fill = "none")


p_TriFunQ_clusters_cluster

# ========================================================================== #
# 8. PLOT CLUSTER VARIABLES BY CLUSTER
# ========================================================================== #

# Clusters
trigeminal_clustered_data <- cbind.data.frame(Cluster = as.factor(TriFunQ_clusters),
                                   trigeminale_data[, names(trigeminale_data) %in%
                                                              c(variables_by_categories$Nasal_chemosensory_perception,
                                                                variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])])

head(trigeminal_clustered_data)
trigeminal_clustered_data$`CO2 threshold`[1:549] <- NA
trigeminal_clustered_data$`CO2 threshold` <- -trigeminal_clustered_data$`CO2 threshold`

trigeminal_clustered_data_long <- reshape2::melt(trigeminal_clustered_data, id.vars = "Cluster")
head(trigeminal_clustered_data_long)

apply(trigeminal_clustered_data[, -1], 2, function(x) wilcox.test(x ~ trigeminal_clustered_data$Cluster, na.rm = TRUE)$p.value)

library(dplyr)
library(ggplot2)

# Dodge width for positioning (same as in jitterdodge and boxplot)
dodge_width <- 0.8

# Calculate the medians for each (Cluster, variable) combination
median_data <- trigeminal_clustered_data_long %>%
  group_by(variable, Cluster) %>%
  dplyr::summarise(median_value = mean(value, na.rm = T), .groups = "drop")
median_data$x_position <- rep(c(0.8, 1.2), nrow(median_data) / 2)

# Plot data and add a line for medians
set.seed(42)
p_trigeminal_clustered_data <- ggplot(trigeminal_clustered_data_long, aes(x = variable, y = value, color = Cluster, fill = Cluster)) +
  geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = dodge_width)) +
  geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = dodge_width), outlier.shape = NA) +
  geom_jitter(alpha = 1, size = 0.3,
              position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = dodge_width)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                             label.y.npc = "top",
                             vjust = -0.2, size = 3) +
# Add the median lines within dodge positions
geom_line(data = median_data, aes(x = x_position, y = median_value, group = variable),
            color = "black", size = 0.8, linetype = "solid", inherit.aes = FALSE) +
  facet_wrap(variable ~ ., nrow = 1, scales = "free", labeller = label_wrap_gen(width = 20)) +
  guides(color = "none") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk4", "cornsilk4")) +
  labs(title = "Raw trigeminal data per variable and cluster", fill = "Cluster", x = NULL) +
  theme_plot() +
  theme(
    legend.direction = "horizontal",
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(colour = "black", size = 6, face = "plain"),
    axis.text.x = element_blank(), # Remove x axis tick labels
    axis.ticks.x = element_blank() # Remove x axis ticks (optional)
  ) +
  guides(color = "none")

p_trigeminal_clustered_data


library(effsize)
library(dplyr)
library(tidyr)
library(purrr)
library(boot)

# Define helper function to compute Cohen's d for two numeric vectors
cohen_d_val <- function(x, y) {
  cohen.d(x, y, hedges.correction = TRUE)$estimate
}

# Bootstrap function for Cohen's d
boot_cohen_d <- function(x, y, R = 1000) {
  data <- data.frame(group = rep(c("x", "y"), times = c(length(x), length(y))),
                     value = c(x, y))

  stat_fun <- function(data, indices) {
    d_x <- data$value[indices][data$group[indices] == "x"]
    d_y <- data$value[indices][data$group[indices] == "y"]
    if (length(d_x) < 2 || length(d_y) < 2) return(NA)
    return(cohen_d_val(d_x, d_y))
  }

  boot_out <- boot(data, stat_fun, R = R)

  ci <- boot.ci(boot_out, type = "perc")
  if (is.null(ci)) {
    return(c(NA, NA))
  } else {
    return(ci$percent[4:5])
  }
}

# List variables
variables <- unique(trigeminal_clustered_data_long$variable)

# Calculate Cohen's d and CIs per variable using map
results_cohen_TriFunQ <- map_df(variables, function(var) {
  data_var <- trigeminal_clustered_data_long %>%
    filter(variable == var) %>%
    filter(!is.na(value))

  x <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[1]]
  y <- data_var$value[data_var$Cluster == levels(data_var$Cluster)[2]]

  # Skip variables with insufficient data
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(variable = var, d = NA, ci_lower = NA, ci_upper = NA))
  }

  d <- cohen_d_val(x, y)
  ci <- boot_cohen_d(x, y)

  tibble(variable = var, d = d, ci_lower = ci[1], ci_upper = ci[2])
})

# Plot Cohen's d with CI error bars
p_cohens_d <- ggplot(results_cohen_TriFunQ, aes(x = reorder(variable, d), y = d)) +
  geom_col(fill = "cornsilk3") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  coord_flip() +
  labs(title = "Cohen's d (95% CI)",
       y = "Cohen's d",
       x = "Variables") +
  theme_plot() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "salmon")

p_cohens_d

# ========================================================================== #
# 9. PLOT CLUSTER VARIABLES HEATMAP BY CLUSTER
# ========================================================================== #

row_order <- order(row_means)
# Reorder the heat matrix rows
heat_matrix_sorted <- heat_matrix[row_order,]

# Reorder the annotation accordingly
pal2 <- colorRampPalette(c("cornsilk", "cornsilk3", "cornsilk4", "grey33"))

row_means_ha <-
  rowAnnotation(bar = anno_barplot(row_means[row_order]),
                show_legend = FALSE,
                annotation_label = "Row means",
                annotation_name_rot = 90,
                show_annotation_name = TRUE,
                width = unit(2, "cm")
  )

#' Create correlation heatmap with ComplexHeatmap
create_heatmap_trig_clustered_data <- function() {
  Heatmap_trigeminal_clustered_data <-
  as.ggplot(
  grid.grabExpr(draw(
  ComplexHeatmap::Heatmap(
  as.matrix(heat_matrix_sorted),
  right_annotation = row_means_ha,
  col = pal2(112), #rev(heat.colors(100)),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  heatmap_legend_param = list(title = "Sensitivity [0,3]", direction = "horizontal", title_position = "lefttop"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_max_height = unit(16, "cm"),
  row_split = -TriFunQ_clusters[row_order]),
  heatmap_legend_side = "bottom")

  )
  )

  grid.newpage()
  draw(ht, heatmap_legend_side = "bottom", newpage = FALSE)
}

# Capture heatmap as graphical object
gp_clustered_variables <- grid.grabExpr(create_heatmap_trig_clustered_data())

# Create separate title plot for consistent alignment
p_title_heat <- ggplot() +
  labs(title = "Clustered variables") +
  theme_minimal(base_family = "Libre Franklin") +
  theme(
    plot.title = element_text(
      face = "plain", size = 12, color = "#222222",
      hjust = 0, margin = margin(b = 10)
    ),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 0, 20)
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "plain", size = 12, color = "#222222",
      hjust = 0, margin = margin(b = 5)
    ),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine title and heatmap
heat_plot_clustered_variables <- plot_grid(
  p_title_heat,
  gp_clustered_variables,
  ncol = 1,
  rel_heights = c(0.05, 1),
  align = "v"
) + coord_fixed(ratio = 1)

print(heat_plot_clustered_variables)
# ggsave("trigeminal_correlation_heatmap_clustered_variables.svg", corr_plot_clustered_variables, width = 8, height = 9)

# ========================================================================== #
# 11. MAKE PCA PLOTS
# ========================================================================== #
# Perform PCA on the imputed clustering data
res.pca_trigeminal_clustered_data <- FactoMineR::PCA(heat_matrix, scale.unit = TRUE, graph = TRUE, ncp = Inf)
res.pca_trigeminal_clustered_data_biplot <- stats::princomp(scale(heat_matrix, center = TRUE, scale = TRUE)) #, scale.unit = TRUE, graph = TRUE, ncp = 8)

# Extract variable coordinates (variables x components)
var_coords <- res.pca_trigeminal_clustered_data$var$coord

# Create a summary table with eigenvalues, standard deviations, variance explained, cumulative variance
eig <- res.pca_trigeminal_clustered_data$eig
std_dev <- sqrt(eig[, 1])
prop_var <- eig[, 2] / 100
cum_var <- eig[, 3] / 100

summary_table <- round(data.frame(
  Eigenvalue = eig[, 1],
  Std_Dev = std_dev,
  Prop_Variance = prop_var,
  Cum_Variance = cum_var
), 4)

rownames(summary_table) <- paste0("Dim.", seq_len(nrow(eig)))

coord_table <- round(res.pca_trigeminal_clustered_data$var$coord, 4)

# Combine variable coordinates and empty rows for summary table variables
combined_table <- rbind(coord_table, t(summary_table))

# Write to CSV
print(combined_table)
write.csv(combined_table[, 1:4], "res.pca_trigeminal_clustered_data_loadings.csv")

# Display eigenvalues and variable contributions
cat("PCA Eigenvalues:\n")
print(res.pca_trigeminal_clustered_data$eig)
cat("PCA Variable Contributions:\n")
print(res.pca_trigeminal_clustered_data$var)

# ========================================================================== #
# PCA SCREE PLOT
# ========================================================================== #

pScree <-
  factoextra::fviz_screeplot(res.pca_trigeminal_clustered_data,
                             fill = "name",
                                     ncp = 100,
                                     choice = "variance",
                                     addlabels = T,
                         barfill = c(rep("cornsilk", sum(eig[, 1] > 1)), rep("ghostwhite", sum(eig[, 1] <= 1))), barcolor = "cornsilk3", baralpha = 0.7) +
  theme_plot() +
  labs(title = "PCA Scree Plot")
#ylim(0, 1.1 * max(res.pca$eig[, 2]))

# eig.val <- get_eig(res.pca_trigeminal_clustered_data)
#
# # Customize labels with more digits
# pScree <- pScree + geom_text(aes(label = sprintf("%.3f", eig.val[, "eigenvalue"])), angle = 90, hjust = -0.2, vjust = 0)


# ========================================================================== #
# PCA FACTOR PLOT
# ========================================================================== #

pBiplot <-
  factoextra::fviz_pca_var(res.pca_trigeminal_clustered_data_biplot,
                           col.var = "contrib", gradient.cols = c("cornsilk3", "cornsilk4", "grey33"), #pal2(200)[150:200],
                           repel = TRUE) +
  theme_plot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position.inside = TRUE, legend.position = c(.2, .2), fill = "Contribution") +
  labs(title = paste0("PCA biplot"))

# ========================================================================== #
# PCA VARIABLE CONTRIBUTION PLOTS
# ========================================================================== #

pContrib_dim1 <-
  factoextra::fviz_contrib(res.pca_trigeminal_clustered_data, choice = "var", axes = 1,
                         fill = "cornsilk1", color = "cornsilk3",
                         ncp = sum(res.pca$eig[, 1] > 1)) +
  theme_plot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = paste0("PCA variable contribution to PC1"), x = NULL)

pContrib_dim2 <-
  factoextra::fviz_contrib(res.pca_trigeminal_clustered_data, choice = "var", axes = 2,
                           fill = "cornsilk1", color = "cornsilk3",
                           ncp = sum(res.pca$eig[, 1] > 1)) +
  theme_plot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = paste0("PCA variable contribution to PC2"), x = NULL)


# ========================================================================== #
# CORRELATION MATRIX CALCULATION AND HEATMAP FOR CLUSTERED VARIABLES
# ========================================================================== #

cat("\nCalculating correlation matrix for clustered variables...\n")


# Calculate Spearman correlations with p-values
cor_results_clustered_variables <- corr.test(heat_matrix, method = cor_method)
print(cor_results_clustered_variables)

corr_mat_clustered_variables <- cor_results_clustered_variables$r # Correlation coefficients
p_mat_clustered_variables <- cor_results_clustered_variables$p # P-values


#' Create correlation heatmap with ComplexHeatmap
create_heatmap_corr_trig_clustered_variables <- function() {
  ht <- Heatmap(
    corr_mat_clustered_variables,
    name = "Correlation",
    col = col_fun,
    na_col = "white",
    cluster_rows = TRUE,
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    cluster_columns = TRUE,
    clustering_method_columns = "ward.D2",
    show_column_dend = FALSE,
    column_names_rot = 90,
    row_dend_width = unit(5, "cm"),
    column_dend_height = unit(6, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- sprintf("%.2f", corr_mat_clustered_variables[i, j])
      pval <- p_mat_clustered_variables[i, j]
      star <- signif_code(pval)
      lbl <- paste0(val, star)
      col_txt <- text_color_fun(fill)
      grid.text(lbl, x, y, gp = gpar(fontsize = 10, col = col_txt))
    },
    rect_gp = gpar(col = NA),
    border = FALSE,
    row_names_side = "left",
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter",
      title_gp = gpar(fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    show_heatmap_legend = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_max_height = unit(10, "cm"),
    row_names_max_width = unit(16, "cm")

  )

  grid.newpage()
  draw(ht, heatmap_legend_side = "bottom", newpage = FALSE)
}

# Capture heatmap as graphical object
gp_corr_clustered_variables <- grid.grabExpr(create_heatmap_corr_trig_clustered_variables())

# Combine title and heatmap
corr_plot_clustered_variables <- plot_grid(
  p_title,
  gp_corr_clustered_variables,
  ncol = 1,
  rel_heights = c(0.05, 1),
  align = "v"
) + coord_fixed(ratio = 1)

print(corr_plot_clustered_variables)
ggsave("trigeminal_correlation_heatmap_clustered_variables.svg", corr_plot_clustered_variables, width = 16, height = 20)


# ========================================================================== #
# 10. COMBINE PCA AND CLUSTER PUBLICATION PLOT
# ========================================================================== #


# Combine all main PCA plots into publication-ready figure
combined_TriFunQ_PCA_plot <- cowplot::plot_grid(
# Upper row: Biplot (A)
  cowplot::plot_grid(
    pBiplot,
    pScree,
    labels = LETTERS[1:2],
    ncol = 1, rel_heights = c(2,1)
  ),
# Lower row: PCA plots (B, C and D)
    cowplot::plot_grid(
      pContrib_dim1,
      pContrib_dim2,
      labels = LETTERS[3:4],
      align = "h", axis = "t",
      nrow = 1
    ),
  ncol = 2,
  rel_widths = c(1, 1)
)

# Display combined figure
print(combined_TriFunQ_PCA_plot)

# Save combined figure
ggsave(
  "combined_TriFunQ_PCA_plot.svg",
  combined_TriFunQ_PCA_plot,
  width = 20,
  height = 16
)



# Combine all main plots into publication-ready figure
combined_TriFunQ_clustering_plot <- cowplot::plot_grid(
# Left column: Cluster analysis plot (A)
  cowplot::plot_grid(
    clustering_heatmap,
    labels = LETTERS[1],
    label_y = 0.98
  ),
# Right column, upper row: Cluster diagnostics (B, C and D)
  cowplot::plot_grid(
    cowplot::plot_grid(
    p_TriFunQ_clusters_dend,
    p_TriFunQ_clusters_silhouette,
    p_TriFunQ_clusters_cluster,
      labels = LETTERS[2:4],
      nrow = 1,
      label_y = 0.98
    ),
# Right column, lower row: Cluster variables statistics (E)
  cowplot::plot_grid(
    p_trigeminal_clustered_data,
    labels = LETTERS[5],
    label_y = 0.98
  ),
  ncol = 1,
  align = "v",
  axis = "lr"
),
cowplot::plot_grid(
  p_cohens_d,
  labels = LETTERS[6],
  label_y = 0.98
),
ncol = 3,
rel_widths = c(1, 5, 2),
align = "h",
axis = "tb"
)

# Display combined figure
print(combined_TriFunQ_clustering_plot)

# Save combined figure
ggsave(
  "combined_TriFunQ_clustering_plot.svg",
  combined_TriFunQ_clustering_plot,
  width = 40,
  height = 12
)


# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #

cat("\nAnalysis complete!\n")
