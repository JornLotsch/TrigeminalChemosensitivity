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
# LOAD REQUIRED LIBRARIES
# ========================================================================== #

library(boot) # Bootstrap methods
library(circlize) # Color mapping for heatmaps
library(cluster) # Clustering algorithms
library(ComplexHeatmap) # Advanced heatmap visualization
library(cowplot) # Plot composition and alignment
library(patchwork) # Plot composition and alignment
library(cvms) # Confusion matrix visualization
library(DataVisualizations) # Pareto Density Estimation
library(dplyr) # Data manipulation
library(effsize) # Effect size calculations
library(factoextra) # Factor visualization
library(FactoMineR) # Factor analysis
library(forcats) # Factor handling
library(ggplot2) # Data visualization
library(ggplotify) # Convert base plots to ggplot
library(ggpmisc) # Additional ggplot2 functionality
library(ggpubr) # Publication-ready plots
library(ggthemes) # Extra themes for ggplot2
library(ggrepel) # For text arrangement in plots
library(grid) # Grid graphics
library(gridExtra) # Arrange multiple grid-based plots
library(Hmisc) # Statistical tools
library(lubridate) # Date/time handling
library(MASS) # Robust statistical methods
library(missForest) # Imputation methods
library(NbClust) # Cluster number detection
library(opGMMassessment) # Gaussian mixture analysis
library(pbmcapply) # Parallel apply functions
library(psych) # Correlation analysis
library(purrr) # Functional programming tools
library(reshape2) # Data reshaping
library(scales) # Scale functions for visualization
library(stringr) # String manipulation
library(tidyr) # Data tidying
library(vcd) # Categorical data visualization
library(viridis) # Color scales for plots
library(stringr) # For strin manipulation

source("globals.R")

# ========================================================================== #
# GLOBAL OPTIONS AND ANALYSIS SWITCHES
# ========================================================================== #

remove_censored <- FALSE
scale_0_100 <- FALSE
analyze_only_untransformed <- FALSE
plot_only_untransformed <- TRUE

# ========================================================================== #
# HELPER FUNCTIONS
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
#'
#' @param y Transformed values
#' @param base Logarithm base used in forward transformation
#' @return Original scale values
inv_slog <- function(y, base = 10) {
  s <- sign(y)
  absY <- abs(y)

  if (base == 0) {
    val <- expm1(absY)
  } else if (base == 2) {
    val <- 2 ^ absY - 1
  } else if (base == 10) {
    val <- 10 ^ absY - 1
  } else {
    val <- expm1(absY / log(base))
  }
  return(s * val)
}

#' Reflected logarithmic transformation
#'
#' @param x Numeric vector
#' @return Reflected and log-transformed values
reflect_slog <- function(x) {
  slog(max(x, na.rm = TRUE) + 1 - x)
}

#' Reflected logarithmic transformation (unflipped)
#'
#' @param x Numeric vector
#' @return Reflected and log-transformed values (negative)
reflect_slog_unflipped <- function(x) {
  -slog(max(x, na.rm = TRUE) + 1 - x)
}

#' Inverse of reflect_slog_unflipped transformation
#'
#' @param y Transformed values
#' @param original_max Maximum value from original data
#' @param base Logarithm base
#' @return Original scale values
inv_reflect_slog_unflipped <- function(y, original_max, base = 10) {
  M <- original_max
  x_original <- M + 1 - inv_slog(-y, base)
  return(x_original)
}

#' Create Pareto Density Estimation plot with ggplot2
#'
#' @param Data Numeric matrix or data frame
#' @return ggplot object with PDE curves
PDEplotGG <- function(Data) {
  Data <- as.matrix(Data)
  m <- matrix(NA, nrow = 0, ncol = 3)

  for (i in seq_len(ncol(Data))) {
    PDE <- ParetoDensityEstimation(as.vector(na.omit(Data[, i])))
    m2 <- PDE$kernels
    m3 <- PDE$paretoDensity
    m1 <- rep(i, length(m2))
    m <- rbind(m, cbind(m1, m2, m3))
  }

  mdf <- data.frame(m)

  p <- ggplot(data = mdf, aes(x = m2, y = m3, colour = factor(m1))) +
    geom_line(aes(linewidth = 1)) +
    guides(linewidth = FALSE)

  return(p)
}

#' Convert p-values to significance stars
#'
#' @param p P-value
#' @return Character string with significance stars
signif_code <- function(p) {
  if (is.na(p)) ""
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

#' Determine text color based on background brightness
#'
#' @param fill_color Background color
#' @return Text color (black or white) for optimal contrast
text_color_fun <- function(fill_color) {
  rgb_val <- col2rgb(fill_color) / 255
  brightness <- 0.299 * rgb_val[1,] + 0.587 * rgb_val[2,] + 0.114 * rgb_val[3,]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}

#' Compute Cohen's d effect size
#'
#' @param x First group values
#' @param y Second group values
#' @return Cohen's d estimate with Hedges correction
cohen_d_val <- function(x, y) {
  cohen.d(x, y, hedges.correction = TRUE)$estimate
}

#' Bootstrap confidence intervals for Cohen's d
#'
#' @param x First group values
#' @param y Second group values
#' @param R Number of bootstrap replicates
#' @return Vector with lower and upper CI bounds
boot_cohen_d <- function(x, y, R = 1000) {
  data <- data.frame(
    group = rep(c("x", "y"), times = c(length(x), length(y))),
    value = c(x, y)
  )

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

#' Bootstrap confidence intervals for correlation
#'
#' @param x First variable
#' @param y Second variable
#' @param R Number of bootstrap replicates
#' @return Vector with median correlation and CI bounds
bootstrap_cor <- function(x, y, R = 1000) {
  data <- data.frame(x = x, y = y)

  boot_cor <- function(data, indices) {
    d <- data[indices,]
    if (sum(complete.cases(d)) < 4) return(NA)
    cor(d$x, d$y, method = "spearman", use = "complete.obs")
  }

  boot_obj <- boot(data, boot_cor, R = R)
  ci <- boot.ci(boot_obj, type = "perc")$percent[4:5]
  median_cor <- median(boot_obj$t[!is.na(boot_obj$t)], na.rm = TRUE)
  c(median_cor, ci[1], ci[2])
}

#' Compute rank-biserial correlation (companion to Wilcoxon test)
#'
#' @param x First group values
#' @param y Second group values
#' @return Rank-biserial correlation coefficient (r)
#' @details
#' Rank-biserial correlation is a non-parametric effect size measure
#' for Mann-Whitney U / Wilcoxon rank-sum test.
#' Formula: r = 1 - (2U)/(n1*n2)
#' Interpretation: |r| = 0.1 (small), 0.3 (medium), 0.5 (large)
rank_biserial_val <- function(x, y) {
  wtest <- wilcox.test(x, y, exact = FALSE)
  U <- as.numeric(wtest$statistic)
  n1 <- length(x)
  n2 <- length(y)
  r <- 1 - (2 * U) / (n1 * n2)
  return(r)
}

#' Bootstrap confidence intervals for rank-biserial correlation
#'
#' @param x First group values
#' @param y Second group values
#' @param R Number of bootstrap replicates (default: 1000)
#' @return Vector with lower and upper CI bounds
boot_rank_biserial <- function(x, y, R = 1000) {
  data <- data.frame(
    group = rep(c("x", "y"), times = c(length(x), length(y))),
    value = c(x, y)
  )

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

# ========================================================================== #
# READ DATA
# ========================================================================== #

trigeminale_data_raw <- read.csv("trigeminale_daten_corrected_translated.csv", check.names = FALSE)
trigeminale_data <- read.csv("analysis_dataset_imputed.csv", check.names = FALSE)

# Check if raw and imputed data are in the same order of cases
plot(as.numeric(trigeminale_data_raw$ID)~ as.numeric(trigeminale_data$ID))


# ========================================================================== #
# INITIAL DATA VISUALIZATION - HEATMAP OVERVIEW
# ========================================================================== #

# Extract and rename relevant trigeminal measures
trigeminal_measures <- c("AmmoLa intensity", "Lateralization (x/20)", "CO2 threshold")

trigeminal_measures_data <- trigeminale_data_raw[, trigeminal_measures]
names(trigeminal_measures_data) <- c("AmmoLa_intensity", "Lateralization", "CO2_threshold")

# Report sample sizes for each measure
cat("\nSample sizes per measure:\n")
print(apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x))))


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

###########################################################################################################################
#####  Analysis on original non-imputed data ##############################################################################
###########################################################################################################################

# ========================================================================== #
# DATA QUALITY ASSESSMENT - CENSORING ANALYSIS
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
# GROUP AGREEMENT ANALYSIS - TOP 10% SENSITIVITY
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
# CONFUSION MATRIX ANALYSIS
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
# VARIABLE TRANSFORMATION AND ALTERNATIVE GMM BASED GROUPING ANALYSIS
# ========================================================================== #

print("Minimum value in AmmoLa_intensity")
print(min(trigeminal_measures_data$AmmoLa_intensity))

if (!analyze_only_untransformed) {
  cat("\nApplying variable transformations...\n")

  # Apply transformations to address skewness
  trigeminal_measures_data <- trigeminal_measures_data %>%
    mutate(
      AmmoLa_intensity_reflect_slog = reflect_slog_unflipped(AmmoLa_intensity),
      Lateralization = Lateralization,
      CO2_threshold_log = -slog(CO2_threshold) # replace regular log by slog for consistency
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
  trigeminal_measures_data$Lateralization,
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

trigeminal_measures_data_grouped_GMM <- trigeminal_measures_data[, c(2,4,5)]
trigeminal_measures_data_grouped_GMM$CO2_threshold_log[1:549] <- NA

# Create binary grouping: Upper mode = most sensitive vs others
trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog <- ifelse(
  trigeminal_measures_data_grouped_GMM$AmmoLa_intensity_reflect_slog > tail(AmmoLa_GMM$Boundaries), 1, 0
)
trigeminal_measures_data_grouped_GMM$Lateralization <- ifelse(
  trigeminal_measures_data_grouped_GMM$Lateralization > tail(Lateralization_GMM$Boundaries), 1, 0
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
  trigeminal_measures_data_grouped_GMM$Lateralization
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
# AGE AND SEX EFFECTS ANALYSIS
# ========================================================================== #

cat("\nAnalyzing age correlations...\n")

# Add age to dataset and calculate correlations
trigeminal_measures_data_age <- trigeminal_measures_data
trigeminal_measures_data_age$age <- as.numeric(trigeminale_data$Age)
corr_mat_age <- Hmisc::rcorr(
  as.matrix(trigeminal_measures_data_age),
  type = "spearman"
)

cat("\nAnalyzing sex differences...\n")

# Add sex to dataset
trigeminal_measures_data_sex <- trigeminal_measures_data
trigeminal_measures_data_sex$sex <- as.factor(trigeminale_data_raw$Gender)

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
# DISTRIBUTION PLOTS - HISTOGRAMS AND DENSITY ESTIMATES
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
# DENSITY PLOTS
# ========================================================================== #

cat("\nGenerating Pareto Density Estimation plots...\n")

# Density plot for AmmoLa intensity

# Overlay density plot for AmmoLa intensity
p_density_AmmoLA <- ggplot(trigeminal_measures_data, aes(x = AmmoLa_intensity)) +
  geom_density(alpha = 0.2, size = 1, fill = "cornsilk4", color = "cornsilk4") +
  labs(title = "Distribution of AmmoLa intensity estimates",
       x = "Rating [%]", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_density_AmmoLA)

# PDE plot for AmmoLa intensity
pPDE_AmmoLA <- PDEplotGG(trigeminal_measures_data$AmmoLa_intensity) +
  theme_plot() +
  labs(title = "Distribution of AmmoLa intensity estimates", x = "Rating [%]", y = "PDE") +
  guides(color = "none") +
  scale_color_manual(values = "cornsilk4")

print(pPDE_AmmoLA)

# Overlay density plot for AmmoLa intensity transformed
p_density_AmmoLA_transformed <- ggplot(trigeminal_measures_data, aes(x = AmmoLa_intensity_reflect_slog)) +
  geom_density(alpha = 0.2, size = 1, fill = "cornsilk4", color = "cornsilk4") +
  labs(title = "Distribution of AmmoLa (transformed)", x = "Reflected log of rating [units]", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_density_AmmoLA_transformed)


# Overlay density plot for Lateralization
p_density_Lateralization <- ggplot(trigeminal_measures_data, aes(x = Lateralization)) +
  geom_density(alpha = 0.2, size = 1, fill = "cornsilk4", color = "cornsilk4") +
  labs(title = "Distribution of lateralization successes",
       x = "Correct [count]", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_density_Lateralization)

# PDE plot for Lateralization
pPDE_Lateralization <- PDEplotGG(trigeminal_measures_data$Lateralization) +
  theme_plot() +
  labs(title = "Distribution of correct lateralizations", x = "Lateralization [n correct]", y = "PDE") +
  guides(color = "none") +
  scale_color_manual(values = "cornsilk4")

print(pPDE_Lateralization)


# Overlay density plot for CO2 thresholds (split by breath protocol)

CO2_threshold_breathing <- trigeminal_measures_data$CO2_threshold[1:549]
CO2_threshold_breath_hold <- trigeminal_measures_data$CO2_threshold[550:length(trigeminal_measures_data$CO2_threshold)]

df_plot_CO2_threshold_breathing <- data.frame(
  value = c(CO2_threshold_breathing, CO2_threshold_breath_hold),
  Breath = rep(c("Uncontrolled", "Hold"), c(length(CO2_threshold_breathing), length(CO2_threshold_breath_hold)))
)

p_breathing_CO2_thresholds <- ggplot(df_plot_CO2_threshold_breathing, aes(x = value, fill = Breath, color = Breath)) +
  geom_density(alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("cornsilk2", "cornsilk4")) +
  scale_color_manual(values = c("cornsilk2", "cornsilk4")) +
  labs(title = "CO2 Threshold: Breath uncontrolled vs hold distr.",
       x = "CO2 threshold (ms)", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_breathing_CO2_thresholds)

# PDE plot for CO2 thresholds (split by breath protocol)

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
  labs(title = "Distribution of CO2 thresholds", color = "Breath", x = "CO2 threshold (ms)", y = "PDE") +
  scale_color_manual(
    values = c("grey83", "cornsilk4"),
    labels = c("Uncontrolled", "Hold")
  ) +
  theme(legend.position.inside = TRUE, legend.position = c(.2, .85))

print(pPDE_CO2)


# Overlay density plot for CO2 thresholds (only breath hold)
p_density_CO2_thresholds_transformed <- ggplot(trigeminal_measures_data[550:length(trigeminal_measures_data$CO2_threshold),], aes(x = CO2_threshold_log)) +
  geom_density(alpha = 0.2, size = 1, fill = "cornsilk4", color = "cornsilk4") +
  labs(title = "Distribution of CO2 thresholds (transformed)", x = "log CO2 threshold (log ms)", y = "Density") +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.7, .8), legend.direction = "horizontal")

print(p_density_CO2_thresholds_transformed)


## Combine distribution plots

library(patchwork)

# 1. Extract all plots with uniform margins
p1 <- p_density_AmmoLA + theme(plot.margin = margin(5, 5, 5, 5))
p2 <- p_density_Lateralization + theme(plot.margin = margin(5, 5, 5, 5))
p3 <- p_breathing_CO2_thresholds + theme(plot.margin = margin(5, 5, 5, 5))

p4 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "No imputation\nneeded\n(complete data)",
           size = 4, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

p5 <- p_observed_imputed_Lateralization + theme(plot.margin = margin(5, 5, 5, 5))
p6 <- p_observed_imputed_CO2_thresholds + theme(plot.margin = margin(5, 5, 5, 5))

p7 <- p_density_AmmoLA_transformed + theme(plot.margin = margin(5, 5, 5, 5))

p8 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "No transformation\n(count data)",
           size = 4, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

p9 <- p_density_CO2_thresholds_transformed + theme(plot.margin = margin(5, 5, 5, 5))

# 2. Build layout with panel labels (skip placeholders p4 and p8)
p_combined_psychophysical_trig_measures_obs <-
  wrap_plots(
    p1, p2, p3,
    p4, p5, p6,
    p7, p8, p9,
    ncol = 3, nrow = 3,
    heights = c(1, 1, 1),
    widths = c(1, 1, 1),
    tag_level = 'new'
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

# Add tags manually to non-placeholder plots
p_combined_psychophysical_trig_measures_obs[[1]] <- p_combined_psychophysical_trig_measures_obs[[1]] + ggtitle("A") + theme(plot.title = element_text(face = "bold", hjust = 0))
p_combined_psychophysical_trig_measures_obs[[2]] <- p_combined_psychophysical_trig_measures_obs[[2]] + ggtitle("B") + theme(plot.title = element_text(face = "bold", hjust = 0))
p_combined_psychophysical_trig_measures_obs[[3]] <- p_combined_psychophysical_trig_measures_obs[[3]] + ggtitle("C") + theme(plot.title = element_text(face = "bold", hjust = 0))
# Skip p4 (placeholder)
p_combined_psychophysical_trig_measures_obs[[5]] <- p_combined_psychophysical_trig_measures_obs[[5]] + ggtitle("D") + theme(plot.title = element_text(face = "bold", hjust = 0))
p_combined_psychophysical_trig_measures_obs[[6]] <- p_combined_psychophysical_trig_measures_obs[[6]] + ggtitle("E") + theme(plot.title = element_text(face = "bold", hjust = 0))
p_combined_psychophysical_trig_measures_obs[[7]] <- p_combined_psychophysical_trig_measures_obs[[7]] + ggtitle("F") + theme(plot.title = element_text(face = "bold", hjust = 0))
# Skip p8 (placeholder)
p_combined_psychophysical_trig_measures_obs[[9]] <- p_combined_psychophysical_trig_measures_obs[[9]] + ggtitle("G") + theme(plot.title = element_text(face = "bold", hjust = 0))

p_combined_psychophysical_trig_measures_obs <- p_combined_psychophysical_trig_measures_obs +
  plot_annotation(
    title = "Trigeminal psychophysics workflow",
    subtitle = "Row 1: Raw | Row 2: Imputation check | Row 3: Transformed",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0)
    )
  )

print(p_combined_psychophysical_trig_measures_obs)

ggsave("p_combined_psychophysical_trig_measures_obs.svg", p_combined_psychophysical_trig_measures_obs, width = 15, height = 12, dpi = 300)

# ========================================================================== #
# CORRELATION MATRIX PREPARATION
# ========================================================================== #

cat("\nPreparing correlation matrix...\n")

# Handle transformed variables if enabled
if (!analyze_only_untransformed) {
  trigeminal_measures_data_transformed <- trigeminal_measures_data[, c("AmmoLa_intensity_reflect_slog", "Lateralization", "CO2_threshold_log")]

  # Invert CO2 so higher = more sensitive
  trigeminal_measures_data_transformed$CO2_threshold_log <- -trigeminal_measures_data_transformed$CO2_threshold_log

  # Remove breath-not-controlled CO2 data
  trigeminal_measures_data_transformed <- trigeminal_measures_data_transformed %>%
    mutate(CO2_threshold_log = ifelse(row_number() <= 549, NA, CO2_threshold_log))

  # Define variable pairs for correlation analysis
  var_pairs <- list(
    c("AmmoLa_intensity_reflect_slog", "Lateralization"),
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
# CORRELATION MATRIX CALCULATION AND HEATMAP
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

# Classes
trigeminal_ammola_high_low_data <- cbind.data.frame(Class = as.factor(trigeminal_measures_data_grouped$AmmoLa_intensity),
                                                    trigeminale_data[, names(trigeminale_data) %in%
                                                                       c(variables_by_categories$Nasal_chemosensory_perception,
                                                                         variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])])

head(trigeminal_ammola_high_low_data)
trigeminal_ammola_high_low_data$`CO2 threshold`[1:549] <- NA
trigeminal_ammola_high_low_data$`CO2 threshold` <- -trigeminal_ammola_high_low_data$`CO2 threshold`

trigeminal_ammola_high_low_data_long <- reshape2::melt(trigeminal_ammola_high_low_data, id.vars = "Class")
head(trigeminal_ammola_high_low_data_long)

round(apply(trigeminal_ammola_high_low_data[, -1], 2, function(x) wilcox.test(x ~ trigeminal_ammola_high_low_data$Class, na.rm = TRUE)$p.value), 2)


# Dodge width for positioning (same as in jitterdodge and boxplot)
dodge_width <- 0.8

# Calculate the medians for each (Class, variable) combination
median_data <- trigeminal_ammola_high_low_data_long %>%
  group_by(variable, Class) %>%
  dplyr::summarise(median_value = mean(value, na.rm = T), .groups = "drop")
median_data$x_position <- rep(c(0.8, 1.2), nrow(median_data) / 2)

# Plot data and add a line for medians
set.seed(42)
p_trigeminal_ammola_high_low_data <- ggplot(trigeminal_ammola_high_low_data_long, aes(x = variable, y = value, color = Class, fill = Class)) +
  geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = dodge_width)) +
  geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = dodge_width), outlier.shape = NA) +
  geom_jitter(alpha = 1, size = 0.2,
              position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = dodge_width)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                             label.y.npc = "top",
                             vjust = -0.2, size = 3) +
# Add the median lines within dodge positions
geom_line(data = median_data, aes(x = x_position, y = median_value, group = variable),
            color = "black", size = 0.8, linetype = "dashed", inherit.aes = FALSE) +
  facet_wrap(variable ~ ., nrow = 3, scales = "free", labeller = label_wrap_gen(width = 20)) +
  guides(color = "none") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk4", "cornsilk4")) +
  labs(title = "Raw trigeminal data per variable and AmmoLa sensitvity class", fill = "Class", x = NULL) +
  theme_plot() +
  theme(
    legend.direction = "horizontal",
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(size = 6, lineheight = 0.9, face = "plain"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5) ,
    axis.text.x = element_blank(), # Remove x axis tick labels
    axis.ticks.x = element_blank() # Remove x axis ticks (optional)
  ) +
  guides(color = "none")

print(p_trigeminal_ammola_high_low_data)



# Define helper function to compute Cohen's d for two numeric vectors
cohen_d_val <- function(x, y) {
  effsize::cohen.d(x, y, hedges.correction = TRUE)$estimate
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

  x <- data_var$value[data_var$Class == levels(data_var$Class)[1]]
  y <- data_var$value[data_var$Class == levels(data_var$Class)[2]]

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
  geom_hline(yintercept = -0.8, linetype = "dashed", color = "salmon") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35))


p_cohens_d_AmmoLa_gigh_low

# Calculate rank-biserial correlation and CIs per variable using map
# Note: rank_biserial_val() and boot_rank_biserial() are defined in the helper functions section
results_rank_biserial_AmmoLa <- map_df(variables, function(var) {
  data_var <- trigeminal_ammola_high_low_data_long %>%
    filter(variable == var) %>%
    filter(!is.na(value))

  x <- data_var$value[data_var$Class == levels(data_var$Class)[1]]
  y <- data_var$value[data_var$Class == levels(data_var$Class)[2]]

  # Skip variables with insufficient data
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(variable = var, r = NA, ci_lower = NA, ci_upper = NA))
  }

  r <- rank_biserial_val(x, y)
  ci <- boot_rank_biserial(x, y)

  tibble(variable = var, r = r, ci_lower = ci[1], ci_upper = ci[2])
})

# Plot rank-biserial correlation with CI error bars
p_rank_biserial_AmmoLa_high_low <- ggplot(results_rank_biserial_AmmoLa, aes(x = reorder(variable, r), y = r)) +
  geom_col(fill = "cornsilk3") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  coord_flip() +
  labs(title = "Rank-biserial correlation (95% CI)",
       y = "Rank-biserial r",
       x = "Variables") +
  theme_plot() +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "salmon") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "salmon") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35))

p_rank_biserial_AmmoLa_high_low


# ========================================================================== #
# COMBINED PUBLICATION FIGURE
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
rel_widths = c(5, 2),
align = "h",
axis = "tb"
)

# Display combined figure
print(combined_trigeminal_analysis_plot)

# Save combined figure
ggsave(
  "combined_trigeminal_analysis_plot.svg",
  combined_trigeminal_analysis_plot,
  width = 30,
  height = 18
)


library(patchwork)
library(dplyr)

# Reorder violin plot to match Cohen's d (largest |d| first)
p_cohens_d_ordered <- p_cohens_d_AmmoLa_gigh_low  # already ordered by effect size

# Get the variable order from Cohen's d plot
variable_order <- results_cohen_AmmoLa %>%
  arrange(desc(abs(d))) %>%
  pull(variable) %>%
  as.character()

# Reorder violin plot factors to match Cohen's d order
trigeminal_ammola_high_low_data_long$variable <- factor(trigeminal_ammola_high_low_data_long$variable,
                                                        levels = variable_order)

# Recreate violin plot with new order
p_trigeminal_ordered <- ggplot(trigeminal_ammola_high_low_data_long, aes(x = variable, y = value, color = Class, fill = Class)) +
  geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(alpha = 1, size = 0.2,
              position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.8)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                             label.y.npc = "top",
                             vjust = -0.2, size = 3) +
  facet_wrap(variable ~ ., nrow = 3, scales = "free", labeller = label_wrap_gen(width = 20)) +
  guides(color = "none") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk4", "cornsilk4")) +
  labs(title = "Distributions by AmmoLa sensitivity class", fill = "Class", x = NULL) +
  theme_plot() +
  theme(
    legend.direction = "horizontal",
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(size = 6, lineheight = 0.9, face = "plain"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5) ,
    axis.text.x = element_blank(), # Remove x axis tick labels
    axis.ticks.x = element_blank() # Remove x axis ticks (optional)
  ) +
  guides(color = "none")

# Combine: Cohen's d (left, flipped) + violin plots (right)
p_combined_effect_sizes <- cowplot::plot_grid(
  p_cohens_d_ordered,
  p_trigeminal_ordered ,
  labels = "AUTO",
  align = "h", axis = "tb",
  rel_widths = c(1,2)
)
print(p_combined_effect_sizes)
ggsave("p_trigeminal_ammola_effect_sizes.svg", p_combined_effect_sizes, width = 17, height = 15, dpi = 300,limitsize = FALSE)
ggsave("p_trigeminal_ammola_effect_sizes.png", p_combined_effect_sizes, width = 17, height = 15, dpi = 300,limitsize = FALSE)

# ========================================================================== #
# RANK-BISERIAL VERSION: Alternative combined figure with rank-biserial
# ========================================================================== #

# Create combined figure with rank-biserial instead of Cohen's d
combined_trigeminal_analysis_plot_rb <- cowplot::plot_grid(
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
  p_rank_biserial_AmmoLa_high_low,
  labels = LETTERS[8]),
ncol = 2,
rel_widths = c(5, 2),
align = "h",
axis = "tb"
)

# Display combined figure with rank-biserial
print(combined_trigeminal_analysis_plot_rb)

# Save combined figure with rank-biserial
ggsave(
  "combined_trigeminal_analysis_plot_rank_biserial.svg",
  combined_trigeminal_analysis_plot_rb,
  width = 30,
  height = 18
)

# Reorder violin plot to match rank-biserial (largest |r| first)
p_rank_biserial_ordered <- p_rank_biserial_AmmoLa_high_low  # already ordered by effect size

# Get the variable order from rank-biserial plot
variable_order_rb <- results_rank_biserial_AmmoLa %>%
  arrange(desc(abs(r))) %>%
  pull(variable) %>%
  as.character()

# Reorder violin plot factors to match rank-biserial order
trigeminal_ammola_high_low_data_long_rb <- trigeminal_ammola_high_low_data_long
trigeminal_ammola_high_low_data_long_rb$variable <- factor(trigeminal_ammola_high_low_data_long_rb$variable,
                                                           levels = variable_order_rb)

# Recreate violin plot with new order (rank-biserial)
p_trigeminal_ordered_rb <- ggplot(trigeminal_ammola_high_low_data_long_rb, aes(x = variable, y = value, color = Class, fill = Class)) +
  geom_violin(alpha = 0.05, width = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.2, width = 0.3, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(alpha = 1, size = 0.2,
              position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.8)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                             label.y.npc = "top",
                             vjust = -0.2, size = 3) +
  facet_wrap(variable ~ ., nrow = 3, scales = "free", labeller = label_wrap_gen(width = 20)) +
  guides(color = "none") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk4", "cornsilk4")) +
  labs(title = "Distributions by AmmoLa sensitivity class", fill = "Class", x = NULL) +
  theme_plot() +
  theme(
    legend.direction = "horizontal",
    legend.position.inside = TRUE, legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.6), color = NA),
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(size = 6, lineheight = 0.9, face = "plain"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5) ,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(color = "none")

# Combine: rank-biserial (left, flipped) + violin plots (right)


# Uniform margins only
p_left <- p_rank_biserial_ordered + theme(plot.margin = margin(5, 5, 5, 5))
p_right <- p_trigeminal_ordered_rb + theme(plot.margin = margin(5, 5, 5, 5))

# patchwork horizontal + BIG A B LABELS + original titles preserved
p_combined_effect_sizes_rb <- (p_left | p_right) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Group differences across AmmoLa sensitivity classes",
    tag_levels = list(c("A", "B")),
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

print(p_combined_effect_sizes_rb)
ggsave("p_trigeminal_ammola_effect_sizes_rank_biserial.svg", p_combined_effect_sizes_rb, width = 17, height = 15, dpi = 300, limitsize = FALSE)
ggsave("p_trigeminal_ammola_effect_sizes_rank_biserial.png", p_combined_effect_sizes_rb, width = 17, height = 15, dpi = 300, limitsize = FALSE)


###########################################################################################################################
#####  Analysis on  imputed complete data ##########################################################################################
###########################################################################################################################

# ========================================================================== #
# PCA
# ========================================================================== #

# Create data matrix
variables_for_clustering_imputed_all <- trigeminale_data[, intersect(names(trigeminale_data), variables_by_categories$Nasal_chemosensory_perception)]

heat_matrix_all <- cbind.data.frame(variables_for_clustering_imputed_all,
                                    AmmoLa_transformed =  scaleRange_01(reflect_slog_unflipped(trigeminale_data$`AmmoLa intensity`)) * 3,
                                    Lateralization =  scaleRange_01(trigeminale_data$`Lateralization (x/20)`) * 3,
                                    CO2_threshold_transformed =  scaleRange_01(-slog(trigeminale_data$`CO2 threshold`)) * 3
                                )


# Perform PCA on the imputed clustering data
res.pca_trigeminal_clustered_data <- FactoMineR::PCA(heat_matrix_all, scale.unit = TRUE, graph = TRUE, ncp = Inf)
res.pca_trigeminal_clustered_data_biplot <- stats::princomp(scale(heat_matrix_all, center = TRUE, scale = TRUE)) #, scale.unit = TRUE, graph = TRUE, ncp = 8)

# Extract variable coordinates (variables x components)
var_coords <- res.pca_trigeminal_clustered_data$var$coord

# Create a summary table with eigenvalues, standard deviations, variance explained, cumulative variance
eig <- res.pca_trigeminal_clustered_data$eig
std_dev <- sqrt(eig[, 1])
prop_var <- eig[, 2] / 100
cum_var <- eig[, 3] / 100
n_comp <- length(which(eig[,1] >= 1))

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
write.csv(combined_table[, 1:n_comp], "res.pca_trigeminal_clustered_data_loadings.csv")

# Display eigenvalues and variable contributions
cat("PCA Eigenvalues:\n")
print(res.pca_trigeminal_clustered_data$eig)
cat("PCA Variable Contributions:\n")
print(res.pca_trigeminal_clustered_data$var)

df_pca_voronoi <- cbind.data.frame(Class = trigeminal_ammola_high_low_data$Class, res.pca_trigeminal_clustered_data$ind$coord[,1:2])

p_pca_voronoi <- VoronoiBiomedPlot::create_voronoi_plot(data = df_pca_voronoi, class_column = "Class") +
  scale_color_manual(values = c("cornsilk3", "cornsilk4")) +
  scale_fill_manual(values = c("cornsilk", "cornsilk4"))

print(p_pca_voronoi)

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

# 1. Extract loadings and compute contribution (squared loadings on Dim1 + Dim2)
loadings <- res.pca_trigeminal_clustered_data_biplot$loadings
var.coord <- as.data.frame(loadings[, 1:2])
names(var.coord) <- c("Dim1", "Dim2")

# contribution per variable (same as fviz_pca_var uses internally)
var.contrib <- rowSums(loadings[, 1:2]^2)
var.coord$contrib <- var.contrib

# 2. Wrap long labels
long_names <- rownames(loadings)
var.coord$label <- str_wrap(long_names, width = 22)

# 3. Build the biplot from scratch (arrows + labels, same gradient)
pBiplot <-
  ggplot(var.coord, aes(Dim1, Dim2, color = contrib)) +
  # Arrows
  geom_segment(
    aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
    arrow = arrow(length = unit(0.2, "cm")),
    alpha = 0.7
  ) +
  # Labels colored by the same contribution
  geom_text_repel(
    aes(label = label),
    segment.alpha = 0.5,
    size = 3,
    direction = "both",
    max.iter = 5000,
    nudge_y = 0.05
  ) +
  # Use exactly your gradient
  scale_color_gradientn(
    colours = c("grey20", "sienna3", "red"),
    name = "Contribution"
  ) +
  theme_plot() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position.inside = TRUE,
    legend.position = c(.2, .2)
  ) +
  labs(title = "PCA biplot")

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
cor_results_clustered_variables <- corr.test(heat_matrix_all, method = cor_method)
print(cor_results_clustered_variables)

corr_mat_clustered_variables <- cor_results_clustered_variables$r # Correlation coefficients
p_mat_clustered_variables <- cor_results_clustered_variables$p # P-values


# Create correlation heatmap with ComplexHeatmap
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


p_title_corr <- ggplot() +
  labs(title = paste0("Correlation matrix (", cor_method, ")")) +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "plain", size = 12, color = "#222222",
      hjust = 0, margin = margin(b = 5, t = 5)
    ),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 20, 0, 10, unit = "pt") # Much smaller margins
  )

# Combine title and heatmap
corr_plot_clustered_variables <- plot_grid(
  p_title_corr,
  gp_corr_clustered_variables,
  ncol = 1,
  rel_heights = c(0.03, 1),
  align = "v"
) #+ coord_fixed(ratio = 1)

print(corr_plot_clustered_variables)
# ggsave("trigeminal_correlation_heatmap_clustered_variables.svg", corr_plot_clustered_variables, width = 16, height = 20)

# Bootstrap function (unchanged)
bootstrap_cor <- function(x, y, R = 1000) {
  data <- data.frame(x = x, y = y)
  boot_cor <- function(data, indices) {
    d <- data[indices,]
    if (sum(complete.cases(d)) < 4) return(NA)
    cor(d$x, d$y, method = "spearman", use = "complete.obs")
  }
  boot_obj <- boot(data, boot_cor, R = R)
  ci <- boot.ci(boot_obj, type = "perc")$percent[4:5]
  median_cor <- median(boot_obj$t[!is.na(boot_obj$t)], na.rm = TRUE)
  c(median_cor, ci[1], ci[2])
}

# Fix column indexing + PARALLEL VERSION
variables <- colnames(heat_matrix_all)
ci_results <- data.frame(
  Var1 = character(), Var2 = character(), Median = numeric(),
  CI_lower = numeric(), CI_upper = numeric(), n_complete = integer(),
  stringsAsFactors = FALSE
)

# Parallel bootstrap (Linux-friendly)
cor_pairs <- combn(seq_along(variables), 2, simplify = FALSE)

ci_results <- pbmclapply(cor_pairs, function(pair_idx) {
  i <- pair_idx[1];
  j <- pair_idx[2]
  x_vals <- heat_matrix_all[[variables[i]]] # [[ ]] fixes column selection
  y_vals <- heat_matrix_all[[variables[j]]]
  boot_result <- bootstrap_cor(x_vals, y_vals)

  data.frame(
    Var1 = variables[i], Var2 = variables[j],
    Median = boot_result[1], CI_lower = boot_result[2],
    CI_upper = boot_result[3],
    n_complete = sum(complete.cases(x_vals, y_vals)),
    stringsAsFactors = FALSE
  )
}, mc.cores = parallel::detectCores() - 1) %>% bind_rows() # Adjust mc.cores as needed

# FIXED plot_df - proper sorting that survives coord_flip()
plot_df <- ci_results %>%
  filter(!is.na(Median)) %>%
  mutate(
    Pair = paste0(Var1, " - ", Var2),
    abs_median = abs(Median)
  ) %>%
  arrange(desc(abs_median)) %>%
# CRITICAL: Set levels in REVERSE order for coord_flip()
mutate(Pair = factor(Pair, levels = rev(Pair))) %>%
  ungroup()

# Verify sorting is correct
print(head(plot_df[, c("Pair", "Median", "abs_median")], 200))

# Plot + Cohen (1988) thresholds
p_cor_ci <- ggplot(plot_df, aes(x = Pair, y = Median)) +
  geom_vline(
    xintercept = seq(1.5, length(unique(plot_df$Pair)) - 0.5, by = 1),
    color = "grey90", linewidth = 0.3, alpha = 0.5
  ) +
  # Cohen (1988) thresholds
  geom_hline(yintercept = 0.5,  color = "salmon", linetype = "dashed",
             linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0.3,  color = "orange", linetype = "dashed",
             linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0.1,  color = "gold",   linetype = "dashed",
             linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -0.1, color = "gold",   linetype = "dashed",
             linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -0.3, color = "orange", linetype = "dashed",
             linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -0.5, color = "salmon", linetype = "dashed",
             linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "solid",
             linewidth = 0.4) +
  geom_rect(
    aes(xmin = as.numeric(Pair) - 0.45, xmax = as.numeric(Pair) + 0.45,
        ymin = CI_lower, ymax = CI_upper),
    fill = "cornsilk3", color = "black", alpha = 0.7, linewidth = 0.2
  ) +
  geom_segment(
    aes(x = as.numeric(Pair) - 0.35, xend = as.numeric(Pair) + 0.35,
        y = Median, yend = Median),
    linewidth = 1, color = "black"
  ) +
  coord_flip(ylim = c(-0.8, 0.8)) +
  labs(
    x = "Variable pairs (sorted by |ρ|, highest on top)",
    y = "Spearman ρ (bootstrap median)",
    title = "Spearman correlation CIs with Cohen (1988) effect size thresholds",
    caption = paste(
      "Dashed lines: |r| < 0.1 very small, 0.1 ≤ |r| < 0.3 small,",
      "0.3 ≤ |r| < 0.5 moderate, |r| ≥ 0.5 large (Cohen, 1988)."
    )
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    plot.margin = margin(l = 80, r = 20, t = 10, b = 40, unit = "pt"),
    plot.caption = element_text(hjust = 0, size = 8)
  )

print(p_cor_ci)

# ggsave("correlation_ci_pairs.svg", p_cor_ci, width = 14, height = 12, dpi = 300)


# Uniform margins only
p_left <- corr_plot_clustered_variables + theme(plot.margin = margin(5, 5, 5, 5))
p_right <- p_cor_ci + theme(plot.margin = margin(5, 5, 5, 5))

# patchwork horizontal + BIG A B LABELS + original titles preserved
correlation_matrix_ci_pairs <- (p_left | p_right) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Correlation of trigeminal measures (all imputed data)",
    tag_levels = list(c("A", "B")),
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))


print(correlation_matrix_ci_pairs)
ggsave("correlation_matrix_ci_pairs.svg", correlation_matrix_ci_pairs,
       width = 40, height = 20, dpi = 300, bg = "white")
ggsave("correlation_matrix_ci_pairs.png", correlation_matrix_ci_pairs,
       width = 40, height = 20, dpi = 300, bg = "white")


# ========================================================================== #
# COMBINE PCA AND CLUSTER PUBLICATION PLOT
# ========================================================================== #


# Combine all main PCA plots into publication-ready figure
combined_TriFunQ_PCA_plot <- cowplot::plot_grid(
# Upper row: Biplot (A)
  cowplot::plot_grid(
    pBiplot,
    pScree,
    labels = LETTERS[1:2],
    ncol = 1, rel_heights = c(2, 1)
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


# Combine all main PCA plots into publication-ready figure

library(patchwork)

# Column 1: p1 ABOVE p2 (vertical, patchwork style)
p_column_1 <- (pBiplot + theme(plot.margin = margin(5, 5, 5, 5))) /
  (pScree + theme(plot.margin = margin(5, 5, 5, 5)))

# Column 2: p3 BESIDE p4 (horizontal, patchwork style)
p_column_2 <- (pContrib_dim1 + theme(plot.margin = margin(5, 5, 5, 5))) |
  (pContrib_dim2 + theme(plot.margin = margin(5, 5, 5, 5)))

# Combine columns side-by-side WITH cowplot-style alignment
combined_TriFunQ_PCA_plot <- cowplot::plot_grid(
  p_column_1, p_column_2,
  labels = LETTERS[1:4],        # A,B,C,D labels
  align = "h", axis = "t",     # EXACT cowplot alignment from template
  rel_widths = c(1, 1),         # Equal column widths
  ncol = 2
)

# Add title above (cowplot style)
combined_TriFunQ_PCA_plot_final <- cowplot::plot_grid(
  ggplot() +
    labs(title = "PCA results") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "plain", hjust = 0, margin = margin(b = 10))),
  combined_TriFunQ_PCA_plot,
  ncol = 1,
  rel_heights = c(0.03, 1),
  align = "v"
)

# Display combined figure
print(combined_TriFunQ_PCA_plot)

# Save combined figure
ggsave("combined_TriFunQ_PCA_plot.svg", combined_TriFunQ_PCA_plot, width = 20, height = 16)
ggsave("combined_TriFunQ_PCA_plot.png", combined_TriFunQ_PCA_plot, width = 20, height = 16, dpi = 300)


# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #

cat("\nAnalysis complete!\n")
