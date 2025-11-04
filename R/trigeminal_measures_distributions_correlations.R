################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Joern Pons-Kuhnemann
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
# 1. LOAD REQUIRED LIBRARIES
# ========================================================================== #

library(circlize)          # Color mapping for heatmaps
library(ComplexHeatmap)    # Advanced heatmap visualization
library(dplyr)             # Data manipulation
library(forcats)           # Factor handling
library(ggpmisc)           # Additional ggplot2 functionality
library(ggplot2)           # Data visualization
library(ggthemes)          # Extra themes for ggplot2
library(grid)              # Grid graphics
library(lubridate)         # Date/time handling
library(MASS)              # Robust statistical methods
library(purrr)             # Functional programming tools
library(readxl)            # Excel file import
library(reshape2)          # Data reshaping
library(scales)            # Scale functions for visualization
library(stringr)           # String manipulation
library(tidyr)             # Data tidying
library(viridis)           # Color scales for plots
library(vcd)               # Categorical data visualization
library(cvms)              # Confusion matrix visualization
library(psych)             # Correlation analysis
library(cowplot)           # Plot composition and alignment
library(Hmisc)             # Statistical tools
library(DataVisualizations) # Pareto Density Estimation

# ========================================================================== #
# 2. GLOBAL OPTIONS AND ANALYSIS SWITCHES
# ========================================================================== #

# Control analysis behavior
remove_censored <- FALSE              # Whether to exclude censored values
scale_0_100 <- FALSE                  # Whether to scale all measures to 0-100
analyze_only_untransformed <- FALSE   # Skip transformation analysis
plot_only_untransformed <- TRUE       # Show only original data in plots

# ========================================================================== #
# 3. HELPER FUNCTIONS
# ========================================================================== #

#' Sign-preserving logarithmic transformation (zero-invariant)
#'
#' @param x Numeric vector to transform
#' @param base Logarithm base (0 = natural log, 2, 10, or custom)
#' @return Transformed vector maintaining sign of original values
#' @details Useful for handling skewed distributions including zero and
#'          negative values. Preserves the sign of the original data.
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

#' Scale vector to specified range
#'
#' @param x Numeric vector to scale
#' @param minX Desired minimum value
#' @param maxX Desired maximum value
#' @return Scaled vector in range [minX, maxX]
scaleRange <- function(x, minX, maxX) {
  x_new <- (maxX - minX) * (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + minX
  return(x_new)
}

#' Min-max normalization to [0,1] with fixed boundaries
#'
#' @param x Numeric vector to normalize
#' @param minX Fixed minimum boundary
#' @param maxX Fixed maximum boundary
#' @return Normalized vector
scale01minmax <- function(x, minX, maxX) {
  (x - minX) / (maxX - minX)
}

#' Create Pareto Density Estimation plot with ggplot2
#'
#' @param Data Matrix or data frame with variables in columns
#' @return ggplot object showing PDE curves
#' @details Uses DataVisualizations package for robust density estimation
PDEplotGG <- function(Data) {
  Data <- as.matrix(Data)
  m <- matrix(NA, nrow = 0, ncol = 3)

  require(DataVisualizations)

  # Calculate PDE for each column
  for (i in 1:ncol(Data)) {
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
# 4. DATA IMPORT AND PREPARATION
# ========================================================================== #

# Import Excel data
cat("Loading data...\n")
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Extract and rename relevant trigeminal measures
trigeminal_measures_vars <- c("R28", "Lateralisierung (x/20)", "CO2-Schwelle")
trigeminal_measures_data <- trigeminale_daten_table1[, trigeminal_measures_vars]
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
    # Create segment identifier: rows 1-549 vs 550+
    # (different CO2 measurement protocols)
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

#' Custom ggplot2 theme for publication-quality plots
#'
#' @return theme object with consistent styling
#' @details Minimal theme with custom colors, fonts, and spacing
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
  scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Value [%]") +
  theme_plot() +
  labs(
    x = "Observation (Subject #)",
    y = "Measure",
    title = "Trigeminal measures",
    fill = "Value [%]"
  ) +
  # Vertical line separating different measurement protocols
  geom_vline(xintercept = 549.5, linetype = "dashed", color = "black", linewidth = 1) +
  # Annotation for first segment (no breath hold)
  annotate(
    "rect", xmin = 0, xmax = 548.5, ymin = 1.8, ymax = 2.2,
    fill = "white", alpha = 0.7
  ) +
  annotate(
    "text", x = 275, y = 2,
    label = "Breath not hold during CO2 threshold measurement",
    size = 5, color = "black"
  ) +
  # Annotation for second segment (breath hold)
  annotate(
    "rect", xmin = 550.5, xmax = 1001, ymin = 1.8, ymax = 2.2,
    fill = "white", alpha = 0.7
  ) +
  annotate(
    "text", x = 775, y = 2,
    label = "Breath hold during CO2 threshold measurement",
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
censored <- sapply(
  seq_along(max_vals),
  function(i) sum(trigeminal_measures_data[[i]] == max_vals[i], na.rm = TRUE)
)
percentage_censored <- (censored / counts) * 100
cat("\nOverall censoring:\n")
print(percentage_censored)

# Censoring in first subset (rows 1-549, breath not controlled)
subset1 <- trigeminal_measures_data[1:549, ]
counts1 <- apply(subset1, 2, function(x) sum(!is.na(x)))
censored1 <- sapply(
  seq_along(max_vals),
  function(i) sum(subset1[[i]] == max_vals[i], na.rm = TRUE)
)
percentage_censored1 <- (censored1 / counts1) * 100
cat("\nCensoring in subset 1 (breath not controlled):\n")
print(percentage_censored1)

# Censoring in second subset (rows 549+, breath hold protocol)
subset2 <- trigeminal_measures_data[549:nrow(trigeminal_measures_data), ]
counts2 <- apply(subset2, 2, function(x) sum(!is.na(x)))
censored2 <- sapply(
  seq_along(max_vals),
  function(i) sum(subset2[[i]] == max_vals[i], na.rm = TRUE)
)
percentage_censored2 <- (censored2 / counts2) * 100
cat("\nCensoring in subset 2 (breath hold):\n")
print(percentage_censored2)

# ========================================================================== #
# 7. GROUP AGREEMENT ANALYSIS - TOP 10% SENSITIVITY
# ========================================================================== #

cat("\nAnalyzing agreement among most sensitive subjects...\n")

# Create binary grouping: most sensitive 10% vs others
trigeminal_measures_data_scaled_grouped <- trigeminal_measures_data_scaled

# For AmmoLa and Lateralization: top 10% = most sensitive
trigeminal_measures_data_scaled_grouped[, 1:2] <- apply(
  trigeminal_measures_data_scaled_grouped[, 1:2],
  2,
  function(x) ifelse(x >= 90, 1, 0)
)

# For CO2: lower threshold = more sensitive (inverted scale)
trigeminal_measures_data_scaled_grouped$CO2_threshold <- ifelse(
  trigeminal_measures_data_scaled_grouped$CO2_threshold <= 10, 1, 0
)

# Summary of sensitivity thresholding
cat("\nDistribution of most sensitive subjects:\n")
print(apply(trigeminal_measures_data_scaled_grouped[, 1:3], 2, table))

# --- Agreement analysis for ALL subjects ---
trigeminal_measures_data_scaled_grouped_all <- trigeminal_measures_data_scaled_grouped

# Fisher's exact test: AmmoLa vs Lateralization
cat("\nFisher's test: AmmoLa vs Lateralization (all subjects)\n")
print(fisher.test(
  trigeminal_measures_data_scaled_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all$Lateralization
))
print(table(
  trigeminal_measures_data_scaled_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all$Lateralization
))

# Fisher's exact test: AmmoLa vs CO2
cat("\nFisher's test: AmmoLa vs CO2 (all subjects)\n")
print(fisher.test(
  trigeminal_measures_data_scaled_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all$CO2_threshold
))
print(table(
  trigeminal_measures_data_scaled_grouped_all$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all$CO2_threshold
))

# --- Agreement analysis ONLY for breath hold subjects ---
trigeminal_measures_data_scaled_grouped_all_breath_hold <-
  trigeminal_measures_data_scaled_grouped

# Exclude CO2 data from breath-not-controlled group
trigeminal_measures_data_scaled_grouped_all_breath_hold$CO2_threshold[1:549] <- NA

# Fisher's tests for breath hold group only
cat("\nFisher's test: AmmoLa vs Lateralization (breath hold only)\n")
print(fisher.test(table(
  trigeminal_measures_data_scaled_grouped_all_breath_hold$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all_breath_hold$Lateralization
)))

cat("\nFisher's test: AmmoLa vs CO2 (breath hold only)\n")
print(fisher.test(table(
  trigeminal_measures_data_scaled_grouped_all_breath_hold$AmmoLa_intensity,
  trigeminal_measures_data_scaled_grouped_all_breath_hold$CO2_threshold
)))

# Store tables and test results for visualization
table_AmmoLa_vs_CO2_most_sensitive <- table(
  trigeminal_measures_data_scaled_grouped_all_breath_hold[, c(1, 3)]
)
ftest_AmmoLa_vs_CO2_most_sensitive <- fisher.test(table_AmmoLa_vs_CO2_most_sensitive)
cat("\nAmmoLa vs CO2 agreement (breath hold):\n")
print(ftest_AmmoLa_vs_CO2_most_sensitive)

table_AmmoLa_vs_lateralization_most_sensitive <- table(
  trigeminal_measures_data_scaled_grouped_all_breath_hold[, c(1, 2)]
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
  AmmoLa = trigeminal_measures_data_scaled_grouped_all_breath_hold$AmmoLa_intensity,
  CO2 = trigeminal_measures_data_scaled_grouped_all_breath_hold$CO2_threshold
)
df_cm_cvms <- df_cm_cvms[complete.cases(df_cm_cvms), ]

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
  geom_text(aes(label = Freq), color = "white", size = 6) +
  scale_fill_gradient(low = "#deebf7", high = "#08519c") +
  theme_minimal() +
  labs(
    title = "Agreement of most sensitive subjects across tests",
    x = "CO2_threshold",
    y = "AmmoLa_intensity",
    fill = "Count"
  ) +
  # Semi-transparent box for statistics
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
  geom_text(aes(label = Freq), color = "white", size = 6) +
  scale_fill_gradient(low = "#deebf7", high = "#08519c") +
  theme_minimal() +
  labs(
    title = "Agreement of most sensitive subjects across tests",
    x = "Lateralization",
    y = "AmmoLa_intensity",
    fill = "Count"
  ) +
  # Semi-transparent box for statistics
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
# 9. VARIABLE TRANSFORMATION ANALYSIS (OPTIONAL)
# ========================================================================== #

if (!analyze_only_untransformed) {
  cat("\nApplying variable transformations...\n")

  # Define transformation functions
  reflect_log <- function(x) slog(max(x, na.rm = TRUE) + 1 - x)
  reflect_log_unflipped <- function(x) -slog(max(x, na.rm = TRUE) + 1 - x)
  square <- function(x) x^2
  log_transform <- function(x) slog(x)

  # Apply transformations to address skewness
  trigeminal_measures_data <- trigeminal_measures_data %>%
    mutate(
      AmmoLa_intensity_reflect_slog = reflect_log_unflipped(AmmoLa_intensity),
      Lateralization_square = square(Lateralization),
      CO2_threshold_slog = log_transform(CO2_threshold)
    )
}

# ========================================================================== #
# 10. AGE AND SEX EFFECTS ANALYSIS
# ========================================================================== #

cat("\nAnalyzing age correlations...\n")

# Add age to dataset and calculate correlations
trigeminal_measures_data_age <- trigeminal_measures_data
trigeminal_measures_data_age$age <- as.numeric(trigeminale_daten_table1$Alter)
corr_mat_age <- Hmisc::rcorr(
  as.matrix(trigeminal_measures_data_age),
  type = "pearson"
)

cat("\nAnalyzing sex differences...\n")

# Add sex to dataset
trigeminal_measures_data_sex <- trigeminal_measures_data
trigeminal_measures_data_sex$sex <- as.factor(trigeminale_daten_table1$Geschlecht)

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
  ss_total <- sum((x - mean(x))^2)
  ss_between <- sum(tapply(x, s, function(g) length(g) * (mean(g) - mean(x))^2))
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
    -grep("slog|square", trigeminal_measures_data_long$variable),
  ]
}

# --- Non-CO2 measures ---
trigeminal_measures_data_long_1 <- trigeminal_measures_data_long %>%
  filter(!str_detect(variable, "CO2"))

p_distribution_tigeminal_nonCO2 <- ggplot(
  trigeminal_measures_data_long_1,
  aes(x = value)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = "lightblue", color = "white", alpha = 0.7
  ) +
  geom_density(color = "darkblue", linewidth = 1) +
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

print(p_distribution_tigeminal_nonCO2)
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
    -grep("slog|square", trigeminal_measures_data_CO2_long$variable),
  ]
}

p_distribution_tigeminal_CO2 <- ggplot(
  trigeminal_measures_data_CO2_long,
  aes(x = value)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    fill = "lightblue", color = "white", alpha = 0.7
  ) +
  geom_density(color = "darkblue", size = 1) +
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

print(p_distribution_tigeminal_CO2)
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
  labs(title = "Distribution of AmmoLa intensity estimates") +
  guides(color = "none") +
  scale_color_manual(values = "dodgerblue")

# PDE plot for Lateralization
pPDE_Lateralization <- PDEplotGG(trigeminal_measures_data$Lateralization) +
  theme_plot() +
  labs(title = "Distribution of the number of correct lateralizations") +
  guides(color = "none") +
  scale_color_manual(values = "dodgerblue")

# PDE plot for CO2 thresholds (split by breath protocol)
CO2_threshold_breathing <- trigeminal_measures_data$CO2_threshold[1:549]
CO2_threshold_breath_hold <- trigeminal_measures_data$CO2_threshold[
  550:length(trigeminal_measures_data$CO2_threshold)
]

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
  labs(title = "Distribution of CO2 thresholds", color = "Breath") +
  scale_color_manual(
    values = c("grey33", "dodgerblue"),
    labels = c("uncontrolled", "hold")
  )

# ========================================================================== #
# 13. CORRELATION MATRIX PREPARATION
# ========================================================================== #

cat("\nPreparing correlation matrix...\n")

# Handle transformed variables if enabled
if (!analyze_only_untransformed) {
  trigeminal_measures_data_transformed <- trigeminal_measures_data[
    , c("AmmoLa_intensity_reflect_slog", "Lateralization_square", "CO2_threshold_slog")
  ]

  # Invert CO2 so higher = more sensitive
  trigeminal_measures_data_transformed$CO2_threshold_slog <-
    -trigeminal_measures_data_transformed$CO2_threshold_slog

  # Remove breath-not-controlled CO2 data
  trigeminal_measures_data_transformed <- trigeminal_measures_data_transformed %>%
    mutate(CO2_threshold_slog = ifelse(row_number() <= 549, NA, CO2_threshold_slog))

  # Define variable pairs for correlation analysis
  var_pairs <- list(
    c("AmmoLa_intensity_reflect_slog", "Lateralization_square"),
    c("AmmoLa_intensity_reflect_slog", "CO2_threshold_slog")
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
all_measures_correlations$CO2_threshold_slog <- -all_measures_correlations$CO2_threshold_slog

# Split CO2 by breath control protocol
all_measures_correlations$CO2_threshold_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_not_hold[1:549] <-
  all_measures_correlations$CO2_threshold[1:549]

all_measures_correlations$CO2_threshold_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_hold[
  550:nrow(all_measures_correlations)
] <- all_measures_correlations$CO2_threshold[550:nrow(all_measures_correlations)]

all_measures_correlations$CO2_threshold_slog_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_not_hold[1:549] <-
  all_measures_correlations$CO2_threshold_slog[1:549]

all_measures_correlations$CO2_threshold_slog_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_hold[
  550:nrow(all_measures_correlations)
] <- all_measures_correlations$CO2_threshold_slog[550:nrow(all_measures_correlations)]

# Filter variables based on analysis settings
if (analyze_only_untransformed | plot_only_untransformed) {
  all_measures_correlations <- all_measures_correlations[
    , -c(grep("slog|square", names(all_measures_correlations)))
  ]
}

# Remove redundant columns
all_measures_correlations <- all_measures_correlations[
  , !names(all_measures_correlations) %in% c("CO2_threshold_breath_not_hold", "CO2_threshold")
]

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
corr_mat <- cor_results$r  # Correlation coefficients
p_mat <- cor_results$p      # P-values

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

# Define color palette for correlation heatmap (NYT-inspired)
breaks <- c(
  -1, -0.9, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.02, -0.01, 0,
  0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 1
)

nyt_colors <- c(
  "ghostwhite", "#fbfbfb", "#e6f0fa", "#c9def9",
  "#add0fa", "#7bb8fa", "dodgerblue2", "#041a58"
)
nyt_colors <- c(rev(nyt_colors), nyt_colors)
color_vec <- colorRampPalette(nyt_colors)(length(breaks))
col_fun <- colorRamp2(breaks, color_vec)

#' Determine text color based on background brightness
#'
#' @param fill_color Background color
#' @return Text color (dark or light) for optimal contrast
text_color_fun <- function(fill_color) {
  rgb_val <- col2rgb(fill_color) / 255
  # Calculate luminance using standard weights
  brightness <- 0.299 * rgb_val[1, ] + 0.587 * rgb_val[2, ] + 0.114 * rgb_val[3, ]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}

#' Create correlation heatmap with ComplexHeatmap
#'
#' @return NULL (draws plot directly)
create_heatmap_trig <- function() {
  ht <- Heatmap(
    corr_mat,
    name = "Correlation",
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE,
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    cluster_columns = TRUE,
    clustering_method_columns = "ward.D2",
    show_column_dend = FALSE,
    column_names_rot = 45,
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(4, "cm"),
    # Custom cell function to add correlation values and significance stars
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
  labs(title = "Correlation matrix (spearman)") +
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

# ========================================================================== #
# 15. COMBINED PUBLICATION FIGURE
# ========================================================================== #

cat("\nCreating combined publication figure...\n")

# Combine all main plots into publication-ready figure
combined_trigeminal_analysis_plot <- cowplot::plot_grid(
  # Top row: Distribution plots (A, B, C)
  cowplot::plot_grid(
    pPDE_AmmoLA,
    pPDE_Lateralization,
    pPDE_CO2,
    labels = LETTERS[1:3],
    nrow = 1,
    align = "h",
    axis = "tb",
    label_y = 0.97
  ),
  # Bottom row: Agreement tables (D, E) and correlation matrix (F)
  cowplot::plot_grid(
    cowplot::plot_grid(
      p_table_AmmoLa_vs_CO2_most_sensitive,
      p_table_AmmoLa_vs_lateralization_most_sensitive,
      labels = LETTERS[4:5],
      nrow = 1,
      align = "h",
      axis = "tb",
    label_y = 0.97
    ),
    cowplot::plot_grid(
      corr_plot,
      labels = LETTERS[6],
      nrow = 1,
      align = "h",
      axis = "tb",
      label_y = 0.97
    ),
    nrow = 1,
    rel_widths = c(2, 1),
    align = "h",
    axis = "tb"
  ),
  nrow = 2,
  rel_heights = c(1, 1),
  align = "v",
  axis = "lr"
)

# Display combined figure
print(combined_trigeminal_analysis_plot)

# Save combined figure
ggsave(
  "combined_trigeminal_analysis_plot.svg",
  combined_trigeminal_analysis_plot,
  width = 20,
  height = 12
)

# ========================================================================== #
# END OF SCRIPT
# ========================================================================== #

cat("\nAnalysis complete!\n")
