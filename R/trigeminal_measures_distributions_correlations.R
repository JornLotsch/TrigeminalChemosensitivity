################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: [YOUR NAME]
# Date: [DATE]
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, transformation, correlation and visualization.
################################################################################

# ======================== #
# 1. Load Required Libraries
# ======================== #
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(forcats)
library(ggpmisc)
library(ggplot2)
library(ggthemes)
library(grid)
library(lubridate)
library(MASS) # geom_smooth (robust fits)
library(purrr)
library(readxl)
library(reshape2)
library(scales)
library(stringr)
library(tidyr)
library(viridis)
library(vcd)
library(cvms)
library(psych)
library(cowplot)

# ================================= #
# 2. Global Options and Switches
# ================================= #
remove_censored <- FALSE
scale_0_100 <- FALSE
analyze_only_untransformed <- FALSE
plot_only_untransformed <- TRUE

# =============================== #
# 3. Helper Functions
# =============================== #

# Zero-invariant log transformation
slog <- function(x, base = 10) {
  absX <- abs(x)
  s <- sign(x)
  if (base == 0) return(s * log1p(absX))
  else if (base == 2) return(s * log2(absX + 1))
  else if (base == 10) return(s * log10(absX + 1))
  else return(s * log1p(absX) * log(base))
}

# Scale to given range [minX, maxX]
scaleRange <- function(x, minX, maxX) {
  x_new <- (maxX - minX) * (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + minX
  return(x_new)
}

# Min-max scaling to 0-1 given a fixed [minX, maxX]
scale01minmax <- function(x, minX, maxX) {
  (x - minX) / (maxX - minX)
}

# =============================== #
# 4. Data Import & Preparation
# =============================== #

# Import Excel data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Extract relevant variables
trigeminal_measures_vars <- c("R28", "Lateralisierung (x/20)", "CO2-Schwelle")
trigeminal_measures_data <- trigeminale_daten_table1[, trigeminal_measures_vars]
names(trigeminal_measures_data) <- c("AmmoLa_intensity", "Lateralization", "CO2_threshold")

# ============================= #
# 5. Initial Data Visualization
# ============================= #

## Scale measures for heatmap and create observation segment
trigeminal_measures_data_scaled <- trigeminal_measures_data %>%
  mutate(
    Lateralization = scale01minmax(Lateralization, minX = 0, maxX = 20) * 100,
    CO2_threshold = scale01minmax(CO2_threshold, minX = 100, maxX = 2000) * 100,
    Segment = if_else(row_number() <= 549, "first_part", "second_part")
  )

## Reshape for heatmap plotting
heatmap_data <- trigeminal_measures_data_scaled %>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols = c(AmmoLa_intensity, Lateralization, CO2_threshold),
    names_to = "Measure",
    values_to = "Value"
  )

## Custom ggplot2 theme for consistency
theme_plot <- function() {
  theme_minimal(base_family = "Libre Franklin") +
    theme(
      plot.title = element_text(face = "plain", size = 18, color = "#222222", hjust = 0, margin = margin(b = 10)),
      axis.title = element_text(face = "plain", size = 12, color = "#444444"),
      axis.text = element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = element_text(size = 8, color = "#888888", hjust = 0, margin = margin(t = 10)),
      panel.grid.major.y = element_line(color = "#dddddd", linetype = "dashed", size = 0.3),
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

## Make and save heatmap plot
p_trigeminal_measures_done <- ggplot(heatmap_data, aes(x = Row, y = Measure, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Value [%]") +
  theme_plot() +
  labs(
    x = "Observation (Subject #)",
    y = "Measure",
    title = "Trigeminal measures",
    fill = "Value [%]"
  ) +
  geom_vline(xintercept = 549.5, linetype = "dashed", color = "black", linewidth = 1) +
  annotate("rect", xmin = 60, xmax = 480, ymin = 1.8, ymax = 2.2, fill = "white", alpha = 0.7) +
  annotate("text", x = 275, y = 2, label = "Breath not hold during CO2 threshold measurement", size = 5, color = "black") +
  annotate("rect", xmin = 575, xmax = 970, ymin = 1.8, ymax = 2.2, fill = "white", alpha = 0.7) +
  annotate("text", x = 775, y = 2, label = "Breath hold during CO2 threshold measurement", size = 5, color = "black")
print(p_trigeminal_measures_done)
ggsave("p_trigeminal_measures_done.svg", p_trigeminal_measures_done, width = 16, height = 4)

# ========================= #
# 6. Data Censoring Checks
# ========================= #
max_vals <- c(100, 20, 2000)
counts <- apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x)))
censored <- sapply(seq_along(max_vals), function(i) sum(trigeminal_measures_data[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored <- (censored / counts) * 100
print(percentage_censored)

subset1 <- trigeminal_measures_data[1:549,]
counts1 <- apply(subset1, 2, function(x) sum(!is.na(x)))
censored1 <- sapply(seq_along(max_vals), function(i) sum(subset1[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored1 <- (censored1 / counts1) * 100
print(percentage_censored1)

subset2 <- trigeminal_measures_data[549:nrow(trigeminal_measures_data),]
counts2 <- apply(subset2, 2, function(x) sum(!is.na(x)))
censored2 <- sapply(seq_along(max_vals), function(i) sum(subset2[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored2 <- (censored2 / counts2) * 100
print(percentage_censored2)

# ====================== #
# 7. Group Agreement Analyses
# ====================== #
trigeminal_measures_data_scaled[, 1:2] <- apply(trigeminal_measures_data_scaled[, 1:2], 2, function(x) ifelse(x >= 90, 1, 0))
trigeminal_measures_data_scaled$CO2_threshold <- ifelse(trigeminal_measures_data_scaled$CO2_threshold <= 10, 1, 0)
apply(trigeminal_measures_data_scaled[, 1:3], 2, table) # Sensitivity thresholding summary

trigeminal_measures_data_scaled1 <- trigeminal_measures_data_scaled

# Fisher's exact tests
fisher.test(trigeminal_measures_data_scaled$AmmoLa_intensity, trigeminal_measures_data_scaled$Lateralization)
fisher.test(trigeminal_measures_data_scaled$AmmoLa_intensity, trigeminal_measures_data_scaled$CO2_threshold)

## Analyze only rows with breath hold
trigeminal_measures_data_scaled <- trigeminal_measures_data_scaled[549:nrow(trigeminal_measures_data),]
fisher.test(trigeminal_measures_data_scaled$AmmoLa_intensity, trigeminal_measures_data_scaled$Lateralization)
fisher.test(trigeminal_measures_data_scaled$AmmoLa_intensity, trigeminal_measures_data_scaled$CO2_threshold)
table2 <- table(trigeminal_measures_data_scaled[, c(1, 3)])

# =========================== #
# 8. Confusion Matrix Analysis
# =========================== #
df_cm_cvms <- cbind.data.frame(
  AmmoLa = trigeminal_measures_data_scaled$AmmoLa_intensity,
  CO2 = trigeminal_measures_data_scaled$CO2_threshold
)
df_cm_cvms <- df_cm_cvms[complete.cases(df_cm_cvms),]

cm <- confusion_matrix(targets = df_cm_cvms$CO2, predictions = df_cm_cvms$AmmoLa)
cm_stats <- cm[1, c("Balanced Accuracy", "F1", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Kappa", "MCC", "Detection Rate")]
vec <- unlist(cm_stats)
formatted <- paste(names(vec), signif(vec, 3), sep = ": ")
midpoint <- ceiling(length(formatted) / 2)
stat_string <- paste(paste(formatted[1:midpoint], collapse = " | "), paste(formatted[(midpoint + 1):length(formatted)], collapse = " | "), sep = "\n")

tbl <- table(df_cm_cvms$AmmoLa, df_cm_cvms$CO2)
fisher_result <- fisher.test(tbl)
fisher_p <- signif(fisher_result$p.value, 3)

p_cm_AmmoLa_vs_CO2 <- plot_confusion_matrix(cm$`Confusion Matrix`[[1]], add_sums = FALSE) +
  ggplot2::labs(
    subtitle = stat_string,
    caption = paste("Fisher's exact test p-value:", fisher_p)
  ) +
  ggplot2::xlab("Predictions: AmmoLa") +
  ggplot2::ylab("Target: CO2") +
  theme_plot()
print(p_cm_AmmoLa_vs_CO2)
ggsave("p_cm_AmmoLa_vs_CO2.svg", p_cm_AmmoLa_vs_CO2, width = 8, height = 8)


# ======================== #
# 9. Optional Scaling and Data Export
# ======================== #
if (scale_0_100) {
  trigeminal_measures_data$Lateralization <- scale01minmax(trigeminal_measures_data$Lateralization, minX = 0, maxX = 20) * 100
  trigeminal_measures_data$CO2_threshold <- scale01minmax(trigeminal_measures_data$CO2_threshold, minX = 0, maxX = 2000) * 100
}
if (remove_censored) {
  for (i in seq_along(max_vals)) {
    trigeminal_measures_data[[i]][trigeminal_measures_data[[i]] == max_vals[i]] <- NA
  }
}
if (scale_0_100) {
  write.csv(trigeminal_measures_data, "trigeminal_measures_scaled_0_100_data.csv")
} else {
  write.csv(trigeminal_measures_data, "trigeminal_measures_data.csv")
}

counts_non_NA <- apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x)))
counts_valid <- apply(trigeminal_measures_data[549:nrow(trigeminal_measures_data),], 2, function(x) sum(!is.na(x)))
counts_non_NA
counts_valid
counts_non_NA - counts_valid

# ================================ #
# 10. Variable Transformations
# ================================ #
if (!analyze_only_untransformed) {
  reflect_log <- function(x) slog(max(x, na.rm = TRUE) + 1 - x)
  reflect_log_unflipped <- function(x) - slog(max(x, na.rm = TRUE) + 1 - x)
  square <- function(x) x ^ 2
  log_transform <- function(x) slog(x)
  trigeminal_measures_data <- trigeminal_measures_data %>%
    mutate(
      AmmoLa_intensity_reflect_slog = reflect_log_unflipped(AmmoLa_intensity),
      Lateralization_square = square(Lateralization),
      CO2_threshold_slog = log_transform(CO2_threshold)
    )
}

# =============================================== #
# 11. Check Correlation with Age and Sex Effects
# =============================================== #
trigeminal_measures_data_age <- trigeminal_measures_data
trigeminal_measures_data_age$age <- as.numeric(trigeminale_daten_table1$Alter)
corr_mat_age <- cor(trigeminal_measures_data_age, use = "pairwise.complete.obs", method = "pearson")

trigeminal_measures_data_sex <- trigeminal_measures_data
trigeminal_measures_data_sex$sex <- as.factor(trigeminale_daten_table1$Geschlecht)
sex_diff_trig <- apply(
  trigeminal_measures_data_sex[, 1:(ncol(trigeminal_measures_data_sex) - 1)], 2,
  function(x) kruskal.test(x ~ as.factor(trigeminal_measures_data_sex$sex))
)
variables <- colnames(trigeminal_measures_data_sex)[colnames(trigeminal_measures_data_sex) != "sex"]
sex_stats <- lapply(variables, function(var) {
  x <- trigeminal_measures_data_sex[[var]]
  s <- trigeminal_measures_data_sex$sex
  idx <- complete.cases(x, s)
  x <- x[idx]
  s <- s[idx]
  aovmod <- aov(x ~ s)
  ss_total <- sum((x - mean(x)) ^ 2)
  ss_between <- sum(tapply(x, s, function(g) length(g) * (mean(g) - mean(x)) ^ 2))
  eta2 <- ss_between / ss_total
  kruskal_p <- tryCatch(kruskal.test(x ~ s)$p.value, error = function(e) NA)
  n <- length(x)
  data.frame(variable = var, eta2 = eta2, kruskal_p = kruskal_p, n = n)
})
sex_stats_df <- bind_rows(sex_stats) %>%
  arrange(desc(eta2))
print(sex_stats_df)

# =============================== #
# 12. Distribution Plots
# =============================== #
trigeminal_measures_data_long <- reshape2::melt(trigeminal_measures_data)
if (plot_only_untransformed) {
  trigeminal_measures_data_long <- trigeminal_measures_data_long[-grep("slog|square", trigeminal_measures_data_long$variable),]
}
## Plot non-CO2 measures
trigeminal_measures_data_long_1 <- trigeminal_measures_data_long %>%
  filter(!str_detect(variable, "CO2"))
p_distribution_tigeminal_nonCO2 <- ggplot(trigeminal_measures_data_long_1, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Distribution of non-CO2 related trigeminal measures", x = "Value", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_plot()
print(p_distribution_tigeminal_nonCO2)
ggsave("p_distribution_tigeminal_nonCO2.svg", p_distribution_tigeminal_nonCO2, width = 10, height = ifelse(plot_only_untransformed, 5, 10))

# Prepare CO2 data/labels
trigeminal_measures_data_CO2 <- trigeminal_measures_data[, str_detect(names(trigeminal_measures_data), "CO2")]
trigeminal_measures_data_CO2 <- trigeminal_measures_data_CO2 %>% mutate(row_id = row_number())
trigeminal_measures_data_CO2 <- trigeminal_measures_data_CO2 %>% mutate(subset_group = ifelse(row_id <= 549, "Normal breath", "Breath hold"))
trigeminal_measures_data_CO2_long <- trigeminal_measures_data_CO2 %>%
  pivot_longer(cols = -c(row_id, subset_group), names_to = "variable", values_to = "value")
if (plot_only_untransformed) trigeminal_measures_data_CO2_long <- trigeminal_measures_data_CO2_long[-grep("slog|square", trigeminal_measures_data_CO2_long$variable),]

# Plot CO2 measures by group
p_distribution_tigeminal_CO2 <- ggplot(trigeminal_measures_data_CO2_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), fill = "lightblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(subset_group * variable ~ ., scales = "free", ncol = 2) +
  labs(title = "Distribution of CO2-related trigeminal measures, split by procedure", x = "Value", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0), strip.text = element_text(face = "bold", size = 10), panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_plot()
print(p_distribution_tigeminal_CO2)
ggsave("p_distribution_tigeminal_CO2.svg", p_distribution_tigeminal_CO2, width = 10, height = ifelse(plot_only_untransformed, 5, 10))

# =============================== #
# 13. Correlation/Heatmap Analysis
# =============================== #
if (!analyze_only_untransformed) {
  trigeminal_measures_data_transformed <- trigeminal_measures_data[, c("AmmoLa_intensity_reflect_slog", "Lateralization_square", "CO2_threshold_slog")]
  trigeminal_measures_data_transformed$CO2_threshold_slog <- -trigeminal_measures_data_transformed$CO2_threshold_slog
  trigeminal_measures_data_transformed <- trigeminal_measures_data_transformed %>%
    mutate(CO2_threshold_slog = ifelse(row_number() <= 549, NA, CO2_threshold_slog))
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
  # Inline functions, stat labels, plotting code as before...
}

# Prepare data for correlation heatmap
all_measures_correlations <- trigeminal_measures_data
all_measures_correlations$CO2_threshold <- -all_measures_correlations$CO2_threshold
all_measures_correlations$CO2_threshold_slog <- -all_measures_correlations$CO2_threshold_slog
all_measures_correlations$CO2_threshold_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold[1:549]
all_measures_correlations$CO2_threshold_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold[550:nrow(all_measures_correlations)]
all_measures_correlations$CO2_threshold_slog_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold_slog[1:549]
all_measures_correlations$CO2_threshold_slog_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold_slog[550:nrow(all_measures_correlations)]

# Filter data
if (analyze_only_untransformed | plot_only_untransformed) all_measures_correlations <- all_measures_correlations[, - c(grep("slog|square", names(all_measures_correlations)))]
all_measures_correlations <- all_measures_correlations[, !names(all_measures_correlations) %in% c("CO2_threshold_breath_not_hold", "CO2_threshold")]

################################################################################
# 14. Correlation Matrix Calculation and Heatmap Visualization
################################################################################

# --- Select only numeric columns for correlation analysis ---
numeric_data <- all_measures_correlations %>%
  dplyr::select(where(is.numeric))

# --- Calculate Spearman correlations; get correlation and p-value matrices ---
cor_method <- "spearman"
cor_results <- corr.test(numeric_data, use = "pairwise", method = cor_method)
corr_mat <- cor_results$r # correlation coefficients
p_mat <- cor_results$p # p-values

# --- Significance code for stars ---
signif_code <- function(p) {
  if (is.na(p)) ""
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}

# --- Define NYT-inspired color palette for correlation heatmap ---
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

# --- Adaptive text color function for heatmaps ---
text_color_fun <- function(fill_color) {
  rgb_val <- col2rgb(fill_color) / 255
  brightness <- 0.299 * rgb_val[1,] + 0.587 * rgb_val[2,] + 0.114 * rgb_val[3,]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}

# --- Plot: correlation heatmap with colored values ---

# Create heatmap without title (remove any column_title argument)
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

# Grab the heatmap plot as grob
gp <- grid.grabExpr(create_heatmap_trig())

# Create the title plot with minimal bottom margin
p_title <- ggplot() +
  labs(title = "Correlation matrix of sensitivity-scaled measures (low value = high sensitivity) (spearman)") +
  theme_minimal(base_family = "Libre Franklin") +
  theme(
    plot.title = element_text(face = "plain", size = 18, color = "#222222", hjust = 0, margin = margin(b = 2)),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 2) # reduce bottom margin
  ) +
  theme_void()


# Combine with ggplot title grob using cowplot
corr_plot <- plot_grid(p_title, gp, ncol = 1, rel_heights = c(0.1, 2), align = "v")

# Save as SVG or display
ggsave("trigeminal_correlation_heatmap.svg", corr_plot, width = 8, height = 9)
print(corr_plot)


# ---- END OF SCRIPT ----
