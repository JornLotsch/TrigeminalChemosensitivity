# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(forcats)
library(readxl)
library(scales)
library(reshape2)
library(ggpmisc)
library(MASS)       # robust linear model for geom_smooth
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggthemes)
library(viridis)
library(opGMMassessment)  # assuming needed for something not shown here

# ----------------------------
# Define zero-invariant log transform function
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
    factor <- log(base)
    return(s * log1p(absX) * factor)
  }
}

# ----------------------------
# Read Excel file with trigeminal sensitivity data
trigeminale_daten_table1 <- read_excel(
  "/home/joern/Aktuell/TrigeminalSensitivity/09Originale/Bormann Trigeminale Studie Daten.xlsx",
  sheet = "Tabelle1"
)

# Variables of interest related to facial pain
trigeminal_measures_vars <- c("R28", "Lateralisierung (x/20)", "CO2-Schwelle")

# Subset data frame and rename columns for clarity
trigeminal_measures_data <- trigeminale_daten_table1[, trigeminal_measures_vars]
names(trigeminal_measures_data) <- c("AmmoLa_intensity", "Lateralization", "CO2_threshold")

# Save raw trigeminal measures data as CSV for record keeping
write.csv(trigeminal_measures_data, "trigeminal_measures_data.csv")

# ----------------------------
# Define transformation functions
reflect_log <- function(x) slog(max(x, na.rm = TRUE) + 1 - x)
square <- function(x) x^2
log_transform <- function(x) slog(x)

# Apply transformations and add new columns
trigeminal_measures_data <- trigeminal_measures_data %>%
  mutate(
    AmmoLa_intensity_reflect_slog = reflect_log(AmmoLa_intensity),
    Lateralization_square = square(Lateralization),
    CO2_threshold_slog = log_transform(CO2_threshold)
  )

# ----------------------------
# Plot distributions of non-CO2 related trigeminal measures
trigeminal_measures_data_long <- reshape2::melt(trigeminal_measures_data)

# Filter to exclude CO2 variables for this plot
trigeminal_measures_data_long_1 <- trigeminal_measures_data_long %>%
  filter(!str_detect(variable, "CO2"))

ggplot(trigeminal_measures_data_long_1, aes(x = value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = "lightblue", color = "white", alpha = 0.7
  ) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(~ variable, scales = "free") +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# ----------------------------
# Prepare CO2-related data and subsets for plotting
trigeminal_measures_data_CO2 <- trigeminal_measures_data[,str_detect(names(trigeminal_measures_data), "CO2") ]

# Add row index to original CO2 data
trigeminal_measures_data_CO2 <- trigeminal_measures_data_CO2 %>%
  mutate(row_id = row_number())

# Create subset group label based on row number
trigeminal_measures_data_CO2 <- trigeminal_measures_data_CO2 %>%
  mutate(subset_group = ifelse(row_id <= 549, "Normal breath", "Breath hold"))

# Pivot to long format for plotting
trigeminal_measures_data_CO2_long <- trigeminal_measures_data_CO2 %>%
  pivot_longer(
    cols = -c(row_id, subset_group),
    names_to = "variable",
    values_to = "value"
  )

# Plot with facets by variable and subset_group
ggplot(trigeminal_measures_data_CO2_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(subset_group * variable ~., scales = "free", ncol = 2) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# ----------------------------
# Prepare transformed data for correlation analysis and plotting
trigeminal_measures_data_transformed <- trigeminal_measures_data[,c("AmmoLa_intensity_reflect_slog", "Lateralization_square", "CO2_threshold_slog")]
# Step 1: Replace CO2_threshold_slog values with NA for first 549 rows
trigeminal_measures_data_transformed <- trigeminal_measures_data_transformed %>%
  mutate(CO2_threshold_slog = ifelse(row_number() <= 549, NA, CO2_threshold_slog))

# Step 2: Define variable pairs for plotting
var_pairs <- list(
  c("AmmoLa_intensity_reflect_slog", "Lateralization_square"),
  c("AmmoLa_intensity_reflect_slog", "CO2_threshold_slog")
)

# Step 3: Create a combined data frame for long plotting of pairs
plot_list <- map(var_pairs, function(vars) {
  df <- trigeminal_measures_data_transformed[, vars]
  names(df) <- c("x", "y")
  df$pair <- paste(vars, collapse = " vs ")
  df
})

plot_df <- bind_rows(plot_list)

# Compute correlation stats for annotation
cor_stats <- map_dfr(var_pairs, function(vars) {
  x <- trigeminal_measures_data_transformed[[vars[1]]]
  y <- trigeminal_measures_data_transformed[[vars[2]]]
  cor_test <- cor.test(x, y, use = "complete.obs")

  data.frame(
    pair = paste(vars, collapse = " vs "),
    label = sprintf("r = %.2f, p = %.3g", cor_test$estimate, cor_test$p.value),
    # Position near top-left inside panel (numeric)
    x = min(x, na.rm = TRUE) + 0.05 * diff(range(x, na.rm = TRUE)),
    y = max(y, na.rm = TRUE) - 0.10 * diff(range(y, na.rm = TRUE))
  )
})

formula_rlm <- y ~ x
# Custom labeller to display y-variable name as strip label
custom_labeller <- function(variable, value) {
  sapply(value, function(x) {
    # Extract y variable from "x_var vs y_var"
    y_var <- strsplit(x, " vs ")[[1]][2]
    return(y_var)
  })
}


p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "rlm",
    formula = formula_rlm,
    se = TRUE,
    fullrange = FALSE,
    color = "darkblue",
    fill = "lightblue",
    alpha = 0.25
  ) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = formula_rlm,
    parse = TRUE,
    size = 4,
    label.x.npc = "left",
    label.y.npc = 0.15,
    color = "red"
  ) +
  geom_text(
    data = cor_stats,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 4.5,
    color = "blue"
  ) +
  facet_grid(y_var ~ x_var, scales = "free", labeller = labeller(pair = custom_labeller)) +
  labs(
    title = "Pairwise Relationships of Trigeminal Measures",
    x = "AmmoLa_intensity_reflect_slog",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    strip.text = element_text(face = "bold", size = 14),
    axis.title.y = element_text(size = 12),
    # Custom y-axis label per facet
    strip.background = element_blank()
  ) +
  # Custom function to set dynamic y-axis label per facet via ggplot2 vignettes workaround
  facet_wrap(~ pair, scales = "free", labeller = labeller(pair = function(x) {
    sapply(x, function(p) {
      # Extract y_var part for label
      y_var <- strsplit(p, " vs ")[[1]][2]
      return(y_var)
    })
  }))

print(p)


# ----------------------------
# Prepare data for correlation heatmap --------------------------------------

# Mask CO2_threshold and transformed CO2_threshold for first 549 rows (set to NA)
all_measurs_correlations <- trigeminal_measures_data %>%
  mutate(CO2_threshold_slog = ifelse(row_number() <= 549, NA, CO2_threshold_slog)) %>%
  mutate(CO2_threshold = ifelse(row_number() <= 549, NA, CO2_threshold))

# Calculate pairwise correlation matrix with pairwise complete observations
corr_mat <- cor(all_measurs_correlations, use = "pairwise.complete.obs", method = "pearson")

# Define NYT-inspired color palette for correlation heatmap
breaks <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 1)
nyt_colors <- c(
  "ghostwhite",
  "#fbfbfb",
  "#e6f0fa",
  "#c9def9",
  "#add0fa",
  "#7bb8fa",
  "dodgerblue2",
  "#041a58",
  "#041a58",
  "#041a58",
  "#041a58",
  "#03124a"
)
color_vec <- colorRampPalette(nyt_colors)(length(breaks))

col_fun <- colorRamp2(breaks, color_vec)

# Function to set text color (dark on light backgrounds, white on dark)
text_color_fun <- function(fill_color) {
  rgb_val <- col2rgb(fill_color) / 255
  brightness <- 0.299 * rgb_val[1,] + 0.587 * rgb_val[2,] + 0.114 * rgb_val[3,]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}

# Create correlation heatmap with values and adaptive text color
ht <- Heatmap(
  corr_mat,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  show_column_dend = FALSE,
  row_dend_width = unit(4, "cm"),
  column_dend_height = unit(4, "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- sprintf("%.2f", corr_mat[i, j])
    col_text <- text_color_fun(fill)
    grid.text(val, x, y, gp = gpar(fontsize = 10, col = col_text))
  },
  rect_gp = gpar(col = NA),
  border = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  top_annotation = NULL,
  show_heatmap_legend = TRUE
)

# Render heatmap on a new page
grid.newpage()

draw(
  ht,
  heatmap_legend_side = "bottom"
)

# Add a bold, left-aligned title manually above heatmap
grid.text(
  "Correlation matrix",
  x = unit(0, "npc") + unit(4, "mm"),   # left margin
  y = unit(1, "npc") - unit(4, "mm"),   # top margin
  just = c("left", "top"),
  gp = gpar(fontsize = 16, fontface = "bold")
)
