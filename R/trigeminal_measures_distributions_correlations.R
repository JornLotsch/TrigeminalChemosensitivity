# Load required libraries
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(forcats)
library(ggpmisc)
library(ggplot2)
library(ggthemes)
library(grid)
library(lubridate)
library(MASS)  # robust linear model for geom_smooth
library(purrr)
library(readxl)
library(reshape2)
library(scales)
library(stringr)
library(tidyr)
library(viridis)
library(grid)

# Switches
remove_censored <- TRUE
scale_0_100 <- FALSE

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

scaleRange <- function(x, minX, maxX) {
  x_new <- (maxX - minX) * (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + minX
  return(x_new)
}

scale01minmax <- function(x, minX, maxX) {
  x_new <- (x - minX) / (maxX - minX)
  return(x_new)
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

# First, plot the data as a heatmap
# Scale relevant columns
trigeminal_measures_data_scaled <- trigeminal_measures_data %>%
  mutate(
    Lateralization = scale01minmax(Lateralization, minX = 0, maxX = 20) * 100,
    CO2_threshold = scale01minmax(CO2_threshold, minX = 0, maxX = 2000) * 100,
    Segment = if_else(row_number() <= 549, "first_part", "second_part")
  )

heatmap_data <- trigeminal_measures_data_scaled %>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols = c(AmmoLa_intensity, Lateralization, CO2_threshold),
    names_to = "Measure",
    values_to = "Value"
  )

# Plot (horizontal heatmap)
# Custom theme for ggplot2
theme_plot <- function() {
  theme_minimal(base_family = "Libre Franklin") +  # Replace font if installed
    theme(
      # Text elements
      plot.title = element_text(face = "plain", size = 18, color = "#222222", hjust = 0, margin = margin(b = 10)),
      axis.title = element_text(face = "plain", size = 12, color = "#444444"),
      axis.text = element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = element_text(size = 8, color = "#888888", hjust = 0, margin = margin(t = 10)),

      # Grid lines
      panel.grid.major.y = element_line(color = "#dddddd", linetype = "dashed", size = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),

      # Axis lines and ticks
      axis.line = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks.length = unit(5, "pt"),

      # Plot background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      # Legend
      legend.position = "right",
      legend.direction = "vertical",

      # Margins and spacing for breathing room
      plot.margin = margin(20, 20, 20, 20),

      # Facet tweaks if needed
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "#222222")
    )
}

# Plot
p_trigeminal_measures_done <- ggplot(heatmap_data, aes(x = Row, y = Measure, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Value [%]") +
  theme_plot() +
  labs(
    x = "Observation (Row)",
    y = "Measure",
    title = "Trigeminal measures",
    fill = "Value [%]"
  ) +
  geom_vline(xintercept = 549.5, linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = 275, y = 0.5,
           label = "Breath not hold during CO2 threshold measurement",
           vjust = -1, size = 5, color = "black") +
  annotate("text", x = 775, y = 0.5,
           label = "Breath hold during CO2 threshold measurement",
           vjust = -1, size = 5, color = "black")

p_trigeminal_measures_done

ggsave(paste0("p_trigeminal_measures_done", ".svg"), p_trigeminal_measures_done, width = 16, height = 4)

# Check how many data is censored
max_vals <- c(100, 20, 2000)
counts <- apply(trigeminal_measures_data, 2, function(x) sum(!is.na(x)))
censored <- sapply(seq_along(max_vals), function(i) sum(trigeminal_measures_data[[i]] == max_vals[i], na.rm = TRUE))
percentage_censored <- (censored / counts) * 100
print(percentage_censored)

# And again, for the change in CO2 measuring
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

# Optionally scale all variables 0 - 100
if (scale_0_100) {
  trigeminal_measures_data$Lateralization <- scale01minmax(trigeminal_measures_data$Lateralization, minX = 0, maxX = 20) * 100
  trigeminal_measures_data$CO2_threshold <- scale01minmax(trigeminal_measures_data$CO2_threshold, minX = 0, maxX = 2000) * 100
}

# Optionally: Replace all censored data with NA
if (remove_censored) {
  for (i in seq_along(max_vals)) {
    trigeminal_measures_data[[i]][trigeminal_measures_data[[i]] == max_vals[i]] <- NA
  }
}

# Save raw trigeminal measures data as CSV for record keeping
if (scale_0_100) {
  write.csv(trigeminal_measures_data, "trigeminal_measures_scaled_0_100_data.csv")
} else {
  write.csv(trigeminal_measures_data, "trigeminal_measures_data.csv")
}

# ----------------------------
# Define transformation functions
reflect_log <- function(x) slog(max(x, na.rm = TRUE) + 1 - x)

reflect_log_unflipped <- function(x) - slog(max(x, na.rm = TRUE) + 1 - x)

square <- function(x) x ^ 2

log_transform <- function(x) slog(x)

# Apply transformations and add new columns
trigeminal_measures_data <- trigeminal_measures_data %>%
  mutate(
    AmmoLa_intensity_reflect_slog = reflect_log_unflipped(AmmoLa_intensity),
    Lateralization_square = square(Lateralization),
    CO2_threshold_slog = log_transform(CO2_threshold)
  )

# Check which one correlates with age
trigeminal_measures_data_age <- trigeminal_measures_data
trigeminal_measures_data_age$age <- as.numeric(trigeminale_daten_table1$Alter)
corr_mat_age <- cor(trigeminal_measures_data_age, use = "pairwise.complete.obs", method = "pearson")

# Check which one is sex different
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
  # Remove missing for fair comparison
  idx <- complete.cases(x, s)
  x <- x[idx]
  s <- s[idx]
  # Parametric ANOVA for eta squared (proportion variance explained by sex)
  aovmod <- aov(x ~ s)
  ss_total <- sum((x - mean(x)) ^ 2)
  ss_between <- sum(tapply(x, s, function(g) length(g) * (mean(g) - mean(x)) ^ 2))
  eta2 <- ss_between / ss_total
  # Non-parametric p-value as reference
  kruskal_p <- tryCatch(kruskal.test(x ~ s)$p.value, error = function(e) NA)
  n <- length(x)
  data.frame(
    variable = var,
    eta2 = eta2,
    kruskal_p = kruskal_p,
    n = n
  )
})

sex_stats_df <- bind_rows(sex_stats) %>%
  arrange(desc(eta2))

print(sex_stats_df)

# ----------------------------
# Plot distributions of non-CO2 related trigeminal measures
trigeminal_measures_data_long <- reshape2::melt(trigeminal_measures_data)

# Filter to exclude CO2 variables for this plot
trigeminal_measures_data_long_1 <- trigeminal_measures_data_long %>%
  filter(!str_detect(variable, "CO2"))

p_distribution_tigeminal_nonCO2 <- ggplot(trigeminal_measures_data_long_1, aes(x = value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = "lightblue", color = "white", alpha = 0.7
  ) +
  geom_density(color = "darkblue", size = 1) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

p_distribution_tigeminal_nonCO2

ggsave(paste0("p_distribution_tigeminal_nonCO2", ".svg"), p_distribution_tigeminal_nonCO2, width = 10, height = 10)

# ----------------------------
# Prepare CO2-related data and subsets for plotting
trigeminal_measures_data_CO2 <- trigeminal_measures_data[, str_detect(names(trigeminal_measures_data), "CO2")]

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
p_distribution_tigeminal_CO2 <- ggplot(trigeminal_measures_data_CO2_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightblue", color = "white", alpha = 0.7) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

p_distribution_tigeminal_CO2

ggsave(paste0("p_distribution_tigeminal_CO2", ".svg"), p_distribution_tigeminal_CO2, width = 10, height = 10)

# ----------------------------
# Analyse and plot correlations as matrix
# ------------ DATA PREP -----------------
trigeminal_measures_data_transformed <- trigeminal_measures_data[,
                                                                 c("AmmoLa_intensity_reflect_slog", "Lateralization_square", "CO2_threshold_slog")]

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

# --- Compose safe plotmath expressions for annotation ---
get_eq_label <- function(x, y) {
  mod <- MASS::rlm(y ~ x)
  coefs <- as.numeric(coef(mod))
  b0 <- formatC(coefs[1], digits = 2, format = "f")
  b1 <- formatC(coefs[2], digits = 2, format = "f")
  signb1 <- ifelse(coefs[2] < 0, "-", "+")
  yhat <- predict(mod)
  r2 <- 1 - sum((y - yhat) ^ 2) / sum((y - mean(y)) ^ 2)
  r2f <- formatC(r2, digits = 2, format = "f")
  # Proper plotmath: don't put *','* before a new == !
  eq <- sprintf(
    "italic(y)==%s %s %s*italic(x) *','~italic(R)^2==%s",
    b0, signb1, abs(as.numeric(b1)), r2f
  )
  eq
}

# Cor label needs same! Use only == ... for r, p
cor_stats <- map_dfr(var_pairs, function(vars) {
  x <- trigeminal_measures_data_transformed[[vars[1]]]
  y <- trigeminal_measures_data_transformed[[vars[2]]]
  cor_test <- cor.test(x, y, use = "complete.obs")
  data.frame(
    pair = paste(vars, collapse = " vs "),
    label = sprintf("italic(r)==%.2f * ','~italic(p)==%.3g",
                    cor_test$estimate, cor_test$p.value),
    stringsAsFactors = FALSE
  )
})

eq_stats <- map2_dfr(var_pairs, seq_along(var_pairs), function(vars, i) {
  x <- trigeminal_measures_data_transformed[[vars[1]]]
  y <- trigeminal_measures_data_transformed[[vars[2]]]
  idx <- complete.cases(x, y)
  eq <- get_eq_label(x[idx], y[idx])
  data.frame(pair = paste(vars, collapse = " vs "), eq_label = eq, stringsAsFactors = FALSE)
})

# --- Combine with ~ for math, no comma except in *','* ---
full_stats_df <- left_join(eq_stats, cor_stats, by = "pair") %>%
  mutate(top_label = paste(eq_label, label, sep = "~','~"))

# --- rest as before ---
facet_tops <- plot_df %>%
  group_by(pair) %>%
  summarise(
    x_label = min(x, na.rm = TRUE),
    y_label = max(y, na.rm = TRUE) + 0.07 * diff(range(y, na.rm = TRUE))
  )

label_df <- left_join(full_stats_df, facet_tops, by = "pair")

p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "rlm",
    formula = y ~ x,
    se = TRUE,
    fullrange = FALSE,
    color = "darkblue",
    fill = "lightblue",
    alpha = 0.22
  ) +
  geom_text(
    data = label_df,
    aes(x = x_label, y = y_label, label = top_label),
    parse = TRUE,
    color = "red",
    size = 4.5,
    hjust = 0,
    vjust = 0.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~pair, scales = "free", labeller = labeller(pair = function(x) {
    sapply(x, function(p) strsplit(p, " vs ")[[1]][2])
  })) +
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
    strip.background = element_blank()
  )

print(p)

# ----------------------------
# Prepare data for correlation heatmap --------------------------------------

# Mask CO2_threshold and transformed CO2_threshold for first 549 rows (set to NA)
# all_measures_correlations <- trigeminal_measures_data %>%
#   mutate(CO2_threshold_slog = ifelse(row_number() <= 549, NA, CO2_threshold_slog)) %>%
#   mutate(CO2_threshold = ifelse(row_number() <= 549, NA, CO2_threshold))

# Split groups, instead of masking
all_measures_correlations <- trigeminal_measures_data

# Fill new columns with NA, then assign only the correct rows
all_measures_correlations$CO2_threshold_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold[1:549]

all_measures_correlations$CO2_threshold_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold[550:nrow(all_measures_correlations)]

all_measures_correlations$CO2_threshold_slog_breath_not_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_not_hold[1:549] <- all_measures_correlations$CO2_threshold_slog[1:549]

all_measures_correlations$CO2_threshold_slog_breath_hold <- NA_real_
all_measures_correlations$CO2_threshold_slog_breath_hold[550:nrow(all_measures_correlations)] <- all_measures_correlations$CO2_threshold_slog[550:nrow(all_measures_correlations)]


# Calculate pairwise correlation matrix with pairwise complete observations
# corr_mat <- cor_stats(all_measures_correlations, use = "pairwise.complete.obs", method = "pearson")

library(psych)

# Select only numeric columns for correlation analysis
numeric_data <- all_measures_correlations %>%
  dplyr::select(where(is.numeric))

# Calculate correlations, p-values, and sample sizes, handling missing data with pairwise deletion
cor_method <- "spearman"
cor_results <- corr.test(numeric_data, use = "pairwise", method = cor_method)

# Extract correlation matrix
corr_mat <- cor_results$r

# Extract p-value matrix
p_mat <- cor_results$p

signif_code <- function(p) {
  if (is.na(p)) {
    ""
  } else if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    ""
  }
}


# Define NYT-inspired color palette for correlation heatmap
breaks <- c(
  -1, -0.9, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.02, -0.01, 0,
  0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 1
)

nyt_colors <- c(
  # reversed original color vector for negative values
  "ghostwhite",
  "#fbfbfb",
  "#e6f0fa",
  "#c9def9",
  "#add0fa",
  "#7bb8fa",
  "dodgerblue2",
  "#041a58"
)

nyt_colors <- c(
  rev(nyt_colors),
  nyt_colors
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
create_heatmap_trig <- function() {
  ht <- Heatmap(
    corr_mat,
    name = "Correlation",
    col = col_fun,
    na_col = "white",
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
      pval <- p_mat[i, j]
      star <- signif_code(pval)
      lbl <- paste0(val, star)
      col_text <- text_color_fun(fill)
      grid.text(lbl, x, y, gp = gpar(fontsize = 10, col = col_text))
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
    show_heatmap_legend = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),

  )

  # Render heatmap on a new page
  grid.newpage()
  draw(
    ht,
    heatmap_legend_side = "bottom"
  )

  # Add a bold, left-aligned title manually above heatmap
  grid.text(
    paste0("Correlation matrix (",cor_method,")"),
    x = unit(0, "npc") + unit(4, "mm"),  # left margin
    y = unit(1, "npc") - unit(4, "mm"),  # top margin
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
}

###############################################################################
# --- Export Plot as SVG ---
###############################################################################

# Capture the heatmap graphic output as a grid object for export
gp <- grid.grabExpr(create_heatmap_trig())
grid.draw(gp)

# Export the plot to an SVG file with specified dimensions
svg(paste0("trigeminal_correlation_heatmap", ".svg"), width = 12, height = 12)
grid.draw(gp)
dev.off()
