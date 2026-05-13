################################################################################
# UTILS.R - Utility Functions for Trigeminal Sensitivity Analysis
################################################################################
# Author: Jorn Lotsch
# Created: 2024
#
# This file contains all reusable utility functions organized by category.
# Load this file via: source("utils.R")
#
# TABLE OF CONTENTS:
#   § 1. Data Manipulation Helpers ............................ Line 30
#   § 2. Statistical Transformations .......................... Line 200
#   § 3. Imputation Functions ................................. Line 400
#   § 4. Date and Text Parsing ................................ Line 500
#   § 5. Statistical Effect Size Functions .................... Line 600
#   § 6. Plotting Helpers ..................................... Line 700
#
# USAGE NOTES:
#   - This file is automatically loaded by globals.R
#   - All functions include roxygen2-style documentation
#   - Functions are organized by category for easy navigation
#
################################################################################


# ============================================================================ #
# § 1. DATA MANIPULATION HELPERS
# ============================================================================ #
# Functions for data wrangling, one-hot encoding, type conversion, etc.
#
# Functions in this section:
#   - one_hot_encode_var()          : One-hot encode categorical variables
#   - convert_jn_to_numeric()       : Convert "j"/"n" to 1/0
#   - convert_truefalse_to_numeric(): Convert "TRUE"/"FALSE" to 1/0
#   - dim_wo_ID()                   : Get dimensions without ID column
#   - nVars_wo_ID()                 : Count variables without ID
#   - var_names_wo_ID()             : Get variable names without ID
#   - left_join_fill()              : Left join with NA filling
#
# ============================================================================ #

#' One-hot encode a nominal vector
#'
#' Converts a factor or character vector to a one-hot encoded data.frame.
#'
#' @param var The vector to encode (factor or character).
#' @param var_name Base name for new columns.
#' @return A data.frame with one column per factor level, values are 0/1.
one_hot_encode_var <- function(var, var_name) {
  factor_var <- factor(var, exclude = NULL)
  levels_var <- levels(factor_var)
  one_hot_list <- lapply(levels_var, function(lev) {
    if (is.na(lev)) {
      as.integer(is.na(var))
    } else {
      as.integer(ifelse(is.na(var), 0, var == lev))
    }
  })
  names(one_hot_list) <- paste0(var_name, "_", ifelse(is.na(levels_var), "NA", levels_var))
  return(as.data.frame(one_hot_list))
}

#' Convert yes/no (j/n) to numeric
#'
#' Recodes "j" to 1, "n" to 0, NA stays NA.
#'
#' @param x Character vector of "j", "n", or NA.
#' @return Numeric vector (1/0/NA).
convert_jn_to_numeric <- function(x) {
  ifelse(is.na(x), NA, ifelse(x == "j", 1, 0))
}

#' Convert TRUE/FALSE logicals to numeric
#'
#' Converts TRUE to 1, FALSE to 0, NA stays NA.
#'
#' @param x Logical vector (TRUE/FALSE/NA).
#' @return Numeric vector (1/0/NA).
convert_truefalse_to_numeric <- function(x) {
  ifelse(is.na(x), NA, ifelse(x == TRUE, 1, 0))
}

#' Get dimension of a data.frame without ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Vector: number of rows and columns (excluding "ID").
dim_wo_ID <- function(df) {
  if ("ID" %in% names(df)) {
    dim(df[, -1])
  } else dim(df)
}

#' Count variables in a data.frame without ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Integer, number of variables not counting "ID".
nVars_wo_ID <- function(df) {
  if ("ID" %in% names(df)) {
    length(names(df)) - 1
  } else length(names(df))
}

#' Get variable names excluding ID column
#'
#' @param df A data.frame or tibble with optional "ID" column.
#' @return Character vector of variable names (excluding "ID").
var_names_wo_ID <- function(df) {
  names(df)[!names(df) %in% "ID"]
}

#' Left join two dataframes and set NA values in new columns to a given value
#'
#' @param df_primary Left (main) dataframe
#' @param df_secondary Right (joined) dataframe
#' @param by Character vector of columns to join by (default: "ID")
#' @param fill_value Value to replace NAs with (default: 0)
#' @return A joined dataframe with NAs in added columns from df_secondary replaced by fill_value
left_join_fill <- function(df_primary, df_secondary, by = "ID", fill_value = 0) {
  # Do join
  df_combined <- dplyr::left_join(df_primary, df_secondary, by = by)

  # Find columns that were added from secondary dataframe
  added_cols <- setdiff(names(df_combined), names(df_primary))

  # For each added column, replace NA with fill_value
  df_combined <- df_combined %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(added_cols), ~ ifelse(is.na(.), fill_value, .)))

  return(df_combined)
}


# ============================================================================ #
# § 2. STATISTICAL TRANSFORMATIONS
# ============================================================================ #
# Functions for data transformations: log, scaling, normalization
#
# Functions in this section:
#   - slog()                        : Sign-preserving log transformation
#   - inv_slog()                    : Inverse of slog
#   - reflect_slog()                : Reflected log transformation
#   - reflect_slog_unflipped()      : Reflected log (unflipped)
#   - inv_reflect_slog_unflipped()  : Inverse reflected log (unflipped)
#   - scaleRange_01()               : Scale to [0, 1]
#   - scale01minmax()               : Min-max normalization with boundaries
#   - to_percent()                  : Scale to [0, 100]
#   - toRange()                     : Scale to arbitrary range
#
# ============================================================================ #

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

#' Scale vector to [0, 1] range
#'
#' @param x Numeric vector to scale
#' @return Scaled vector in range [0, 1]
scaleRange_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
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

#' Scale data to percentage range [0, 100]
#'
#' @param data Numeric vector or matrix to scale.
#' @return Scaled data where each column has minimum 0 and maximum 100.
#' @details Uses \code{toRange} internally.
to_percent <- function(data) {
  toRange(data, 0, 100)
}

#' Scale data to a specified range
#'
#' @param data Numeric vector or matrix to scale.
#' @param lower Numeric scalar specifying the lower bound of the target range.
#' @param upper Numeric scalar specifying the upper bound of the target range.
#' @return Scaled data where each column is transformed to lie within [lower, upper].
#' @details Each column is independently rescaled using min-max normalization.
#' Missing values are ignored when computing minima and maxima.
#' If a column has constant values (min == max), it is left unchanged.
#' If \code{lower > upper}, the values are swapped internally.
#' @examples
#' toRange(c(1, 2, 3), 0, 1)
#' toRange(matrix(1:6, ncol = 2), -1, 1)
toRange <- function(data, lower, upper){
  data <- as.matrix(data)
  if(lower==upper){
    error('interval width can not be 0!')
  }
  if (lower > upper){
    temp <- upper;
    upper <- lower;
    lower <- upper;
  }
  range <- upper - lower

  n <- dim(data)[1]
  d <- dim(data)[2]

  if ((n==1) & (d > 1)){ # row vector to colum vector
    data <- t(data)
    wasRowVector <- 1
  }
  else{
    wasRowVector <- 0
  }

  nRow <- dim(data)[1]
  nCol <- dim(data)[2]

  # Min = ones(Rows,1)*nanmin(data);
  min <-apply(data,2,min,na.rm=TRUE)
  min <- matrix(min,nRow,nCol,byrow=TRUE)
  # Max = ones(Rows,1)*nanmax(data);
  max <- apply(data,2,max,na.rm=TRUE)
  max <- matrix(max,nRow,nCol,byrow=TRUE)

  # Range = Max-Min;
  range <- max-min

  # Range(Range==0) =1; % falls Min==Max lass Daten in Ruhe
  range[range==0]<-1

  # ScaledData = (data-Min)./Range;          % scale to [0,1]

  scaleData <- (data-min)/range
  scaleData <- lower + scaleData * (upper-lower)

  if(wasRowVector==1){
    scaleData = t(scaleData)
  }

  return(scaleData)
}


# ============================================================================ #
# § 3. IMPUTATION FUNCTIONS
# ============================================================================ #
# Functions for handling missing data
#
# Functions in this section:
#   - detect_variable_type()        : Detect variable type for imputation
#   - fallback_impute()             : Simple fallback imputation
#   - impute_dataset()              : Main imputation function
#
# ============================================================================ #

# Functions will be added here in next step


# ============================================================================ #
# § 4. DATE AND TEXT PARSING
# ============================================================================ #
# Functions for parsing dates, years, and text descriptions
#
# Functions in this section:
#   - convert_mixed_dates()         : Standardize mixed date formats
#   - extract_years_custom()        : Extract years from text
#   - strip_years_time()            : Remove year/time from text
#
# ============================================================================ #

# Functions will be added here in next step


# ============================================================================ #
# § 5. STATISTICAL EFFECT SIZE FUNCTIONS
# ============================================================================ #
# Functions for calculating effect sizes and confidence intervals
#
# Functions in this section:
#   - signif_code()                 : Convert p-values to stars
#   - cohen_d_val()                 : Calculate Cohen's d
#   - boot_cohen_d()                : Bootstrap CI for Cohen's d
#   - bootstrap_cor()               : Bootstrap CI for correlation
#   - rank_biserial_val()           : Calculate rank-biserial correlation
#   - boot_rank_biserial()          : Bootstrap CI for rank-biserial
#
# ============================================================================ #

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

#' Compute Cohen's d effect size
#'
#' @param x First group values
#' @param y Second group values
#' @return Cohen's d estimate with Hedges correction
cohen_d_val <- function(x, y) {
  effsize::cohen.d(x, y, hedges.correction = TRUE)$estimate
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

  boot_out <- boot::boot(data, stat_fun, R = R)
  ci <- boot::boot.ci(boot_out, type = "perc")

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

  boot_obj <- boot::boot(data, boot_cor, R = R)
  ci <- boot::boot.ci(boot_obj, type = "perc")$percent[4:5]
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

  boot_out <- boot::boot(data, stat_fun, R = R)
  ci <- boot::boot.ci(boot_out, type = "perc")

  if (is.null(ci)) {
    return(c(NA, NA))
  } else {
    return(ci$percent[4:5])
  }
}


# ============================================================================ #
# § 6. PLOTTING HELPERS
# ============================================================================ #
# Functions for visualization and plot customization
#
# Functions in this section:
#   - PDEplotGG()                   : Pareto Density Estimation plot
#   - theme_plot()                  : Custom publication theme
#   - text_color_fun()              : Determine text color for contrast
#
# ============================================================================ #

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
    PDE <- DataVisualizations::ParetoDensityEstimation(as.vector(na.omit(Data[, i])))
    m2 <- PDE$kernels
    m3 <- PDE$paretoDensity
    m1 <- rep(i, length(m2))
    m <- rbind(m, cbind(m1, m2, m3))
  }

  mdf <- data.frame(m)
  require(ggplot2)

  p <- ggplot2::ggplot(data = mdf, ggplot2::aes(x = m2, y = m3, colour = factor(m1))) +
    ggplot2::geom_line(ggplot2::aes(linewidth = 1)) +
    ggplot2::guides(linewidth = "none")

  return(p)
}

#' Custom publication-quality ggplot2 theme
#'
#' @return ggplot2 theme object
theme_plot <- function() {
  ggplot2::theme_minimal(base_family = "Libre Franklin") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain", size = 12, color = "#222222",
        hjust = 0, margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(face = "plain", size = 10, color = "#444444"),
      axis.text = ggplot2::element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = ggplot2::element_text(
        size = 8, color = "#888888",
        hjust = 0, margin = ggplot2::margin(t = 10)
      ),
      panel.grid.major.y = ggplot2::element_line(
        color = "#dddddd", linetype = "dashed", size = 0.3
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks = ggplot2::element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks.length = grid::unit(5, "pt"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12, color = "#222222")
    )
}

#' Determine text color based on background brightness
#'
#' @param fill_color Background color
#' @return Text color (black or white) for optimal contrast
text_color_fun <- function(fill_color) {
  rgb_val <- grDevices::col2rgb(fill_color) / 255
  brightness <- 0.299 * rgb_val[1,] + 0.587 * rgb_val[2,] + 0.114 * rgb_val[3,]
  ifelse(brightness > 0.6, "#111111", "#FFFFFF")
}




# ============================================================================
#   § 7. PENALIZED REGRESSION MODELS
# ============================================================================
#   Combined framework for fitting and comparing penalized linear models
# (ridge, lasso, and elastic net), alongside standard linear regression
# Features:
#   - Automatic design matrix construction (handles factor encoding)
# - Cross-validated model tuning (lambda selection via cv.glmnet)
# - Parallel comparison of ridge (L2), lasso (L1), and elastic net
# - Integration of classical GLM estimates and p-values
# - Unified coefficient table with selection indicators
# Functions in this section:
#   - run_penalized_regression_all() : Fit models and summarize coefficients
# ============================================================================

run_penalized_regression_all <- function(train_data,
                                         train_target,
                                         alpha_elastic = 0.5,
                                         nfolds = 5,
                                         seed = 42,
                                         ridge_threshold = 0.05) {

  cat(sprintf("Penalized Linear Regression (ridge, lasso, elastic net) ===\n"))
  cat(sprintf("Dataset: %d features, %d samples\n",
              ncol(train_data), nrow(train_data)))

  if (ncol(train_data) == 0 || nrow(train_data) == 0) {
    cat("No data available - skipping\n")
    return(NULL)
  }

  ## ---------- 1. Prepare outcome and design matrix ----------
  y <- as.numeric(train_target)

  # model.matrix will create dummy variables if needed
  X <- model.matrix(~ ., data = train_data)[, -1, drop = FALSE]  # drop intercept column

  ## ---------- 2. Standard linear regression (for p values) ----------
  lr_train_data <- train_data
  lr_train_data$target <- y

  if (ncol(train_data) == 1) {
    formula_str <- paste("target ~", names(train_data)[1])
  } else {
    formula_str <- "target ~ ."
  }

  glm_fit <- glm(as.formula(formula_str), data = lr_train_data)

  # Get coefficient table, including p values, as a data frame
  glm_coef_df <- broom::tidy(glm_fit) %>%
    dplyr::filter(term != "(Intercept)") %>%        # drop intercept
    dplyr::select(variable = term,
                  glm_estimate = estimate,
                  glm_p = p.value)

  ## ---------- 3. Fit penalized models ----------
  set.seed(seed)

  # Helper to fit cv.glmnet and extract coefficients at lambda.min
  fit_penalized <- function(alpha_value) {
    cv_fit <- cv.glmnet(
      x = X,
      y = y,
      family = "gaussian",
      alpha = alpha_value,
      nfolds = nfolds
    )

    best_lambda <- cv_fit$lambda.min

    final_model <- glmnet(
      x = X,
      y = y,
      family = "gaussian",
      alpha = alpha_value,
      lambda = best_lambda
    )

    list(cv_fit = cv_fit,
         model = final_model,
         lambda = best_lambda)
  }

  ridge_res   <- fit_penalized(alpha_value = 0)
  lasso_res   <- fit_penalized(alpha_value = 1)
  elastic_res <- fit_penalized(alpha_value = alpha_elastic)

  cat(sprintf("ridge  lambda.min = %g\n", ridge_res$lambda))
  cat(sprintf("lasso  lambda.min = %g\n", lasso_res$lambda))
  cat(sprintf("elastic lambda.min = %g (alpha = %.2f)\n",
              elastic_res$lambda, alpha_elastic))

  ## ---------- 4. Extract coefficient vectors ----------
  get_coef_df <- function(res, label) {
    cf <- as.matrix(coef(res$model))   # includes intercept
    tibble(
      variable = rownames(cf),
      coef = as.numeric(cf)
    ) %>%
      filter(variable != "(Intercept)") %>%
      rename_with(~ paste0(label, "_coef"), .cols = coef)
  }

  ridge_coef_df   <- get_coef_df(ridge_res,   "ridge")
  lasso_coef_df   <- get_coef_df(lasso_res,   "lasso")
  elastic_coef_df <- get_coef_df(elastic_res, "elastic")

  ## ---------- 5. Merge all coefficients ----------
  coef_table <- glm_coef_df %>%
    full_join(ridge_coef_df,   by = "variable") %>%
    full_join(lasso_coef_df,   by = "variable") %>%
    full_join(elastic_coef_df, by = "variable")

  # Make sure we have all variables that appear in X, even if dropped in glm
  # (e.g., due to singularities)
  all_vars <- setdiff(colnames(X), "(Intercept)")
  coef_table <- coef_table %>%
    right_join(tibble(variable = all_vars), by = "variable") %>%
    arrange(variable)

  ## ---------- 6. Selection indicators ----------
  coef_table <- coef_table %>%
    mutate(
      ridge_selected   = if_else(!is.na(ridge_coef)   & abs(ridge_coef)   > ridge_threshold, TRUE, FALSE),
      lasso_selected   = if_else(!is.na(lasso_coef)   & lasso_coef   != 0, TRUE, FALSE),
      elastic_selected = if_else(!is.na(elastic_coef) & elastic_coef != 0, TRUE, FALSE)
    )

  ## ---------- 7. Print a compact table ----------
  cat("\n=== Variable selection summary ===\n")
  print(
    coef_table %>%
      dplyr::select(variable,
                    glm_p,
                    ridge_coef, ridge_selected,
                    lasso_coef, lasso_selected,
                    elastic_coef, elastic_selected)
  )

  ## ---------- 8. Return everything for further use ----------
  invisible(list(
    glm_fit = glm_fit,
    ridge   = ridge_res,
    lasso   = lasso_res,
    elastic = elastic_res,
    coef_table = coef_table
  ))
}


run_penalized_regression_remove_colinear_all <- function(train_data,
                                         train_target,
                                         alpha_elastic = 0.5,
                                         nfolds = 5,
                                         seed = 42,
                                         ridge_threshold = 0.05,
                                         vif_limit = 10) {

  cat(sprintf("Penalized Linear Regression (ridge, lasso, elastic net) ===\n"))
  cat(sprintf("Dataset: %d features, %d samples\n",
              ncol(train_data), nrow(train_data)))

  if (ncol(train_data) == 0 || nrow(train_data) == 0) {
    cat("No data available - skipping\n")
    return(NULL)
  }

  ## ---------- 1. Prepare outcome and design matrix ----------
  y <- as.numeric(train_target)

  # model.matrix will create dummy variables if needed
  X <- model.matrix(~ ., data = train_data)[, -1, drop = FALSE]  # drop intercept column

  ## ---------- 2. Standard linear regression (for p values) with remval of colinear or aliased vars ----------


  actual_data_base <- cbind.data.frame(
    train_data,
    target = y
  )

  current_data <- actual_data_base
  i <- 0

  repeat {
    i <- i + 1
    cat("Iteration:", i, "\n")

    lm_res <- lm(target ~ ., data = current_data)

    aliased_mat <- alias(lm_res)$Complete
    aliased_vars <- character(0)

    if (!is.null(aliased_mat)) {
      aliased_vars <- rownames(aliased_mat)
      aliased_vars <- gsub("`", "", aliased_vars)
      aliased_vars <- intersect(aliased_vars, names(current_data))
      cat("Aliased:", paste(aliased_vars, collapse = ", "), "\n")
    }

    if (length(aliased_vars) > 0) {
      before <- ncol(current_data)
      current_data <- current_data[, !names(current_data) %in% aliased_vars, drop = FALSE]
      after <- ncol(current_data)

      if (before == after) stop("Aliased terms found, but no matching columns were removed.")
      next
    }

    vif_vals <- car::vif(lm_res)
    if (is.matrix(vif_vals)) vif_vals <- vif_vals[, 1]

    max_vif <- max(vif_vals, na.rm = TRUE)
    cat("Max VIF:", max_vif, "\n")

    if (is.finite(max_vif) && max_vif <= vif_limit) break

    var_to_remove <- names(which.max(vif_vals))
    if (length(var_to_remove) == 0 || is.na(var_to_remove)) stop("Could not identify a variable to remove.")

    current_data <- current_data[, !names(current_data) %in% var_to_remove, drop = FALSE]
  }

  lr_train_data <- current_data

  if (ncol(current_data) == 2) {
    formula_str <- paste("target ~", names(current_data)[2])
  } else {
    formula_str <- "target ~ ."
  }

  glm_fit <- glm(as.formula(formula_str), data = lr_train_data)

  # Get coefficient table, including p values, as a data frame
  glm_coef_df <- broom::tidy(glm_fit) %>%
    dplyr::filter(term != "(Intercept)") %>%        # drop intercept
    dplyr::select(variable = term,
                  glm_estimate = estimate,
                  glm_p = p.value)

  ## ---------- 3. Fit penalized models ----------
  set.seed(seed)

  # Helper to fit cv.glmnet and extract coefficients at lambda.min
  fit_penalized <- function(alpha_value) {
    cv_fit <- cv.glmnet(
      x = X,
      y = y,
      family = "gaussian",
      alpha = alpha_value,
      nfolds = nfolds
    )

    best_lambda <- cv_fit$lambda.min

    final_model <- glmnet(
      x = X,
      y = y,
      family = "gaussian",
      alpha = alpha_value,
      lambda = best_lambda
    )

    list(cv_fit = cv_fit,
         model = final_model,
         lambda = best_lambda)
  }

  ridge_res   <- fit_penalized(alpha_value = 0)
  lasso_res   <- fit_penalized(alpha_value = 1)
  elastic_res <- fit_penalized(alpha_value = alpha_elastic)

  cat(sprintf("ridge  lambda.min = %g\n", ridge_res$lambda))
  cat(sprintf("lasso  lambda.min = %g\n", lasso_res$lambda))
  cat(sprintf("elastic lambda.min = %g (alpha = %.2f)\n",
              elastic_res$lambda, alpha_elastic))

  ## ---------- 4. Extract coefficient vectors ----------
  get_coef_df <- function(res, label) {
    cf <- as.matrix(coef(res$model))   # includes intercept
    tibble(
      variable = rownames(cf),
      coef = as.numeric(cf)
    ) %>%
      filter(variable != "(Intercept)") %>%
      rename_with(~ paste0(label, "_coef"), .cols = coef)
  }

  ridge_coef_df   <- get_coef_df(ridge_res,   "ridge")
  lasso_coef_df   <- get_coef_df(lasso_res,   "lasso")
  elastic_coef_df <- get_coef_df(elastic_res, "elastic")

  ## ---------- 5. Merge all coefficients ----------
  coef_table <- glm_coef_df %>%
    full_join(ridge_coef_df,   by = "variable") %>%
    full_join(lasso_coef_df,   by = "variable") %>%
    full_join(elastic_coef_df, by = "variable")

  # Make sure we have all variables that appear in X, even if dropped in glm
  # (e.g., due to singularities)
  all_vars <- setdiff(colnames(X), "(Intercept)")
  coef_table <- coef_table %>%
    right_join(tibble(variable = all_vars), by = "variable") %>%
    arrange(variable)

  ## ---------- 6. Selection indicators ----------
  coef_table <- coef_table %>%
    mutate(
      ridge_selected   = if_else(!is.na(ridge_coef)   & abs(ridge_coef)   > ridge_threshold, TRUE, FALSE),
      lasso_selected   = if_else(!is.na(lasso_coef)   & lasso_coef   != 0, TRUE, FALSE),
      elastic_selected = if_else(!is.na(elastic_coef) & elastic_coef != 0, TRUE, FALSE)
    )

  ## ---------- 7. Print a compact table ----------
  cat("\n=== Variable selection summary ===\n")
  print(
    coef_table %>%
      dplyr::select(variable,
                    glm_p,
                    ridge_coef, ridge_selected,
                    lasso_coef, lasso_selected,
                    elastic_coef, elastic_selected)
  )

  ## ---------- 8. Return everything for further use ----------
  invisible(list(
    glm_fit = glm_fit,
    ridge   = ridge_res,
    lasso   = lasso_res,
    elastic = elastic_res,
    coef_table = coef_table
  ))
}


run_penalized_multinomial_all <- function(train_data,
                                          train_target,
                                          alpha_elastic = 0.5,
                                          nfolds = 5,
                                          seed = 42,
                                          ridge_threshold = 0.05,
                                          any_class_selection = TRUE) {
  cat(sprintf("Penalized Multinomial Regression (ridge, lasso, elastic net) ===\n"))
  cat(sprintf(
    "Dataset: %d features, %d samples\n",
    ncol(train_data), nrow(train_data)
  ))

  if (ncol(train_data) == 0 || nrow(train_data) == 0) {
    cat("No data available - skipping\n")
    return(NULL)
  }

  ## 1. Outcome and design matrix
  y <- as.factor(train_target)
  X <- model.matrix(~., data = train_data)[, -1, drop = FALSE]

  ## ---------- NEW: Unpenalized multinomial model ----------
  df_glm <- train_data
  df_glm$target <- y

  multinom_fit <- nnet::multinom(target ~ ., data = df_glm, trace = FALSE)

  # Extract coefficients (matrix: classes x variables)
  coef_mat <- coef(multinom_fit)

  # Handle binary case (returns vector instead of matrix)
  if (is.vector(coef_mat)) {
    coef_mat <- matrix(coef_mat, nrow = 1)
    rownames(coef_mat) <- levels(y)[2]
  }

  glm_long <- as.data.frame(coef_mat) %>%
    tibble::rownames_to_column("class") %>%
    tidyr::pivot_longer(-class, names_to = "variable", values_to = "glm_coef")

  s <- summary(multinom_fit)
  z_mat <- s$coefficients / s$standard.errors
  p_mat <- (1 - pnorm(abs(z_mat), 0, 1)) * 2
  if (is.vector(p_mat)) {
    p_mat <- matrix(p_mat, nrow = 1)
    rownames(p_mat) <- levels(y)[2]
  }
  glm_pval_long <- as.data.frame(p_mat) %>%
    tibble::rownames_to_column("class") %>%
    tidyr::pivot_longer(-class, names_to = "variable", values_to = "glm_pval")

  ## --------------------------------------------------------

  set.seed(seed)

  fit_penalized_multinom <- function(alpha_value) {
    cv_fit <- glmnet::cv.glmnet(
      x = X,
      y = y,
      family = "multinomial",
      alpha = alpha_value,
      nfolds = nfolds
    )
    best_lambda <- cv_fit$lambda.min
    final_model <- glmnet::glmnet(
      x = X,
      y = y,
      family = "multinomial",
      alpha = alpha_value,
      lambda = best_lambda
    )
    list(
      cv_fit = cv_fit,
      model = final_model,
      lambda = best_lambda
    )
  }

  ridge_res <- fit_penalized_multinom(alpha_value = 0)
  lasso_res <- fit_penalized_multinom(alpha_value = 1)
  elastic_res <- fit_penalized_multinom(alpha_value = alpha_elastic)

  cat(sprintf("ridge   lambda.min = %g\n", ridge_res$lambda))
  cat(sprintf("lasso   lambda.min = %g\n", lasso_res$lambda))
  cat(sprintf(
    "elastic lambda.min = %g (alpha = %.2f)\n",
    elastic_res$lambda, alpha_elastic
  ))

  ## Extract penalized coefficients
  get_coef_long <- function(res, label) {
    cf_list <- coef(res$model)
    purrr::map_dfr(names(cf_list), function(cls) {
      mat <- as.matrix(cf_list[[cls]])
      tibble::tibble(
        class    = cls,
        variable = rownames(mat),
        coef     = as.numeric(mat)
      )
    }) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      dplyr::rename(!!paste0(label, "_coef") := coef)
  }

  ridge_long <- get_coef_long(ridge_res, "ridge")
  lasso_long <- get_coef_long(lasso_res, "lasso")
  elastic_long <- get_coef_long(elastic_res, "elastic")

  ## ---------- UPDATED MERGE: include glm ----------
  coef_table <- glm_long %>%
    dplyr::full_join(glm_pval_long, by = c("variable", "class")) %>%
    dplyr::full_join(ridge_long, by = c("variable", "class")) %>%
    dplyr::full_join(lasso_long, by = c("variable", "class")) %>%
    dplyr::full_join(elastic_long, by = c("variable", "class")) %>%
    dplyr::arrange(variable, class)

  ## Selection indicators
  coef_table <- coef_table %>%
    dplyr::mutate(
      glm_selected     = dplyr::if_else(!is.na(glm_pval) & glm_pval < 0.05, TRUE, FALSE),
      ridge_selected   = dplyr::if_else(!is.na(ridge_coef) & abs(ridge_coef) > ridge_threshold, TRUE, FALSE),
      lasso_selected   = dplyr::if_else(!is.na(lasso_coef) & lasso_coef != 0, TRUE, FALSE),
      elastic_selected = dplyr::if_else(!is.na(elastic_coef) & elastic_coef != 0, TRUE, FALSE)
    )

  cat("\n=== Variable selection summary (per class) ===\n")
  print(
    coef_table %>%
      dplyr::select(
        variable, class,
        glm_coef, glm_selected,
        ridge_coef, ridge_selected,
        lasso_coef, lasso_selected,
        elastic_coef, elastic_selected
      )
  )

  if (any_class_selection) {
    var_level <- coef_table %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        ridge_selected_any = any(ridge_selected),
        lasso_selected_any = any(lasso_selected),
        elastic_selected_any = any(elastic_selected),
        .groups = "drop"
      )
    cat("\n=== Variable-level selection (any class) ===\n")
    print(var_level)
  }

  invisible(list(
    multinom_fit = multinom_fit, # NEW
    ridge = ridge_res,
    lasso = lasso_res,
    elastic = elastic_res,
    coef_table = coef_table
  ))
}


run_penalized_multinomial_remove_colinear_all <- function(train_data,
                                          train_target,
                                          alpha_elastic = 0.5,
                                          nfolds = 5,
                                          seed = 42,
                                          ridge_threshold = 0.05,
                                          any_class_selection = TRUE,
                                          vif_limit = 10) {
  cat(sprintf("Penalized Multinomial Regression (ridge, lasso, elastic net) ===\n"))
  cat(sprintf(
    "Dataset: %d features, %d samples\n",
    ncol(train_data), nrow(train_data)
  ))

  if (ncol(train_data) == 0 || nrow(train_data) == 0) {
    cat("No data available - skipping\n")
    return(NULL)
  }

  ## 1. Outcome and design matrix
  y <- as.factor(train_target)
  X <- model.matrix(~., data = train_data)[, -1, drop = FALSE]  # full matrix for penalized models

  ## ---------- 2. Remove aliased and collinear predictors (lm proxy) ----------
  current_data <- train_data
  i <- 0

  repeat {
    i <- i + 1
    cat("Iteration:", i, "\n")

    lm_proxy <- lm(as.numeric(y) ~ ., data = cbind.data.frame(target = as.numeric(y), current_data))

    aliased_mat <- alias(lm_proxy)$Complete
    aliased_vars <- character(0)

    if (!is.null(aliased_mat)) {
      aliased_vars <- rownames(aliased_mat)
      aliased_vars <- gsub("`", "", aliased_vars)
      aliased_vars <- intersect(aliased_vars, names(current_data))
      cat("Aliased:", paste(aliased_vars, collapse = ", "), "\n")
    }

    if (length(aliased_vars) > 0) {
      before <- ncol(current_data)
      current_data <- current_data[, !names(current_data) %in% aliased_vars, drop = FALSE]
      after <- ncol(current_data)
      if (before == after) stop("Aliased terms found, but no matching columns were removed.")
      next
    }

    vif_vals <- car::vif(lm_proxy)
    if (is.matrix(vif_vals)) vif_vals <- vif_vals[, 1]

    max_vif <- max(vif_vals, na.rm = TRUE)
    cat("Max VIF:", max_vif, "\n")

    if (is.finite(max_vif) && max_vif <= vif_limit) break

    var_to_remove <- names(which.max(vif_vals))
    if (length(var_to_remove) == 0 || is.na(var_to_remove)) stop("Could not identify a variable to remove.")

    current_data <- current_data[, !names(current_data) %in% var_to_remove, drop = FALSE]
  }

  ## ---------- 3. Unpenalized multinomial on collinearity-cleaned data ----------
  df_glm <- current_data
  df_glm$target <- y

  multinom_fit <- nnet::multinom(target ~ ., data = df_glm, trace = FALSE)

  # Extract coefficients (matrix: classes x variables)
  coef_mat <- coef(multinom_fit)

  # Handle binary case (returns vector instead of matrix)
  if (is.vector(coef_mat)) {
    coef_mat <- matrix(coef_mat, nrow = 1)
    rownames(coef_mat) <- levels(y)[2]
  }

  glm_long <- as.data.frame(coef_mat) %>%
    tibble::rownames_to_column("class") %>%
    tidyr::pivot_longer(-class, names_to = "variable", values_to = "glm_coef")

  s <- summary(multinom_fit)
  z_mat <- s$coefficients / s$standard.errors
  p_mat <- (1 - pnorm(abs(z_mat), 0, 1)) * 2
  if (is.vector(p_mat)) {
    p_mat <- matrix(p_mat, nrow = 1)
    rownames(p_mat) <- levels(y)[2]
  }
  glm_pval_long <- as.data.frame(p_mat) %>%
    tibble::rownames_to_column("class") %>%
    tidyr::pivot_longer(-class, names_to = "variable", values_to = "glm_pval")

  ## --------------------------------------------------------

  set.seed(seed)

  fit_penalized_multinom <- function(alpha_value) {
    cv_fit <- glmnet::cv.glmnet(
      x = X,
      y = y,
      family = "multinomial",
      alpha = alpha_value,
      nfolds = nfolds
    )
    best_lambda <- cv_fit$lambda.min
    final_model <- glmnet::glmnet(
      x = X,
      y = y,
      family = "multinomial",
      alpha = alpha_value,
      lambda = best_lambda
    )
    list(
      cv_fit = cv_fit,
      model = final_model,
      lambda = best_lambda
    )
  }

  ridge_res <- fit_penalized_multinom(alpha_value = 0)
  lasso_res <- fit_penalized_multinom(alpha_value = 1)
  elastic_res <- fit_penalized_multinom(alpha_value = alpha_elastic)

  cat(sprintf("ridge   lambda.min = %g\n", ridge_res$lambda))
  cat(sprintf("lasso   lambda.min = %g\n", lasso_res$lambda))
  cat(sprintf(
    "elastic lambda.min = %g (alpha = %.2f)\n",
    elastic_res$lambda, alpha_elastic
  ))

  ## Extract penalized coefficients
  get_coef_long <- function(res, label) {
    cf_list <- coef(res$model)
    purrr::map_dfr(names(cf_list), function(cls) {
      mat <- as.matrix(cf_list[[cls]])
      tibble::tibble(
        class    = cls,
        variable = rownames(mat),
        coef     = as.numeric(mat)
      )
    }) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      dplyr::rename(!!paste0(label, "_coef") := coef)
  }

  ridge_long <- get_coef_long(ridge_res, "ridge")
  lasso_long <- get_coef_long(lasso_res, "lasso")
  elastic_long <- get_coef_long(elastic_res, "elastic")

  ## ---------- UPDATED MERGE: include glm ----------
  coef_table <- glm_long %>%
    dplyr::full_join(glm_pval_long, by = c("variable", "class")) %>%
    dplyr::full_join(ridge_long, by = c("variable", "class")) %>%
    dplyr::full_join(lasso_long, by = c("variable", "class")) %>%
    dplyr::full_join(elastic_long, by = c("variable", "class")) %>%
    dplyr::arrange(variable, class)

  ## Selection indicators
  coef_table <- coef_table %>%
    dplyr::mutate(
      glm_selected     = dplyr::if_else(!is.na(glm_pval) & glm_pval < 0.05, TRUE, FALSE),
      ridge_selected   = dplyr::if_else(!is.na(ridge_coef) & abs(ridge_coef) > ridge_threshold, TRUE, FALSE),
      lasso_selected   = dplyr::if_else(!is.na(lasso_coef) & lasso_coef != 0, TRUE, FALSE),
      elastic_selected = dplyr::if_else(!is.na(elastic_coef) & elastic_coef != 0, TRUE, FALSE)
    )

  cat("\n=== Variable selection summary (per class) ===\n")
  print(
    coef_table %>%
      dplyr::select(
        variable, class,
        glm_coef, glm_selected,
        ridge_coef, ridge_selected,
        lasso_coef, lasso_selected,
        elastic_coef, elastic_selected
      )
  )

  if (any_class_selection) {
    var_level <- coef_table %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        ridge_selected_any = any(ridge_selected),
        lasso_selected_any = any(lasso_selected),
        elastic_selected_any = any(elastic_selected),
        .groups = "drop"
      )
    cat("\n=== Variable-level selection (any class) ===\n")
    print(var_level)
  }

  invisible(list(
    multinom_fit = multinom_fit,
    ridge = ridge_res,
    lasso = lasso_res,
    elastic = elastic_res,
    coef_table = coef_table
  ))
}

################################################################################
# END OF UTILS.R
################################################################################
# All utility functions have been successfully migrated and organized.
# This file is sourced automatically by globals.R.
################################################################################
