################################################################################
# Trigeminal Sensitivity Analysis - Modulator Regression Analyses
#
# Author: Jorn Lotsch
# Date: 2026-04-29
#
# Description:
#   Regression analyses linking candidate modulator variables (demographics,
#   clinical, and other non-trigeminal variables) to three trigeminal
#   psychophysical outcomes (Lateralization, AmmoLa intensity, CO2 threshold).
#   Combines Boruta feature selection with penalized regression (ridge, lasso,
#   elastic net) and Random Forest regression. Models are trained on an imputed
#   training dataset and validated on an independent dataset. Confirmed
#   modulators are further characterised with Wilcoxon or Spearman tests in
#   both training and validation data.
#
# Input Files:
#   - analysis_dataset_training_imputed.csv: Training data with imputed values
#   - analysis_dataset_validation_imputed.csv: Validation data with imputed values
#   - globals.R: Global variables, color palettes, and transformation functions
#   - utils.R: Utility functions
#
# Output Files:
#   - regression_modulators_coef_tables_combined.csv: Coefficients per outcome
#   - regression_modulators_final_table_modulators.csv: Validation performance
#   - regression_modulators_sig_features_stats.csv: Stats for confirmed modulators
#   - combined_plots_reg_modulators_psy.svg/.png: Predicted vs. observed plots
#   - p_ammoLa_Boruta.svg/.png: Boruta importance plot for AmmoLa intensity
#
# Key Analyses:
#   1. Data loading, variable selection, and transformation
#   2. Boruta feature selection (Random Forest-based importance)
#   3. Penalized linear regression (ridge, lasso, elastic net via glmnet)
#   4. Standard GLM for coefficient p-values
#   5. RF regression with tuned hyperparameters (mtry, ntree)
#   6. Model validation: predicted vs. observed on independent dataset
#   7. Post-hoc characterisation of confirmed modulators (Wilcoxon/Spearman)
#   8. Exploratory sex-difference analysis
#
# Statistical Methods:
#   - Boruta algorithm (RF-based feature selection)
#   - Penalized regression: ridge (alpha=0), lasso (alpha=1), elastic net (alpha=0.5)
#   - Cross-validation for lambda tuning (cv.glmnet, 5-fold)
#   - NLS for prediction assessment; ANOVA for model comparison
#   - Wilcoxon rank-sum test (binary modulators); Spearman correlation (continuous)
################################################################################

# ============================================================================ #
# 1. LOAD REQUIRED LIBRARIES
# ============================================================================ #

# Load custom functions and global variables
source("globals.R")
source("utils.R")

# ============================================================================ #
# 2. GLOBAL OPTIONS AND ANALYSIS PARAMETERS
# ============================================================================ #

# Computational resources
nProc_possible <- parallel::detectCores() - 1
nProc_desired <- 24

# Color palettes
original_colors <- c(actual_palette[1], actual_palette[2], actual_palette[3], actual_palette[4])
dark_colors <- darken(original_colors, amount = 0.3)

# ============================================================================ #
# 3. HELPER FUNCTIONS
# ============================================================================ #

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


prepare_boruta_plot_data <- function(boruta_res, decision_col = "Boruta.res") {
  # Melt importance history
  imp_long <- reshape2::melt(boruta_res$ImpHistory)
  colnames(imp_long) <- c("Iteration", "Feature", "Importance")

  # Extract final decisions
  decisions <- data.frame(Decision = boruta_res$finalDecision)
  decisions$Feature <- rownames(decisions)

  # Assign color categories
  imp_long$Color <- ifelse(
    imp_long$Feature %in% decisions$Feature[decisions$Decision == "Confirmed"], "Chosen",
    ifelse(
      imp_long$Feature %in% decisions$Feature[decisions$Decision %in% c("Tentative")],
      "Tentative",
      ifelse(
        imp_long$Feature %in% decisions$Feature[decisions$Decision %in% c("Rejected")],
        "Rejected", "Shadow")
    )
  )

  # Ensure consistent factor levels for colors
  imp_long$Color <- factor(imp_long$Color, levels = c("Chosen", "Tentative", "Rejected", "Shadow"))

  list(importance = imp_long, decisions = decisions)
}


extract_anova <- function(anova_table) {
  anova_df <- as.data.frame(anova_table)
  f_value <- anova_df$`F value`[2]
  p_value <- anova_df$`Pr(>F)`[2]
  return(list(f_value = f_value, p_value = p_value))
}


extract_full_model_coefs <- function(full_model) {
  model_df <- as.data.frame(full_model$coefficients)
  list(
    const_est = model_df["const", "Estimate"],
    const_p   = model_df["const", "Pr(>|t|)"],
    slope_est = model_df["slope", "Estimate"],
    slope_p   = model_df["slope", "Pr(>|t|)"]
  )
}

get_selected_vars_mod <- function(fit_type, coef_table, ref_data) {
  if (fit_type == "glm_fit") {
    mm_vars <- unique(coef_table$variable[!is.na(coef_table$glm_p) & coef_table$glm_p < 0.05])
  } else {
    col     <- switch(fit_type,
      "elastic" = "elastic_selected",
      "lasso"   = "lasso_selected",
      "ridge"   = "ridge_selected"
    )
    mm_vars <- unique(coef_table$variable[coef_table[[col]] == TRUE])
  }
  mm_vars <- trimws(gsub("`", "", mm_vars))
  # Map model.matrix column names back to original data column names
  names(ref_data)[
    names(ref_data) %in% mm_vars |
    make.names(names(ref_data)) %in% make.names(mm_vars)
  ]
}

retrain_glm <- function(train_data, train_target) {
  lr_data        <- train_data
  lr_data$target <- as.numeric(train_target)
  formula_str <- if (ncol(train_data) == 1) paste("target ~", names(train_data)[1]) else "target ~ ."
  glm(as.formula(formula_str), data = lr_data)
}

retrain_penalized_linear <- function(fit_type, train_data, train_target, nfolds = 5, seed = 42) {
  alpha_val   <- switch(fit_type, "elastic" = 0.5, "lasso" = 1, "ridge" = 0)
  y           <- as.numeric(train_target)
  X           <- model.matrix(~ ., data = train_data)[, -1, drop = FALSE]
  set.seed(seed)
  cv_fit      <- glmnet::cv.glmnet(x = X, y = y, family = "gaussian", alpha = alpha_val, nfolds = nfolds)
  final_model <- glmnet::glmnet(x = X, y = y, family = "gaussian", alpha = alpha_val, lambda = cv_fit$lambda.min)
  list(model = final_model, lambda = cv_fit$lambda.min)
}

# ============================================================================ #
# 4. LOAD AND PREPARE DATA
# ============================================================================ #

cat("\n=== Loading Data ===\n")

# Load training data (imputed)
trigeminale_training_data <- read.csv("analysis_dataset_training_imputed.csv",
                                      check.names = FALSE)
cat("Training data loaded:", nrow(trigeminale_training_data), "samples\n")

# Extract variables for regression
variables_for_fs_imputed <- trigeminale_training_data[,!names(trigeminale_training_data) %in% c(variables_by_categories$Nasal_chemosensory_perception,
                                                                                                variables_by_categories$Psychophysical_measurements[c(1, 2, 4)], "ID")]
targets_for_fs_imputed <- trigeminale_training_data[,variables_by_categories$Psychophysical_measurements[c(1, 2, 4)]]

rownames(variables_for_fs_imputed) <- trigeminale_training_data$ID
rownames(targets_for_fs_imputed) <- trigeminale_training_data$ID


# Load validation data (imputed)
trigeminale_validation_data <- read.csv("analysis_dataset_validation_imputed.csv",
                                        check.names = FALSE)
cat("Validation data loaded:", nrow(trigeminale_validation_data), "samples\n")

# Extract variables for regression (validation)
variables_for_fs_imputed_validation <- trigeminale_validation_data[,!names(trigeminale_validation_data) %in% c(variables_by_categories$Nasal_chemosensory_perception,
                                                                                                variables_by_categories$Psychophysical_measurements[c(1, 2, 4)], "ID")]
targets_for_fs_imputed_validation <- trigeminale_validation_data[,variables_by_categories$Psychophysical_measurements[c(1, 2, 4)]]

rownames(variables_for_fs_imputed_validation) <- trigeminale_validation_data$ID
rownames(targets_for_fs_imputed_validation) <- trigeminale_validation_data$ID

# ============================================================================ #
# 5. DATA TRANSFORMATION AND PREPROCESSING
# ============================================================================ #

cat("\n=== Transforming Variables ===\n")

# Transform psychophysical targets (training)
# AmmoLa intensity: reflected slog (handles ceiling effects)
targets_for_fs_imputed$`AmmoLa intensity` <-
  reflect_slog_unflipped(targets_for_fs_imputed$`AmmoLa intensity`)

# CO2 threshold: negative slog (lower = more sensitive)
targets_for_fs_imputed$`CO2 threshold` <-
  -slog(targets_for_fs_imputed$`CO2 threshold`)

# Lateralization: none
targets_for_fs_imputed$`Lateralization (x/20)` <-
  targets_for_fs_imputed$`Lateralization (x/20)`

# Transform psychophysical targets (validation)
targets_for_fs_imputed_validation$`AmmoLa intensity` <-
  reflect_slog_unflipped(targets_for_fs_imputed_validation$`AmmoLa intensity`)

targets_for_fs_imputed_validation$`CO2 threshold` <-
  -slog(targets_for_fs_imputed_validation$`CO2 threshold`)

targets_for_fs_imputed_validation$`Lateralization (x/20)` <-
  targets_for_fs_imputed_validation$`Lateralization (x/20)`

cat("Targets transformed to [0, 3] scale\n")


# Inspect dataset structure
cat("\nDataset structure:\n")
print(names(targets_for_fs_imputed))

# ============================================================================ #
# 6. BORUTA FEATURE SELECTION AND PENALIZED REGRESSION
# ============================================================================ #

reg_modulators_resultat <- pbmcapply::pbmclapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  actual_data <-  variables_for_fs_imputed
  set.seed(42)

  # Quick tune RF for REGRESSION
  mtry_values <- unique(round(c(2, ncol(actual_data) / 3, ncol(actual_data) / 2, sqrt(ncol(actual_data)), log2(ncol(actual_data)), ncol(actual_data))))  # regression default is ncol/3, not sqrt
  ntree_values <- c(100, 200, 500, 1000)
  grid <- expand.grid(mtry = mtry_values, ntree = ntree_values)
  grid$median_mse <- NA
  mses_matrix <- matrix(NA, nrow = nrow(grid), ncol = 10)

  for (run in 1:10) {
    set.seed(42 + run)
    for (i in 1:nrow(grid)) {
      model <- randomForest(
        x = actual_data,
        y = targets_for_fs_imputed[[reg_target_varname]],
        mtry = grid$mtry[i],
        ntree = grid$ntree[i]
      )
      mses_matrix[i, run] <- model$mse[length(model$mse)]  # final OOB MSE
    }
  }

  grid$median_mse <- apply(mses_matrix, 1, median)
  best <- grid[which.min(grid$median_mse), ]

  Boruta_res <- Boruta::Boruta(x = actual_data,
                               y = targets_for_fs_imputed[[reg_target_varname]], ntree = best[["ntree"]], mtry = best[["mtry"]], maxRuns = 1000, doTrace = TRUE)

  penalized_regression_res <- run_penalized_regression_all(
    train_data   = actual_data,
    train_target = targets_for_fs_imputed[[reg_target_varname]],
    alpha_elastic = 0.5,
    nfolds = 5,
    ridge_threshold = 0.05  # or something small that makes sense in your scale
  )
  return(list(Boruta_res = Boruta_res, rf_best_params = best, penalized_regression_res = penalized_regression_res))
}, mc.cores = nProc_desired)

names(reg_modulators_resultat) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

# ============================================================================ #
# 7. INITIAL RESULT INSPECTION
# ============================================================================ #

par(mfrow=c(2,3))
plot(reg_modulators_resultat[[1]]$Boruta_res, las = 2)
plot(reg_modulators_resultat[[2]]$Boruta_res, las = 2)
plot(reg_modulators_resultat[[3]]$Boruta_res, las = 2)
par(mfrow=c(1,1))

print(reg_modulators_resultat[["Lateralization (x/20)"]]$penalized_regression_res)
print(reg_modulators_resultat[["AmmoLa intensity"]]$penalized_regression_res)
print(reg_modulators_resultat[["CO2 threshold"]]$penalized_regression_res)

# ============================================================================ #
# 8. REGRESSION MODEL PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

reg_modulators_predictions <- lapply(c("all", "reduced"), function(feature_set) {
  result <- lapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
    print(paste(feature_set, reg_target_varname))
    ct <- reg_modulators_resultat[[reg_target_varname]]$penalized_regression_res$coef_table

    fits <- lapply(c("glm_fit", "elastic", "lasso", "ridge"), function(fit_type) {
      print(fit_type)

      if (feature_set == "all") {
        val_data <- variables_for_fs_imputed_validation[, !names(variables_for_fs_imputed_validation) %in% reg_target_varname]
        if (fit_type == "glm_fit") {
          y_new <- stats::predict(
            object  = reg_modulators_resultat[[reg_target_varname]]$penalized_regression_res$glm_fit,
            newdata = val_data
          )
        } else {
          y_new <- stats::predict(
            object = reg_modulators_resultat[[reg_target_varname]]$penalized_regression_res[[fit_type]]$model,
            newx   = as.matrix(val_data)
          )
        }
      } else {
        sel_vars  <- get_selected_vars_mod(fit_type, ct, variables_for_fs_imputed)
        val_data  <- variables_for_fs_imputed_validation[, sel_vars, drop = FALSE]
        train_sub <- variables_for_fs_imputed[, sel_vars, drop = FALSE]
        if (fit_type == "glm_fit") {
          model <- retrain_glm(train_sub, targets_for_fs_imputed[[reg_target_varname]])
          y_new <- stats::predict(object = model, newdata = val_data)
        } else {
          model <- retrain_penalized_linear(fit_type, train_sub, targets_for_fs_imputed[[reg_target_varname]])
          y_new <- stats::predict(object = model$model, newx = as.matrix(val_data))
        }
      }

      y_obs         <- targets_for_fs_imputed_validation[[reg_target_varname]]
      nls_pred      <- nls(y_new ~ const + slope * y_obs, start = list(const = median(y_new), slope = 1))
      nls_pred_null <- nls(y_new ~ const + 0 * y_obs, start = list(const = median(y_new)))

      var_anaylsis  <- anova(nls_pred, nls_pred_null)
      full_model    <- summary(nls_pred)
      reduced_model <- summary(nls_pred_null)
      df_pred       <- cbind.data.frame(y_new, y_obs)

      p <- ggplot(df_pred, aes(x = y_obs, y = y_new)) +
        geom_point(color = actual_palette[4]) +
        geom_smooth(method = lm, color = "red", se = TRUE) +
        stat_poly_eq(
          formula = y ~ x,
          aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
          parse = TRUE
        ) +
        geom_hline(yintercept = median(df_pred$y_new), linetype = "dashed", color = "grey60") +
        theme_plot() +
        labs(title = paste0("Prediction of ", reg_target_varname, ": Model = ", fit_type))

      return(list(var_anaylsis = var_anaylsis, full_model = full_model, reduced_model = reduced_model,
                  df_pred = df_pred, pred_plot = p))
    })
    names(fits) <- c("glm_fit", "elastic", "lasso", "ridge")
    fits
  })
  names(result) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")
  result
})
names(reg_modulators_predictions) <- c("all", "reduced")

# ============================================================================ #
# 9. BORUTA IMPORTANCE VISUALIZATION
# ============================================================================ #

vars_to_update <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

for (v in vars_to_update) {
  reg_modulators_resultat[[v]]$penalized_regression_res$coef_table <-
    reg_modulators_resultat[[v]]$penalized_regression_res$coef_table %>%
    mutate(variable = variable %>%
             str_replace_all("`", "") %>%
             str_squish()) %>%
    left_join(
      enframe(reg_modulators_resultat[[v]]$Boruta_res$finalDecision,
              name = "variable",
              value = "Boruta") %>%
        mutate(variable = str_squish(variable)),
      by = "variable"
    )
}


selected_vars <- names(reg_modulators_resultat[[1]]$Boruta_res$finalDecision)[reg_modulators_resultat[[1]]$Boruta_res$finalDecision == "Confirmed"]
selected_vars


Boruta_data <- prepare_boruta_plot_data(boruta_res = reg_modulators_resultat$`AmmoLa intensity`$Boruta_res)
Boruta_data$importance <- Boruta_data$importance %>%
  group_by(Feature) %>%
  mutate(median = median(Importance[is.finite(Importance)], na.rm = TRUE)) %>%
  ungroup()


p_ammoLa_Boruta <- ggplot(Boruta_data$importance, aes(x = reorder(Feature, median), y = Importance, color = Color)) +
  geom_boxplot() +
  theme_plot() +
  labs(
    title = "AmmoLa rating: Relevant features according to 'Boruta' analysis",
    fill = "Decision", color = "Decision",
    x = "Features", y = "Importance"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    legend.position = c(.1, .8),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
    legend.key       = element_rect(fill = ggplot2::alpha("white", 0.5))
  ) +
  scale_color_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  ) +
  scale_fill_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  )

print(p_ammoLa_Boruta)
ggsave("p_ammoLa_Boruta.svg", p_ammoLa_Boruta,
       width = 25, height = 8, dpi = 300, limitsize = FALSE)
ggsave("p_ammoLa_Boruta.png", p_ammoLa_Boruta,
       width = 25, height = 8, dpi = 300, limitsize = FALSE)


# ============================================================================ #
# 10. RANDOM FOREST PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

# With all features

rf_modulators_predictions <- lapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  actual_data <-  variables_for_fs_imputed
  set.seed(42)

  rf_model <- randomForest(
    x = actual_data,
    y = targets_for_fs_imputed[[reg_target_varname]],
    mtry = reg_modulators_resultat[[reg_target_varname]]$rf_best_params$mtry,
    ntree = reg_modulators_resultat[[reg_target_varname]]$rf_best_params$ntree
  )

  y_new <- predict(rf_model, newdata = variables_for_fs_imputed_validation)

  y_obs <- targets_for_fs_imputed_validation[[reg_target_varname]]
  nls_pred <- nls(y_new ~ const + slope*y_obs, start = list(const = median(y_new), slope = 1))
  nls_pred_null  <- nls(y_new ~ const + 0 * y_obs, start = list(const = median(y_new)))

  var_anaylsis <- anova(nls_pred,nls_pred_null)
  full_model <- summary(nls_pred)
  reduced_model <- summary(nls_pred_null)

  df_pred <- cbind.data.frame(y_new, y_obs)

  p <- ggplot(df_pred, aes(x = y_obs, y = y_new)) +
    geom_point(color = actual_palette[4]) +
    geom_smooth(method=lm , color="red", se=TRUE) +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(after_stat(eq.label),
                                   after_stat(rr.label),
                                   sep = "~~~")),
                 parse = TRUE) +
    geom_hline(yintercept = median(df_pred$y_new), linetype = "dashed", color = "grey60") +
    theme_plot() +
    labs(title = paste0("Prediction of ", reg_target_varname, ": Model = ", "RF"))


  return(list(var_anaylsis = var_anaylsis, full_model = full_model, reduced_model = reduced_model, df_pred = df_pred, pred_plot = p))
})

names(rf_modulators_predictions) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

# With selected features only
rf_modulators_predictions_selected_features <- lapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  actual_data <-  variables_for_fs_imputed[
    names(reg_modulators_resultat[[reg_target_varname]]$Boruta_res$finalDecision[
      reg_modulators_resultat[[reg_target_varname]]$Boruta_res$finalDecision == "Confirmed"])]

  set.seed(42)

  # Quick tune RF for REGRESSION
  mtry_values <- unique(round(c(2, ncol(actual_data) / 3, ncol(actual_data) / 2, sqrt(ncol(actual_data)), log2(ncol(actual_data)), ncol(actual_data))))  # regression default is ncol/3, not sqrt
  ntree_values <- c(100, 200, 500, 1000)
  grid <- expand.grid(mtry = mtry_values, ntree = ntree_values)
  grid$median_mse <- NA
  mses_matrix <- matrix(NA, nrow = nrow(grid), ncol = 10)

  for (run in 1:10) {
    set.seed(42 + run)
    for (i in 1:nrow(grid)) {
      model <- randomForest(
        x = actual_data,
        y = targets_for_fs_imputed[[reg_target_varname]],
        mtry = grid$mtry[i],
        ntree = grid$ntree[i]
      )
      mses_matrix[i, run] <- model$mse[length(model$mse)]  # final OOB MSE
    }
  }

  grid$median_mse <- apply(mses_matrix, 1, median)
  best <- grid[which.min(grid$median_mse), ]

  set.seed(42)

  rf_model <- randomForest(
    x = actual_data,
    y = targets_for_fs_imputed[[reg_target_varname]],
    mtry = best[["mtry"]],
    ntree = best[["ntree"]]
  )

  y_new <- predict(rf_model, newdata = variables_for_fs_imputed_validation)

  y_obs <- targets_for_fs_imputed_validation[[reg_target_varname]]
  nls_pred <- nls(y_new ~ const + slope*y_obs, start = list(const = median(y_new), slope = 1))
  nls_pred_null  <- nls(y_new ~ const + 0 * y_obs, start = list(const = median(y_new)))

  var_anaylsis <- anova(nls_pred,nls_pred_null)
  full_model <- summary(nls_pred)
  reduced_model <- summary(nls_pred_null)

  df_pred <- cbind.data.frame(y_new, y_obs)

  p <- ggplot(df_pred, aes(x = y_obs, y = y_new)) +
    geom_point(color = actual_palette[4]) +
    geom_smooth(method=lm , color="red", se=TRUE) +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(after_stat(eq.label),
                                   after_stat(rr.label),
                                   sep = "~~~")),
                 parse = TRUE) +
    geom_hline(yintercept = median(df_pred$y_new), linetype = "dashed", color = "grey60") +
    theme_plot() +
    labs(title = paste0("Prediction of ", reg_target_varname, ": Model = ", "RF"))


  return(list(var_anaylsis = var_anaylsis, full_model = full_model, reduced_model = reduced_model, df_pred = df_pred, pred_plot = p))
})

names(rf_modulators_predictions_selected_features) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

rf_modulators_predictions$`AmmoLa intensity`$var_anaylsis
rf_modulators_predictions_selected_features$`AmmoLa intensity`$var_anaylsis


# ============================================================================ #
# 11. COMPILE AND SUMMARIZE RESULTS
# ============================================================================ #

# Extract and combine plots

plot_list_reg_modulators_psy <- lapply(
  lapply(reg_modulators_predictions[["all"]], "[[", "glm_fit"),
  "[[", "pred_plot"
)
# Add title above (cowplot style)
combined_plots_reg_modulators_psy <-
  cowplot::plot_grid(cowplot::plot_grid(
    ggplot() +
      labs(title = "Prediction of psychophysical trigeminal measures from candidate modulators of trigeminal sensitivity",
           subtitle = "Predicted versus measured using standard linear regression on new data") +
      theme(plot.title = element_text(size = 25, face = "plain", margin = ggplot2::margin(b = 4)),
            plot.subtitle = element_text(size = 20, face = "plain", margin = ggplot2::margin(b = 1))) +
      theme_void()
  ),
  plot_grid(
    plotlist = plot_list_reg_modulators_psy,
    nrow = 1),
  ncol = 1,
  rel_heights = c(0.1, 1)
  )

print(combined_plots_reg_modulators_psy)


# Define the penalized regression types
penal_types <- c("glm_fit", "elastic", "lasso", "ridge")

# Build the long-format summary table
make_validation_table_modulators <- function() {
  outcomes   <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")
  reg_models <- c("glm_fit", "elastic", "lasso", "ridge")

  rows <- lapply(c("all", "reduced"), function(fs) {
    lapply(outcomes, function(outcome) {
      reg_rows <- lapply(reg_models, function(model) {
        res     <- reg_modulators_predictions[[fs]][[outcome]][[model]]
        anova_f <- round(res$var_anaylsis$F[2], 3)
        anova_p <- signif(res$var_anaylsis$`Pr(>F)`[2], 3)
        slope   <- round(res$full_model$coefficients["slope", "Estimate"], 3)
        slope_p <- signif(res$full_model$coefficients["slope", "Pr(>|t|)"], 3)
        data.frame(Feature_set = fs, Outcome = outcome, Model = model,
                   ANOVA_F = anova_f, ANOVA_p = anova_p, Slope = slope, Slope_p = slope_p)
      })

      rf_res  <- if (fs == "all") rf_modulators_predictions[[outcome]] else rf_modulators_predictions_selected_features[[outcome]]
      anova_f <- round(rf_res$var_anaylsis$F[2], 3)
      anova_p <- signif(rf_res$var_anaylsis$`Pr(>F)`[2], 3)
      slope   <- round(rf_res$full_model$coefficients["slope", "Estimate"], 3)
      slope_p <- signif(rf_res$full_model$coefficients["slope", "Pr(>|t|)"], 3)
      rf_row  <- data.frame(Feature_set = fs, Outcome = outcome, Model = "RF",
                            ANOVA_F = anova_f, ANOVA_p = anova_p, Slope = slope, Slope_p = slope_p)

      do.call(rbind, c(reg_rows, list(rf_row)))
    })
  })

  do.call(rbind, do.call(c, rows))
}

validation_table_modulators <- make_validation_table_modulators()

# Clean up model names for publication
validation_table_modulators$Model <- recode(validation_table_modulators$Model,
  "glm_fit" = "OLS",
  "elastic" = "Elastic net",
  "lasso"   = "LASSO",
  "ridge"   = "Ridge",
  "RF"      = "Random forest"
)

validation_table_modulators$Feature_set <- recode(validation_table_modulators$Feature_set,
  "all"     = "All features",
  "reduced" = "Method-selected features"
)

# Clean up outcome names
validation_table_modulators$Outcome <- recode(validation_table_modulators$Outcome,
  "Lateralization (x/20)" = "Lateralization",
  "AmmoLa intensity"      = "AmmoLa intensity",
  "CO2 threshold"         = "CO\u2082 threshold"
)

print(validation_table_modulators)

# ============================================================================ #
# 13. ANALYZE SIGNIFICANT MODULATORS
# ============================================================================ #

outcomes <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")
# Combine penalized regression coefficient tables across all outcomes
coef_tables_combined <- dplyr::bind_rows(
  lapply(outcomes, function(o) {
    reg_modulators_resultat[[o]]$penalized_regression_res$coef_table %>%
      dplyr::mutate(outcome = o, .before = 1)
  })
)

sig_outcome <- validation_table_modulators[
  validation_table_modulators$Feature_set == "All features" &
  validation_table_modulators$ANOVA_p < 0.05 &
  validation_table_modulators$Slope_p  < 0.05 &
  validation_table_modulators$Slope    > 0,
  c("Outcome", "Model")
]
# Reverse the outcome recode so names match coef_tables_combined
sig_outcome$Outcome <- dplyr::recode(sig_outcome$Outcome,
  "Lateralization"      = "Lateralization (x/20)",
  "CO\u2082 threshold"  = "CO2 threshold"
)
names(sig_outcome) <- c("outcome", "Model")

filtered_rows <- data.frame()

filtered_rows <- coef_tables_combined %>%
  inner_join(sig_outcome, by = "outcome") %>%
  filter(
    (Model == "OLS" & glm_p < 0.05) |
      (Model == "Ridge" & ridge_selected) |
      (Model == "LASSO" & lasso_selected) |
      (Model == "Elastic net" & elastic_selected) |
      (Model == "Random forest" & Boruta == "Confirmed")
  )

filtered_rows$Boruta <- ifelse(filtered_rows$Boruta == "Confirmed", TRUE, FALSE)
filtered_rows[,1:2]

# Statistical tests for confirmed modulators (Wilcoxon for binary, Spearman for continuous)
unique_pairs <- unique(filtered_rows[, c("outcome", "variable")])

stats_results <- lapply(seq_len(nrow(unique_pairs)), function(i) {
  outcome_name  <- unique_pairs$outcome[i]
  variable_name <- unique_pairs$variable[i]

  y <- targets_for_fs_imputed[[outcome_name]]
  x <- variables_for_fs_imputed[[variable_name]]

  n_unique <- length(unique(na.omit(x)))

  if (n_unique <= 2) {
    test <- wilcox.test(y ~ x, exact = FALSE)
    data.frame(
      outcome   = outcome_name,
      variable  = variable_name,
      test      = "Wilcoxon",
      statistic = unname(test$statistic),
      p_value   = test$p.value,
      signif = symnum(test$p.value, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols   = c("***", "**", "*", ".", "")),
      effsize = rank_biserial_val(y = y, x = x)
    )
  } else {
    test <- cor.test(y, x, method = "spearman", exact = FALSE)
    data.frame(
      outcome   = outcome_name,
      variable  = variable_name,
      test      = "Spearman",
      statistic = unname(test$estimate),
      p_value   = test$p.value,
      signif = symnum(test$p.value, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols   = c("***", "**", "*", ".", "")),
      effsize = unname(test$estimate)
    )
  }
})

stats_results_df <- dplyr::bind_rows(stats_results)
print(stats_results_df)


stats_results_validation <- lapply(seq_len(nrow(unique_pairs)), function(i) {
  outcome_name  <- unique_pairs$outcome[i]
  variable_name <- unique_pairs$variable[i]

  y <- targets_for_fs_imputed_validation[[outcome_name]]
  x <- variables_for_fs_imputed_validation[[variable_name]]

  n_unique <- length(unique(na.omit(x)))

  if (n_unique <= 2) {
    if (n_unique == 2) {
    test <- wilcox.test(y ~ x, exact = FALSE)
    data.frame(
      outcome   = outcome_name,
      variable  = variable_name,
      test      = "Wilcoxon",
      statistic = unname(test$statistic),
      p_value   = test$p.value,
      signif = symnum(test$p.value, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols   = c("***", "**", "*", ".", "")),
      effsize = rank_biserial_val(y = y, x = x)
    ) }
  } else {
    test <- cor.test(y, x, method = "spearman", exact = FALSE)
    data.frame(
      outcome   = outcome_name,
      variable  = variable_name,
      test      = "Spearman",
      statistic = unname(test$estimate),
      p_value   = test$p.value,
      signif = symnum(test$p.value, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols   = c("***", "**", "*", ".", "")),
      effsize = unname(test$estimate)
    )
  }
})

stats_results_validation_df <- dplyr::bind_rows(stats_results_validation)
print(stats_results_validation_df)

stats_results_all <- rbind.data.frame(cbind.data.frame(Data = "Training", stats_results_df),
                                      cbind.data.frame(Data = "Validation", stats_results_validation_df))

joined_df <- stats_results_df %>%
  left_join(
    stats_results_validation_df,
    by = c("outcome", "variable"),
    suffix = c("_train", "_validation")
  )

joined_df$mean_effect <- rowMeans(joined_df[, c("effsize_validation", "effsize_train")], na.rm = TRUE)

cor.test(joined_df$effsize_validation, joined_df$effsize_train, method = "kendall")

joined_df_not1 <- joined_df[abs(joined_df$mean_effect) < 0.99,]
cor.test(joined_df_not1$effsize_validation, joined_df_not1$effsize_train, method = "kendall")


set.seed(42)
p_modulator_cor_Effsizes_tau <-
  ggplot(joined_df, aes(x = effsize_train, y = effsize_validation, color = as.factor(mean_effect))) +
  geom_point(size = 4, show.legend = FALSE) +
  ggpubr::stat_cor(data = joined_df, aes(x = effsize_train, y = effsize_validation, color = mean_effect), method = "kendall", cor.coef.name = "tau", label.x.npc = "left", label.y.npc = "top") +
  ggrepel::geom_text_repel(aes(x = effsize_train, y = effsize_validation, label = variable), force = TRUE, force_pull = TRUE, size = 1.5, show.legend = FALSE) +
  guides(color = "none") +
  theme_plot()




# ============================================================================ #
# 12. SEX DIFFERENCES STATISTICS # FOR PURE EXPLORATION
# ============================================================================ #

# Training set
# x = your 0/1 vector, same length as nrow(df)
# df = data frame with numeric columns

res_sex_training <- lapply(targets_for_fs_imputed, function(col) {
  wilcox.test(col ~ as.factor(variables_for_fs_imputed$Gender_w))
})
res_sex_validation<- lapply(targets_for_fs_imputed_validation, function(col) {
  wilcox.test(col ~ as.factor(variables_for_fs_imputed_validation$Gender_w))
})

res_sex_training
res_sex_validation

# ============================================================================ #
# 13. SAVE RESULTS
# ============================================================================ #

write.csv(coef_tables_combined, "regression_modulators_coef_tables_combined.csv", row.names = FALSE)

write.csv(validation_table_modulators, "regression_modulators_final_table_modulators.csv", row.names = FALSE)

write.csv(stats_results_all, "regression_modulators_sig_features_stats.csv")

cat("Saved: regression_modulators_coef_tables_combined.csv\n")
cat("Saved: regression_modulators_final_table_modulators.csv\n")

ggsave("combined_plots_reg_modulators_psy.svg", combined_plots_reg_modulators_psy,
       width = 15, height = 6, dpi = 300, limitsize = FALSE)
ggsave("combined_plots_reg_modulators_psy.png", combined_plots_reg_modulators_psy,
       width = 15, height = 6, dpi = 300, limitsize = FALSE)


# ============================================================================ #
# END OF SCRIPT
# ============================================================================ #
