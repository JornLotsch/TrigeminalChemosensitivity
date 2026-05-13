################################################################################
# Trigeminal Sensitivity Analysis - Within-Trigeminal Regression
#
# Author: Jorn Lotsch
# Date: 2026-04-29
#
# Description:
#   Regression analyses predicting each of three trigeminal psychophysical
#   outcomes (Lateralization, AmmoLa intensity, CO2 threshold) from the
#   remaining variables in the trigeminal dataspace (nasal chemosensory
#   perception and the other two psychophysical trigeminal measures). Combines
#   Boruta feature selection with penalized regression (ridge, lasso, elastic
#   net) and Random Forest regression. Models are trained on an imputed training
#   dataset and validated on an independent dataset.
#
# Input Files:
#   - analysis_dataset_training_imputed.csv: Training data with imputed values
#   - analysis_dataset_validation_imputed.csv: Validation data with imputed values
#   - globals.R: Global variables, color palettes, and transformation functions
#   - utils.R: Utility functions
#
# Output Files:
#   - regression_coef_tables_combined.csv: Coefficients for all models and outcomes
#   - regression_final_table.csv: Validation performance (slope, ANOVA F) per model
#   - combined_plots_reg_psy.svg/.png: Predicted vs. observed plots
#
# Key Analyses:
#   1. Data loading, variable selection, and transformation
#   2. Boruta feature selection (Random Forest-based importance)
#   3. Penalized linear regression (ridge, lasso, elastic net via glmnet)
#   4. Standard GLM for coefficient p-values
#   5. RF regression with tuned hyperparameters (mtry, ntree)
#   6. Model validation: predicted vs. observed on independent dataset
#   7. Compilation of variable selection and prediction results
#
# Statistical Methods:
#   - Boruta algorithm (RF-based feature selection)
#   - Penalized regression: ridge (alpha=0), lasso (alpha=1), elastic net (alpha=0.5)
#   - Cross-validation for lambda tuning (cv.glmnet, 5-fold)
#   - NLS for prediction assessment; ANOVA for model comparison
################################################################################

# ============================================================================ #\n# PART 2: MACHINE LEARNING FRAMEWORK — Psychophysical Regression
# ============================================================================ #\n# Regression analyses predicting each psychophysical outcome from remaining\n# trigeminal variables (20 predictors per outcome).\n#\n# Model Development: Training set only (n=800)\n# Validation: Independent held-out set (n=201)\n#\n# Regression Methods: OLS, Ridge, Lasso, Elastic Net (α=0.5), Random Forest\n# Hyperparameter Tuning: 5-fold cross-validation for λ; grid search for RF\n# Performance Assessment: Nonlinear LS (predicted vs. observed) on validation set\n#\n# See README.md \"PART 2\" > \"Prediction of psychophysical measures (from trigeminal variables)\"\n# and manuscript Methods for full technical details\n# ============================================================================ #

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

# ============================================================================ #
# 4. LOAD AND PREPARE DATA
# ============================================================================ #

cat("\n=== Loading Data ===\n")

# Load training data (imputed)
trigeminale_training_data <- read.csv("analysis_dataset_training_imputed.csv",
                                      check.names = FALSE)
cat("Training data loaded:", nrow(trigeminale_training_data), "samples\n")

# Extract variables for regression
variables_for_regression_imputed <- trigeminale_training_data[, c(variables_by_categories$Nasal_chemosensory_perception,
                                                                variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
rownames(variables_for_regression_imputed) <- trigeminale_training_data$ID

# Load validation data (imputed)
trigeminale_validation_data <- read.csv("analysis_dataset_validation_imputed.csv",
                                        check.names = FALSE)
cat("Validation data loaded:", nrow(trigeminale_validation_data), "samples\n")

# Extract variables for regression (validation)
variables_for_regression_imputed_validation <- trigeminale_validation_data[, c(variables_by_categories$Nasal_chemosensory_perception,
                                                                             variables_by_categories$Psychophysical_measurements[c(1, 2, 4)])]
rownames(variables_for_regression_imputed_validation) <- trigeminale_validation_data$ID

# ============================================================================ #
# 5. DATA TRANSFORMATION AND PREPROCESSING
# ============================================================================ #

cat("\n=== Transforming Variables ===\n")

# Backup original variable names
original_names <- names(variables_for_regression_imputed)

# Transform psychophysical variables (training)
# AmmoLa intensity: reflected slog (handles ceiling effects)
variables_for_regression_imputed$`AmmoLa intensity` <-
  reflect_slog_unflipped(variables_for_regression_imputed$`AmmoLa intensity`)

# CO2 threshold: negative slog (lower = more sensitive)
variables_for_regression_imputed$`CO2 threshold` <-
  -slog(variables_for_regression_imputed$`CO2 threshold`)

# Lateralization: none
variables_for_regression_imputed$`Lateralization (x/20)` <-
  variables_for_regression_imputed$`Lateralization (x/20)`

# Transform psychophysical variables (validation)
variables_for_regression_imputed_validation$`AmmoLa intensity` <-
  reflect_slog_unflipped(variables_for_regression_imputed_validation$`AmmoLa intensity`)

variables_for_regression_imputed_validation$`CO2 threshold` <-
  -slog(variables_for_regression_imputed_validation$`CO2 threshold`)

variables_for_regression_imputed_validation$`Lateralization (x/20)` <-
  variables_for_regression_imputed_validation$`Lateralization (x/20)`

cat("Variables transformed to [0, 3] scale\n")


# Inspect dataset structure
cat("\nDataset structure:\n")
print(names(variables_for_regression_imputed))

# ============================================================================ #
# 6. BORUTA FEATURE SELECTION AND PENALIZED REGRESSION
# ============================================================================ #

reg_resultat <- pbmcapply::pbmclapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  actual_data <-  variables_for_regression_imputed[,!names(variables_for_regression_imputed) %in% reg_target_varname]
  names(actual_data) <- make.names(names(actual_data))
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
        y = variables_for_regression_imputed[[reg_target_varname]],
        mtry = grid$mtry[i],
        ntree = grid$ntree[i]
      )
      mses_matrix[i, run] <- model$mse[length(model$mse)]  # final OOB MSE
    }
  }

  grid$median_mse <- apply(mses_matrix, 1, median)
  best <- grid[which.min(grid$median_mse), ]

  # bestmtry <- randomForest::tuneRF(x = actual_data,
  #                    y = variables_for_regression_imputed_1[[reg_target_varname]], improve=1e-5, ntree=500)
  # best_mtry <- bestmtry[which.min(bestmtry[, 2]), 1]

  Boruta_res <- Boruta::Boruta(x = actual_data,
                               y = variables_for_regression_imputed[[reg_target_varname]], ntree = best[["ntree"]], mtry = best[["mtry"]], maxRuns = 1000, doTrace = TRUE)

  penalized_regression_res <- run_penalized_regression_all(
    train_data   = actual_data,
    train_target = variables_for_regression_imputed[[reg_target_varname]],
    alpha_elastic = 0.5,
    nfolds = 5,
    ridge_threshold = 0.05  # or something small that makes sense in your scale
  )
  return(list(Boruta_res = Boruta_res, rf_best_params = best, penalized_regression_res = penalized_regression_res))
}, mc.cores = nProc_desired)

names(reg_resultat) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

# ============================================================================ #
# 7. INITIAL RESULT INSPECTION
# ============================================================================ #

par(mfrow=c(2,3))
plot(reg_resultat[[1]]$Boruta_res, las = 2)
plot(reg_resultat[[2]]$Boruta_res, las = 2)
plot(reg_resultat[[3]]$Boruta_res, las = 2)
par(mfrow=c(1,1))

print(reg_resultat[["Lateralization (x/20)"]]$penalized_regression_res)
print(reg_resultat[["AmmoLa intensity"]]$penalized_regression_res)
print(reg_resultat[["CO2 threshold"]]$penalized_regression_res)

# ============================================================================ #
# 8. REGRESSION MODEL PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

reg_predictions <- lapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  print(reg_target_varname)
  fits <- lapply(c( "glm_fit", "elastic", "lasso", "ridge"), function(fit_type) {
    print(fit_type)

    if(fit_type =="glm_fit") {
      y_new <-stats::predict(object = reg_resultat[[reg_target_varname]]$penalized_regression_res$glm_fit,
                             newdata = (variables_for_regression_imputed_validation[,!names(variables_for_regression_imputed_validation) %in% reg_target_varname]))
    } else {
      y_new <-stats::predict(object = reg_resultat[[reg_target_varname]]$penalized_regression_res[[fit_type]]$model,
                             newx = as.matrix(variables_for_regression_imputed_validation[,!names(variables_for_regression_imputed_validation) %in% reg_target_varname]))
    }
    y_obs <- variables_for_regression_imputed_validation[[reg_target_varname]]
    nls_pred <- nls(y_new ~ const + slope*y_obs, start = list(const = median(y_new), slope = 1))
    nls_pred_null  <- nls(y_new ~ const + 0 * y_obs, start = list(const = median(y_new)))

    var_anaylsis <- anova(nls_pred,nls_pred_null)
    full_model <- summary(nls_pred)
    reduced_model <- summary(nls_pred_null)

    df_pred <- cbind.data.frame(y_new, y_obs)
    rng <- range(c(df_pred$y_obs, df_pred$y_new), na.rm = TRUE)

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
      labs(title = paste0("Prediction of ", reg_target_varname, ": Model = ", fit_type))


    return(list(var_anaylsis = var_anaylsis, full_model = full_model, reduced_model = reduced_model, df_pred = df_pred, pred_plot = p))
  })
  names(fits) <- c("glm_fit", "elastic", "lasso", "ridge")
  return(fits)
})

names(reg_predictions) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

# ============================================================================ #
# 9. BORUTA IMPORTANCE VISUALIZATION
# ============================================================================ #

vars_to_update <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

for (v in vars_to_update) {
  reg_resultat[[v]]$penalized_regression_res$coef_table <-
    reg_resultat[[v]]$penalized_regression_res$coef_table %>%
    mutate(variable = variable %>%
             str_replace_all("`", "") %>%
             str_squish()) %>%
    left_join(
      enframe(reg_resultat[[v]]$Boruta_res$finalDecision,
              name = "variable",
              value = "Boruta") %>%
        mutate(variable = str_squish(variable)),
      by = "variable"
    )
}


selected_vars <- names(reg_resultat[[1]]$Boruta_res$finalDecision)[reg_resultat[[1]]$Boruta_res$finalDecision == "Confirmed"]
selected_vars


Boruta_data <- prepare_boruta_plot_data(boruta_res = reg_resultat[[1]]$Boruta_res)
Boruta_data$importance <- Boruta_data$importance %>%
  group_by(Feature) %>%
  mutate(median = median(Importance[is.finite(Importance)], na.rm = TRUE)) %>%
  ungroup()


ggplot(Boruta_data$importance, aes(x = reorder(Feature, median), y = Importance, color = Color)) +
  geom_boxplot() +
  theme_light() +
  labs(
    title = "title",
    fill = "Decision", color = "Decision",
    x = "Features", y = "Importance"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    legend.position = c(.1, .8),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
    legend.key       = element_rect(fill = ggplot2::alpha("white", 0.5)),
    plot.title = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  ) +
  scale_fill_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  )

# ============================================================================ #
# 10. RANDOM FOREST PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

rf_predictions <- lapply(c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold"), function(reg_target_varname) {
  actual_data <-  variables_for_regression_imputed[,!names(variables_for_regression_imputed) %in% reg_target_varname]
  names(actual_data) <- make.names(names(actual_data))
  set.seed(42)

  rf_model <- randomForest(
    x = actual_data,
    y = variables_for_regression_imputed[[reg_target_varname]],
    mtry = reg_resultat[[reg_target_varname]]$rf_best_params$mtry,
    ntree = reg_resultat[[reg_target_varname]]$rf_best_params$ntree
  )

  y_new <- predict(rf_model, newdata = variables_for_regression_imputed_validation)

  y_obs <- variables_for_regression_imputed_validation[[reg_target_varname]]
  nls_pred <- nls(y_new ~ const + slope*y_obs, start = list(const = median(y_new), slope = 1))
  nls_pred_null  <- nls(y_new ~ const + 0 * y_obs, start = list(const = median(y_new)))

  var_anaylsis <- anova(nls_pred,nls_pred_null)
  full_model <- summary(nls_pred)
  reduced_model <- summary(nls_pred_null)

  df_pred <- cbind.data.frame(y_new, y_obs)
  rng <- range(c(df_pred$y_obs, df_pred$y_new), na.rm = TRUE)

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

names(rf_predictions) <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")

rf_predictions$`CO2 threshold`$full_model$coefficients

# ============================================================================ #
# 11. COMPILE AND SUMMARIZE RESULTS
# ============================================================================ #

# Extract and combine plots

plot_list_reg_psy <- lapply(lapply(reg_predictions, "[[", 1), "[[", "pred_plot")
# Add title above (cowplot style)
combined_plots_reg_psy <-
  cowplot::plot_grid(cowplot::plot_grid(
    ggplot() +
      labs(title = "Prediction of psychophysical trigeminal measures from the respective other trigeminal measures",
           subtitle = "Predicted versus measured using standard linear regression on new data") +
      theme(plot.title = element_text(size = 25, face = "plain", margin = ggplot2::margin(b = 4)),
            plot.subtitle = element_text(size = 20, face = "plain", margin = ggplot2::margin(b = 1))) +
      theme_void()
  ),
  plot_grid(
    plotlist = plot_list_reg_psy,
    nrow = 1),
  ncol = 1,
  rel_heights = c(0.1, 1)
  )

print(combined_plots_reg_psy)


# Build the long-format summary table
make_validation_table <- function(reg_predictions) {

  outcomes <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")
  models <- c("glm_fit", "elastic", "lasso", "ridge", "RF")

  rows <- lapply(outcomes, function(outcome) {
    lapply(models, function(model) {

      if (model == "RF") {
        res <- rf_predictions[[which(outcomes == outcome)]]
      } else {
        res <- reg_predictions[[which(outcomes == outcome)]][[model]]
      }

      anova_f   <- round(res$var_anaylsis$F[2], 3)
      anova_p   <- signif(res$var_anaylsis$`Pr(>F)`[2], 3)
      slope     <- round(res$full_model$coefficients["slope", "Estimate"], 3)
      slope_p   <- signif(res$full_model$coefficients["slope", "Pr(>|t|)"], 3)

      data.frame(
        Outcome = outcome,
        Model   = model,
        ANOVA_F = anova_f,
        ANOVA_p = anova_p,
        Slope   = slope,
        Slope_p = slope_p
      )
    })
  })

  do.call(rbind, do.call(c, rows))
}

validation_table <- make_validation_table(reg_predictions)

# Clean up model names for publication
validation_table$Model <- recode(validation_table$Model,
                                 "glm_fit"  = "OLS",
                                 "elastic" = "Elastic net",
                                 "lasso"    = "LASSO",
                                 "ridge"    = "Ridge",
                                 "RF"       = "Random forest"
)

# Clean up outcome names
validation_table$Outcome <- recode(validation_table$Outcome,
                                   "Lateralization (x/20)" = "Lateralization",
                                   "AmmoLa intensity"      = "AmmoLa intensity",
                                   "CO2 threshold"         = "CO\u2082 threshold"
)

print(validation_table)

# ============================================================================ #
# 12. SAVE RESULTS
# ============================================================================ #

# Combine penalized regression coefficient tables across all outcomes
outcomes <- c("Lateralization (x/20)", "AmmoLa intensity", "CO2 threshold")
coef_tables_combined <- dplyr::bind_rows(
  lapply(outcomes, function(o) {
    reg_resultat[[o]]$penalized_regression_res$coef_table %>%
      dplyr::mutate(outcome = o, .before = 1)
  })
)
write.csv(coef_tables_combined, "regression_coef_tables_combined.csv", row.names = FALSE)

write.csv(validation_table, "regression_final_table.csv", row.names = FALSE)

cat("Saved: regression_coef_tables_combined.csv\n")
cat("Saved: regression_final_table.csv\n")

ggsave("combined_plots_reg_psy.svg", combined_plots_reg_psy,
       width = 15, height = 6, dpi = 300, limitsize = FALSE)
ggsave("combined_plots_reg_psy.png", combined_plots_reg_psy,
       width = 15, height = 6, dpi = 300, limitsize = FALSE)


# ============================================================================ #
# END OF SCRIPT
# ============================================================================ #
