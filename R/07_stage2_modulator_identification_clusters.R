################################################################################
# Trigeminal Cluster Classification - Modulators Analysis
#
# Author: Jorn Lotsch
# Date: 2026-04-29
#
# Description:
#   Classification of trigeminal cluster membership using the modulators dataset.
#   Combines Boruta feature selection with penalized multinomial regression
#   (ridge, lasso, elastic net) and Random Forest classification. Models are
#   trained on the imputed training dataset and evaluated on the validation set.
#
# Input Files:
#   - trigeminal_clustered_data.csv
#   - analysis_dataset_training_imputed.csv
#   - analysis_dataset_validation_imputed.csv
#   - globals.R
#   - utils.R
#
# Output Files:
#   - clusters_trig_modulators_sig_predictors.csv
#   - clusters_trig_modulators_classification_reg.csv
#   - clusters_trig_modulators_classification_rf.csv
#   - clusters_trig_modulators_coef_table_full.csv
#   - p_clusters_trig_modulators_Boruta.svg
#   - p_clusters_trig_modulators_Boruta.png
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
        "Rejected", "Shadow"
      )
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

get_selected_vars_reg <- function(fit_type, coef_table, ref_data) {
  col <- switch(fit_type,
    "multinom_fit" = "glm_selected",
    "elastic"      = "elastic_selected",
    "lasso"        = "lasso_selected",
    "ridge"        = "ridge_selected"
  )
  mm_vars <- unique(coef_table$variable[coef_table[[col]] == TRUE])
  mm_vars <- trimws(gsub("`", "", mm_vars))
  names(ref_data)[
    names(ref_data) %in% mm_vars |
    make.names(names(ref_data)) %in% make.names(mm_vars)
  ]
}

retrain_multinom <- function(train_data, train_target) {
  df <- train_data
  df$target <- as.factor(train_target)
  nnet::multinom(target ~ ., data = df, trace = FALSE)
}

retrain_penalized_multinom <- function(fit_type, train_data, train_target, nfolds = 5, seed = 42) {
  alpha_val <- switch(fit_type, "elastic" = 0.5, "lasso" = 1, "ridge" = 0)
  y <- as.factor(train_target)
  X <- model.matrix(~ ., data = train_data)[, -1, drop = FALSE]
  set.seed(seed)
  cv_fit      <- glmnet::cv.glmnet(x = X, y = y, family = "multinomial", alpha = alpha_val, nfolds = nfolds)
  final_model <- glmnet::glmnet(x = X, y = y, family = "multinomial", alpha = alpha_val, lambda = cv_fit$lambda.min)
  list(model = final_model, lambda = cv_fit$lambda.min)
}

# ============================================================================ #
# 4. LOAD AND PREPARE DATA
# ============================================================================ #
cat("\n=== Loading Data ===\n")

# Read clustered data
all_trig_clustered_data_and_clusters <- read.csv(
  "trigeminal_clustered_data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

all_trig_clustered_data_and_clusters$ID <- as.character(all_trig_clustered_data_and_clusters$ID)
all_trig_clustered_data_and_clusters$Set <- as.character(all_trig_clustered_data_and_clusters$Set)
all_trig_clustered_data_and_clusters$Cluster <- as.factor(all_trig_clustered_data_and_clusters$Cluster)

trig_clusters_training_data <- subset(all_trig_clustered_data_and_clusters, Set == "Training")
trig_clusters_validation_data <- subset(all_trig_clustered_data_and_clusters, Set == "Validation")

# Read imputed datasets
trigeminale_training_data <- read.csv(
  "analysis_dataset_training_imputed.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

trigeminale_validation_data <- read.csv(
  "analysis_dataset_validation_imputed.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Ensure IDs are character
trigeminale_training_data$ID <- as.character(trigeminale_training_data$ID)
trigeminale_validation_data$ID <- as.character(trigeminale_validation_data$ID)

# Match imputed data to clustered data by ID
idx_train <- match(trig_clusters_training_data$ID, trigeminale_training_data$ID)
idx_valid <- match(trig_clusters_validation_data$ID, trigeminale_validation_data$ID)

if (any(is.na(idx_train))) stop("Some training IDs from clustered data were not found in imputed training data.")
if (any(is.na(idx_valid))) stop("Some validation IDs from clustered data were not found in imputed validation data.")

trigeminale_training_data <- trigeminale_training_data[idx_train, ]
trigeminale_validation_data <- trigeminale_validation_data[idx_valid, ]

# Build predictor sets by removing ID and any variables you do not want as predictors
drop_vars <- unique(c(
  "ID",
  variables_by_categories$Nasal_chemosensory_perception,
  variables_by_categories$Psychophysical_measurements[c(1, 2, 4)]
))

trig_train_x <- trigeminale_training_data[, !names(trigeminale_training_data) %in% drop_vars, drop = FALSE]
trig_valid_x <- trigeminale_validation_data[, !names(trigeminale_validation_data) %in% drop_vars, drop = FALSE]

# Add the target variable from clustered data
trig_train_final <- cbind(
  ID = trig_clusters_training_data$ID,
  Cluster = factor(trig_clusters_training_data$Cluster),
  trig_train_x
)

trig_valid_final <- cbind(
  ID = trig_clusters_validation_data$ID,
  Cluster = factor(trig_clusters_validation_data$Cluster),
  trig_valid_x
)

# Optional checks
stopifnot(identical(trig_train_final$ID, trig_clusters_training_data$ID))
stopifnot(identical(trig_valid_final$ID, trig_clusters_validation_data$ID))

# After matching rows by ID
rownames(trig_train_final) <- trig_train_final$ID
rownames(trig_valid_final) <- trig_valid_final$ID

# Remove ID columns
trig_train_final$ID <- NULL
trig_valid_final$ID <- NULL


# ============================================================================ #
# 6. BORUTA FEATURE SELECTION AND PENALIZED REGRESSION
# ============================================================================ #

clusters_trig_modulators_resultat <- pbmcapply::pbmclapply(c("regression", "RF"), function(method) {
  rf_best_params <- NA

  actual_data <- trig_train_final[, !names(trig_train_final) %in% c("Cluster")]

  switch(method,
    "RF" = {
      # Quick tune RF for CLASSIFICATION
      mtry_values <- unique(c(1, floor(sqrt(ncol(actual_data))), ceiling(ncol(actual_data) / 3), floor(ncol(actual_data) / 2)))
      ntree_values <- c(100, 200, 500, 1000)

      grid <- expand.grid(mtry = mtry_values, ntree = ntree_values)
      grid$error <- NA
      errors_matrix <- matrix(NA, nrow = nrow(grid), ncol = 10)
      for (run in 1:10) {
        set.seed(42 + run)
        for (i in 1:nrow(grid)) {
          model <- randomForest(
            x = actual_data,
            y = as.factor(clusters_training),
            mtry = grid$mtry[i],
            ntree = grid$ntree[i]
          )
          errors_matrix[i, run] <- mean(model$err.rate[, 1])
        }
      }
      grid$median_error <- apply(errors_matrix, 1, median)
      best <- grid[which.min(grid$median_error), ]
      rf_best_params <- best

      set.seed(42)
      res_classification <- Boruta::Boruta(
        x = actual_data,
        y = as.factor(clusters_training), ntree = best[["ntree"]], mtry = best[["mtry"]], maxRuns = 1000, doTrace = TRUE
      )
    },
    "regression" = {
      res_classification <- run_penalized_multinomial_all(
        train_data = actual_data,
        train_target = as.factor(clusters_training),
        alpha_elastic = 0.5,
        nfolds = 5,
        ridge_threshold = 0.05 # or something small that makes sense in your scale
      )
    }
  )

  return(list(res_classification = res_classification, rf_best_params = rf_best_params, actual_data = actual_data))
}, mc.cores = nProc_desired)

names(clusters_trig_modulators_resultat) <- c("regression", "RF")

# ============================================================================ #
# 7. INITIAL RESULT INSPECTION
# ============================================================================ #

par(mfrow = c(2, 1))
plot(clusters_trig_modulators_resultat$RF$res_classification, las = 2)
par(mfrow = c(1, 1))

print(clusters_trig_modulators_resultat$regression$res_classification)

# ============================================================================ #
# 8. CLASSIFCATION MODEL PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

clusters_trig_modulators_predictions_reg <- lapply(c("all", "reduced"), function(feature_set) {
  ct         <- clusters_trig_modulators_resultat$regression$res_classification$coef_table
  train_full <- trig_train_final[, !names(trig_train_final) %in% c("Cluster")]
  val_full   <- trig_valid_final[, !names(trig_valid_final) %in% c("Cluster")]

  fits <- lapply(c("multinom_fit", "elastic", "lasso", "ridge"), function(fit_type) {
    print(paste(feature_set, fit_type))

    if (feature_set == "all") {
      val_data <- val_full
      model    <- clusters_trig_modulators_resultat$regression$res_classification[[fit_type]]
    } else {
      sel_vars  <- get_selected_vars_reg(fit_type, ct, train_full)
      val_data  <- val_full[, sel_vars, drop = FALSE]
      train_sub <- train_full[, sel_vars, drop = FALSE]
      model <- if (fit_type == "multinom_fit") {
        retrain_multinom(train_sub, clusters_training)
      } else {
        retrain_penalized_multinom(fit_type, train_sub, clusters_training)
      }
    }

    reg_repeated <- lapply(1:100, function(s) {
      set.seed(41 + s)
      rows_to_sample <- sample(nrow(val_data), replace = TRUE)
      sampled_df     <- val_data[rows_to_sample, ]
      sampled_cls    <- clusters_validation[rows_to_sample]

      if (fit_type == "multinom_fit") {
        y_new <- stats::predict(object = model, newdata = sampled_df)
      } else {
        y_new <- as.factor(as.numeric(stats::predict(
          object = model$model,
          newx   = as.matrix(sampled_df),
          type   = "class"
        )))
      }
      y_obs  <- sampled_cls
      cm_raw <- caret::confusionMatrix(
        factor(y_new, levels = unique(clusters_training)),
        factor(y_obs, levels = unique(clusters_training))
      )$byClass
      cm_reg <- if (length(unique(clusters_training)) > 2) colMeans(cm_raw, na.rm = TRUE) else cm_raw
      list(cm_reg = cm_reg, df_pred = cbind.data.frame(y_new, y_obs))
    })

    cm_reg_rep  <- do.call(rbind, lapply(reg_repeated, function(x) x$cm_reg))
    df_pred_rep <- do.call(cbind, lapply(reg_repeated, function(x) x$df_pred))
    list(cm_reg_rep = cm_reg_rep, df_pred_rep = df_pred_rep)
  })
  names(fits) <- c("multinom_fit", "elastic", "lasso", "ridge")
  fits
})
names(clusters_trig_modulators_predictions_reg) <- c("all", "reduced")

# ============================================================================ #
# 9. BORUTA IMPORTANCE VISUALIZATION
# ============================================================================ #

clusters_trig_modulators_resultat$regression$res_classification$coef_table <-
  clusters_trig_modulators_resultat$regression$res_classification$coef_table %>%
  mutate(variable = variable %>% str_replace_all("`", "") %>% str_squish()) %>%
  left_join(
    enframe(clusters_trig_modulators_resultat$RF$res_classification$finalDecision,
      name = "variable",
      value = "Boruta"
    ) %>%
      mutate(variable = str_squish(variable)),
    by = "variable"
  )

selected_vars_clusters_trig_modulators <- names(clusters_trig_modulators_resultat$RF$res_classification$finalDecision)[clusters_trig_modulators_resultat$RF$res_classification$finalDecision == "Confirmed"]

selected_vars_clusters_trig_modulators

Boruta_data_clusters_trig_modulators <- prepare_boruta_plot_data(boruta_res = clusters_trig_modulators_resultat$RF$res_classification)
Boruta_data_clusters_trig_modulators$importance <- Boruta_data_clusters_trig_modulators$importance %>%
  group_by(Feature) %>%
  mutate(median = median(Importance[is.finite(Importance)], na.rm = TRUE)) %>%
  ungroup()

p_clusters_trig_modulators_Boruta <- ggplot(Boruta_data_clusters_trig_modulators$importance, aes(x = reorder(Feature, median), y = Importance, color = Color)) +
  geom_boxplot() +
  labs(
    title = "Relevant features according to 'Boruta' analysis",
    fill = "Decision", color = "Decision",
    x = "Variables", y = "Importance"
  ) +
  theme_plot() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.position = c(.1, .8),
    legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
    legend.key       = element_rect(fill = ggplot2::alpha("white", 0.5)),
    axis.text.y = element_text(size = 8),
    plot.margin = ggplot2::margin(5.5, 25, 5.5, 5.5)
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

clusters_trig_modulators_predictions_rf <- lapply(c("all", "reduced"), function(feature_set) {
  actual_data <- trig_train_final[, !names(trig_train_final) %in% c("Cluster")]
  if (feature_set == "all") {
    mtry <- clusters_trig_modulators_resultat$RF$rf_best_params$mtry
    ntree <- clusters_trig_modulators_resultat$RF$rf_best_params$ntree
  } else {
    actual_data <- actual_data[
      names(clusters_trig_modulators_resultat$RF$res_classification$finalDecision)[clusters_trig_modulators_resultat$RF$res_classification$finalDecision == "Confirmed"]
    ]
    set.seed(42)

    # Quick tune RF for CLASSIFICATION
    mtry_values <- unique(c(1, floor(sqrt(ncol(actual_data))), ceiling(ncol(actual_data) / 3), floor(ncol(actual_data) / 2)))
    ntree_values <- c(100, 200, 500, 1000)

    grid <- expand.grid(mtry = mtry_values, ntree = ntree_values)
    grid$error <- NA
    errors_matrix <- matrix(NA, nrow = nrow(grid), ncol = 10)
    for (run in 1:10) {
      set.seed(42 + run)
      for (i in 1:nrow(grid)) {
        model <- randomForest(
          x = actual_data,
          y = as.factor(clusters_training),
          mtry = grid$mtry[i],
          ntree = grid$ntree[i]
        )
        errors_matrix[i, run] <- mean(model$err.rate[, 1])
      }
    }
    grid$median_error <- apply(errors_matrix, 1, median)
    best <- grid[which.min(grid$median_error), ]

    ntree <- best[["ntree"]]
    mtry <- best[["mtry"]]
  }
  set.seed(42)

  rf_model <- randomForest(
    x = actual_data,
    y = as.factor(clusters_training),
    mtry = mtry,
    ntree = ntree
  )

  rf_repeated <- lapply(1:100, function(s) {
    set.seed(41 + s)
    rows_to_sample <- sample(nrow(trig_valid_final[, -1]), replace = TRUE)
    sampled_df <- trig_valid_final[rows_to_sample, -1]
    sampled_cls <- clusters_validation[rows_to_sample]

    y_new <- predict(rf_model, newdata = sampled_df)
    y_obs <- sampled_cls

    cm_raw <- caret::confusionMatrix(
      factor(y_new, levels = unique(clusters_training)),
      factor(y_obs, levels = unique(clusters_training))
    )$byClass
    cm_rf <- if (length(unique(clusters_training)) > 2) colMeans(cm_raw, na.rm = TRUE) else cm_raw

    df_pred <- cbind.data.frame(y_new, y_obs)

    return(list(cm_rf = cm_rf, df_pred = df_pred))
  })

  cm_rf_rep <- do.call(rbind, lapply(rf_repeated, function(x) x$cm_rf))
  df_pred_rep <- do.call(cbind, lapply(rf_repeated, function(x) x$df_pred))
  return(list(cm_rf_rep = cm_rf_rep, df_pred_rep = df_pred_rep))
})

names(clusters_trig_modulators_predictions_rf) <- c("all", "reduced")

# ============================================================================ #
# 11. COMPILE AND SUMMARIZE RESULTS
# ============================================================================ #

# --- 11a. Significant predictors from penalized multinomial regression ---

coef_table_clusters_trig_modulators <- clusters_trig_modulators_resultat$regression$res_classification$coef_table

sig_predictors_clusters_trig_modulators <- coef_table_clusters_trig_modulators %>%
  dplyr::filter(lasso_selected | elastic_selected | ridge_selected) %>%
  dplyr::arrange(variable, class) %>%
  dplyr::select(
    variable, class,
    glm_coef, glm_selected,
    ridge_coef, ridge_selected,
    lasso_coef, lasso_selected,
    elastic_coef, elastic_selected,
    dplyr::any_of("Boruta")
  )

cat("\n=== Significant predictors (selected by any penalized method) ===\n")
print(sig_predictors_clusters_trig_modulators)

# --- shared helper: bootstrap summary from a rep x metric matrix ---
# cm_mat: 100 x n_metrics matrix; each row is one rep's macro-averaged metrics

make_bootstrap_summary <- function(cm_mat) {
  list(
    mean  = colMeans(cm_mat, na.rm = TRUE),
    ci_lo = apply(cm_mat, 2, quantile, 0.025, na.rm = TRUE),
    ci_hi = apply(cm_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

key_metrics <- c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy")

build_classification_table <- function(predictions_list, label_col) {
  do.call(rbind, lapply(names(predictions_list), function(nm) {
    bs <- make_bootstrap_summary(predictions_list[[nm]]$cm_reg_rep)
    row <- data.frame(setNames(list(nm), label_col), stringsAsFactors = FALSE)
    for (m in key_metrics) {
      safe <- make.names(m)
      row[[paste0(safe, "_mean")]] <- round(bs$mean[m], 3)
      row[[paste0(safe, "_CI_lo")]] <- round(bs$ci_lo[m], 3)
      row[[paste0(safe, "_CI_hi")]] <- round(bs$ci_hi[m], 3)
    }
    row
  }))
}

# --- 11b. Classification performance: penalized regression (bootstrap mean [95% CI]) ---

classification_table_clusters_trig_modulators_reg <- do.call(rbind, lapply(names(clusters_trig_modulators_predictions_reg), function(fs) {
  tbl <- build_classification_table(clusters_trig_modulators_predictions_reg[[fs]], "Model")
  tbl$Model <- dplyr::recode(tbl$Model,
    "multinom_fit" = "Multinomial logistic",
    "elastic"      = "Elastic net",
    "lasso"        = "LASSO",
    "ridge"        = "Ridge"
  )
  tbl$Feature_set <- dplyr::recode(fs,
    "all"     = "All features",
    "reduced" = "Method-selected features"
  )
  tbl[, c("Feature_set", "Model", setdiff(names(tbl), c("Feature_set", "Model")))]
}))

cat("\n=== Classification performance: penalized regression (bootstrap mean [95% CI], 100 reps) ===\n")
print(classification_table_clusters_trig_modulators_reg)

# --- 11c. Classification performance: Random Forest (bootstrap mean [95% CI]) ---

classification_table_clusters_trig_modulators_rf <- do.call(rbind, lapply(names(clusters_trig_modulators_predictions_rf), function(fs) {
  bs <- make_bootstrap_summary(clusters_trig_modulators_predictions_rf[[fs]]$cm_rf_rep)
  row <- data.frame(Feature_set = fs, stringsAsFactors = FALSE)
  for (m in key_metrics) {
    nm <- make.names(m)
    row[[paste0(nm, "_mean")]] <- round(bs$mean[m], 3)
    row[[paste0(nm, "_CI_lo")]] <- round(bs$ci_lo[m], 3)
    row[[paste0(nm, "_CI_hi")]] <- round(bs$ci_hi[m], 3)
  }
  row
}))

classification_table_clusters_trig_modulators_rf$Feature_set <- dplyr::recode(classification_table_clusters_trig_modulators_rf$Feature_set,
  "all"     = "All features",
  "reduced" = "Boruta-selected features"
)

cat("\n=== Classification performance: Random Forest (bootstrap mean [95% CI], 100 reps) ===\n")
print(classification_table_clusters_trig_modulators_rf)

# ============================================================================ #
# 12. SAVE RESULTS
# ============================================================================ #

write.csv(sig_predictors_clusters_trig_modulators, "clusters_trig_modulators_sig_predictors.csv", row.names = FALSE)
write.csv(classification_table_clusters_trig_modulators_reg, "clusters_trig_modulators_classification_reg.csv", row.names = FALSE)
write.csv(classification_table_clusters_trig_modulators_rf, "clusters_trig_modulators_classification_rf.csv", row.names = FALSE)
write.csv(coef_table_clusters_trig_modulators, "clusters_trig_modulators_coef_table_full.csv", row.names = FALSE)
write.csv(Boruta_data_clusters_trig_modulators$importance[!duplicated(Boruta_data_clusters_trig_modulators$importance$Feature), ], "Boruta_importance_trigeminal_sensitivity_clusters_modulators.csv")

cat("Saved: clusters_trig_modulators_sig_predictors.csv\n")
cat("Saved: clusters_trig_modulators_classification_reg.csv\n")
cat("Saved: clusters_trig_modulators_classification_rf.csv\n")
cat("Saved: clusters_trig_modulators_coef_table_full.csv\n")
cat("Saved: Boruta_importance_trigeminal_sensitivity_clusters_modulators.csv\n")

ggsave("p_clusters_trig_modulators_Boruta.svg", p_clusters_trig_modulators_Boruta, width = 16, height = 10, dpi = 300, limitsize = FALSE)
ggsave("p_clusters_trig_modulators_Boruta.png", p_clusters_trig_modulators_Boruta, width = 16, height = 10, dpi = 300, limitsize = FALSE)

# ============================================================================ #
# END OF SCRIPT
# ============================================================================ #
