################################################################################
# Trigeminal Cluster Classification Analysis
#
# Author: Jorn Lotsch
# Date: 2026-04-29
#
# Description:
#   Classification of trigeminal cluster membership using clustered training and
#   validation datasets. Combines Boruta feature selection with penalized
#   multinomial regression (ridge, lasso, elastic net) and Random Forest
#   classification. Models are trained on the training set and evaluated on the
#   validation set.
#
# Input Files:
#   - trigeminal_clustered_data.csv
#   - globals.R
#   - utils.R
#
# Output Files:
#   - clusters_trig_sig_predictors.csv
#   - clusters_trig_classification_reg.csv
#   - clusters_trig_classification_rf.csv
#   - clusters_trig_coef_table_full.csv
#   - p_clusters_Boruta.svg
#   - p_clusters_Boruta.png
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
  # Map model.matrix column names back to original data column names
  # (handles both direct matches and make.names-encoded names)
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

# Load all data and cluster assignments
all_trig_clustered_data_and_clusters <- read.csv("trigeminal_clustered_data.csv", row.names = 1, check.names = FALSE)
rownames(all_trig_clustered_data_and_clusters) <- all_trig_clustered_data_and_clusters$ID


trig_clusters_training_data <- all_trig_clustered_data_and_clusters[
  all_trig_clustered_data_and_clusters$Set == "Training",
  !names(all_trig_clustered_data_and_clusters) %in% c("Set", "ID")
]
trig_clusters_validation_data <- all_trig_clustered_data_and_clusters[
  all_trig_clustered_data_and_clusters$Set == "Validation",
  !names(all_trig_clustered_data_and_clusters) %in% c("Set", "ID")
]
trig_clusters_training_data$Cluster <- as.factor(trig_clusters_training_data$Cluster)
trig_clusters_validation_data$Cluster <- as.factor(trig_clusters_validation_data$Cluster)

clusters_training <- as.factor(trig_clusters_training_data$Cluster)
clusters_validation <- as.factor(trig_clusters_validation_data$Cluster)


# ============================================================================ #
# 6. BORUTA FEATURE SELECTION AND PENALIZED REGRESSION
# ============================================================================ #

clusters_trig_resultat <- pbmcapply::pbmclapply(c("regression", "RF"), function(method) {
  rf_best_params <- NA
  actual_data <- trig_clusters_training_data[, !names(trig_clusters_training_data) %in% c("Cluster")]
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

names(clusters_trig_resultat) <- c("regression", "RF")

# ============================================================================ #
# 7. INITIAL RESULT INSPECTION
# ============================================================================ #

par(mfrow = c(2, 1))
plot(clusters_trig_resultat$RF$res_classification, las = 2)
par(mfrow = c(1, 1))

print(clusters_trig_resultat$regression$res_classification)

# ============================================================================ #
# 8. CLASSIFCATION MODEL PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

clusters_trig_predictions_reg <- lapply(c("all", "reduced"), function(feature_set) {
  ct         <- clusters_trig_resultat$regression$res_classification$coef_table
  train_full <- trig_clusters_training_data[, !names(trig_clusters_training_data) %in% c("Cluster")]
  val_full   <- trig_clusters_validation_data[, !names(trig_clusters_validation_data) %in% c("Cluster")]

  fits <- lapply(c("multinom_fit", "elastic", "lasso", "ridge"), function(fit_type) {
    print(paste(feature_set, fit_type))

    if (feature_set == "all") {
      val_data <- val_full
      model    <- clusters_trig_resultat$regression$res_classification[[fit_type]]
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
names(clusters_trig_predictions_reg) <- c("all", "reduced")


# ============================================================================ #
# 9. BORUTA IMPORTANCE VISUALIZATION
# ============================================================================ #

clusters_trig_resultat$regression$res_classification$coef_table <-
  clusters_trig_resultat$regression$res_classification$coef_table %>%
  mutate(variable = variable %>% str_replace_all("`", "") %>% str_squish()) %>%
  left_join(
    enframe(clusters_trig_resultat$RF$res_classification$finalDecision,
      name = "variable",
      value = "Boruta"
    ) %>%
      mutate(variable = str_squish(variable)),
    by = "variable"
  )


selected_vars <- names(clusters_trig_resultat$RF$res_classification$finalDecision)[clusters_trig_resultat$RF$res_classification$finalDecision == "Confirmed"]
selected_vars


Boruta_data <- prepare_boruta_plot_data(boruta_res = clusters_trig_resultat$RF$res_classification)
Boruta_data$importance <- Boruta_data$importance %>%
  group_by(Feature) %>%
  mutate(median = median(Importance[is.finite(Importance)], na.rm = TRUE)) %>%
  ungroup()


p_clusters_Boruta <- ggplot(Boruta_data$importance, aes(x = reorder(Feature, median), y = Importance, color = Color)) +
  geom_boxplot() +
  labs(
    title = "Relevant features according to 'Boruta' analysis",
    fill = "Decision", color = "Decision",
    x = "Variables", y = "Importance"
  ) +
  theme_plot() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.margin = ggplot2::margin(5.5, 25, 5.5, 5.5)
  ) +
  scale_color_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  ) +
  scale_fill_manual(
    values = c("Chosen" = "chartreuse4", "Tentative" = "gold", "Rejected" = "salmon", "Shadow" = "grey50")
  )


# SECOND VERSION OF COMBINED PLOTS FOR TRAINING DATA

cat("\n=== Creating Combined Training Plot ===\n")

# Combine effect size plot (left) + violin plot (right)
p_left_clusters_train <- (p_effect_training + theme(plot.margin = ggplot2::margin(5, 5, 5, 5)) + labs(title = "Kruskal-Wallis η²")) |
  (p_clusters_Boruta + theme(plot.margin = ggplot2::margin(5, 5, 5, 5)) + coord_flip(clip = "off") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) + labs(title = "Boruta (RF) importance"))
p_right_clusters_train <- p_trigeminal_clustered_training_data + theme(plot.margin = ggplot2::margin(5, 5, 5, 5))

p_combined_clusters_train_sup <- (p_left_clusters_train | p_right_clusters_train) +
  plot_layout(widths = c(1, 1, 5)) +
  plot_annotation(
    title = "Cluster differences in trigeminal measures: Training data",
    tag_levels = LETTERS,
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  ) &
  theme(plot.margin = ggplot2::margin(5, 5, 5, 5))

print(p_combined_clusters_train_sup)

# Save combined training plot
ggsave("p_combined_clusters_train_sup.svg", p_combined_clusters_train_sup,
  width = 22, height = 15, dpi = 300, limitsize = FALSE
)
ggsave("p_combined_clusters_train_sup.png", p_combined_clusters_train_sup,
  width = 22, height = 15, dpi = 300, limitsize = FALSE
)

cat("Combined training plot saved\n")


# ============================================================================ #
# 10. RANDOM FOREST PREDICTIONS ON VALIDATION DATA
# ============================================================================ #

clusters_trig_predictions_rf <- lapply(c("all", "reduced"), function(feature_set) {
  actual_data <- trig_clusters_training_data[, !names(trig_clusters_training_data) %in% c("Cluster")]
  if (feature_set == "all") {
    mtry <- clusters_trig_resultat$RF$rf_best_params$mtry
    ntree <- clusters_trig_resultat$RF$rf_best_params$ntree
  } else {
    actual_data <- actual_data[
      names(clusters_trig_resultat$RF$res_classification$finalDecision)[clusters_trig_resultat$RF$res_classification$finalDecision == "Confirmed"]
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
    rows_to_sample <- sample(nrow(trig_clusters_validation_data[, -1]), replace = TRUE)
    sampled_df <- trig_clusters_validation_data[rows_to_sample, -1]
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


names(clusters_trig_predictions_rf) <- c("all", "reduced")


# ============================================================================ #
# 11. COMPILE AND SUMMARIZE RESULTS
# ============================================================================ #

# --- 11a. Significant predictors from penalized multinomial regression ---

coef_table_clusters <- clusters_trig_resultat$regression$res_classification$coef_table

sig_predictors_clusters <- coef_table_clusters %>%
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
print(sig_predictors_clusters)


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

classification_table_reg <- do.call(rbind, lapply(names(clusters_trig_predictions_reg), function(fs) {
  tbl <- build_classification_table(clusters_trig_predictions_reg[[fs]], "Model")
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
print(classification_table_reg)


# --- 11c. Classification performance: Random Forest (bootstrap mean [95% CI]) ---

classification_table_rf <- do.call(rbind, lapply(names(clusters_trig_predictions_rf), function(fs) {
  bs <- make_bootstrap_summary(clusters_trig_predictions_rf[[fs]]$cm_rf_rep)
  row <- data.frame(Feature_set = fs, stringsAsFactors = FALSE)
  for (m in key_metrics) {
    nm <- make.names(m)
    row[[paste0(nm, "_mean")]] <- round(bs$mean[m], 3)
    row[[paste0(nm, "_CI_lo")]] <- round(bs$ci_lo[m], 3)
    row[[paste0(nm, "_CI_hi")]] <- round(bs$ci_hi[m], 3)
  }
  row
}))

classification_table_rf$Feature_set <- dplyr::recode(classification_table_rf$Feature_set,
  "all"     = "All features",
  "reduced" = "Boruta-selected features"
)

cat("\n=== Classification performance: Random Forest (bootstrap mean [95% CI], 100 reps) ===\n")
print(classification_table_rf)


# ============================================================================ #
# 12. SAVE RESULTS
# ============================================================================ #

write.csv(sig_predictors_clusters, "clusters_trig_sig_predictors.csv", row.names = FALSE)
write.csv(classification_table_reg, "clusters_trig_classification_reg.csv", row.names = FALSE)
write.csv(classification_table_rf, "clusters_trig_classification_rf.csv", row.names = FALSE)
write.csv(coef_table_clusters, "clusters_trig_coef_table_full.csv", row.names = FALSE)
write.csv(Boruta_data$importance[!duplicated(Boruta_data$importance$Feature), ], "Boruta_importance_trigeminal_sensitivity_clusters.csv")

cat("Saved: clusters_trig_sig_predictors.csv\n")
cat("Saved: clusters_trig_classification_reg.csv\n")
cat("Saved: clusters_trig_classification_rf.csv\n")
cat("Saved: clusters_trig_coef_table_full.csv\n")
cat("Saved: Boruta_importance_trigeminal_sensitivity_clusters.csv\n")

ggsave("p_clusters_Boruta.svg", p_clusters_Boruta +
         theme(
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
           legend.position = c(.1, .8),
           legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
           legend.key       = element_rect(fill = ggplot2::alpha("white", 0.5))
         ), width = 12, height = 12, dpi = 300, limitsize = FALSE)
ggsave("p_clusters_Boruta.png", p_clusters_Boruta +
         theme(
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
           legend.position = c(.1, .8),
           legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
           legend.key       = element_rect(fill = ggplot2::alpha("white", 0.5))
         ), width = 12, height = 12, dpi = 300, limitsize = FALSE)



# ============================================================================ #
# 12. TEST WHETHER CLUSTERS CAN BE IDENTIFED BY ONE VARIABLE ONLY
# ============================================================================ #

clusters_trig_bf <- lapply(names(trig_clusters_training_data)[2:ncol(trig_clusters_training_data)], function(feature) {
  print(feature)
  actual_data <- trig_clusters_training_data[, !names(trig_clusters_training_data) %in% c("Cluster")]
  actual_data <- as.data.frame(trig_clusters_training_data[, feature])
  names(actual_data) <- feature
  set.seed(42)

  # # Quick tune RF for CLASSIFICATION
  mtry_values <- sqrt(ncol(actual_data))
  ntree_values <- c(200, 500, 1000)

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

  set.seed(42)

  rf_model <- randomForest(
    x = actual_data,
    y = as.factor(clusters_training),
    mtry = mtry,
    ntree = ntree
  )

  rf_repeated <- lapply(1:100, function(s) {
    set.seed(41 + s)
    rows_to_sample <- sample(nrow(trig_clusters_validation_data[, -1]), replace = TRUE)
    sampled_df <- trig_clusters_validation_data[rows_to_sample, -1]
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
})#, mc.cores = 6)


names(clusters_trig_bf) <- names(trig_clusters_training_data)[2:ncol(trig_clusters_training_data)]
res_one_feature_classification <- lapply(lapply(clusters_trig_bf, "[[", 1), function(df) {
  quantile(df[, "Balanced Accuracy"], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
})

tab_one_feature_classification <- do.call(rbind, lapply(names(res_one_feature_classification), function(n) {
  cbind(varname = n, as.data.frame(t(res_one_feature_classification[[n]])))
}))

rownames(tab_one_feature_classification) <- NULL
tab_one_feature_classification

tab_filteredone_feature_classification <- tab_one_feature_classification[as.numeric(tab_one_feature_classification$`2.5%`) > 0.5, ]
tab_filteredone_feature_classification


# ============================================================================ #
# 13. INTERPREST CLUSTERS BASED ON FACTOR ANALYIS
# ============================================================================ #

# Load coordinates for later cluster interpretation
FA_coordinates <- read.csv("FA_PA_coordinates.csv", row.names = 1)
FA_coordinates$ID <- rownames(FA_coordinates)

FA_coordinates_training <- merge(
  FA_coordinates,
  trig_clusters_training_data[, c("ID", "Cluster")],
  by = "ID"
)

FA_coordinates_training$ID <- NULL

FA_coordinates_training_means <- FA_coordinates_training %>%
  group_by(Cluster) %>%
  summarise(
    mean_col1 = mean(PA1, na.rm = TRUE),
    mean_col2 = mean(PA2, na.rm = TRUE),
    .groups = "drop"
  )


FA_coordinates_validation <- merge(
  FA_coordinates,
  trig_clusters_validation_data[, c("ID", "Cluster")],
  by = "ID"
)

FA_coordinates_validation$ID <- NULL

FA_coordinates_validation_means <- FA_coordinates_validation %>%
  group_by(Cluster) %>%
  summarise(
    mean_col1 = mean(PA1, na.rm = TRUE),
    mean_col2 = mean(PA2, na.rm = TRUE),
    .groups = "drop"
  )


print(FA_coordinates_training_means)
print(FA_coordinates_validation_means)



# ============================================================================ #
# END OF SCRIPT
# ============================================================================ #
