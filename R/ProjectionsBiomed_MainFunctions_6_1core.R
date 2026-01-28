############### Functions ProjectionsBiomed new #########################################################
# Extended colorblind palette
cb_palette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

safe_extract_cluster_data <- function(cluster_list) {
  if (is.null(cluster_list)) stop("Error: `cluster_list` is NULL.")
  extracted_data <- tryCatch(
    extract_cluster_data(cluster_list),
    error = function(e) stop("Error while extracting cluster data:", conditionMessage(e))
  )
  return(extracted_data)
}

# Helper function for row medians
rowMedians <- function(x, na.rm = TRUE) {
  if (!is.data.frame(x) | !is.matrix(x)) {
    x
  }
  apply(x, 1, median, na.rm = na.rm)
}

############### Libraries ###############

# Package loading function, avoiding installation during the function execution
# Constants for required package names
REQUIRED_PACKAGES <- c(
  "pbmcapply", "deldir", "ggplot2", "ggrepel", "ggplotify",
  "grid", "cowplot", "cluster", "FactoMineR", "ABCanalysis",
  "fastICA", "MASS", "Rtsne", "RDRToolbox", "mixOmics",
  "umap", "uwot", "caret", "pracma", "combinat", "NbClust", "parallel",
  "dplyr", "ape", "combinat", "NMF", "RANN", "pls", "Matrix",
  "twosamples", "car", "stringr"
)

############### Main function ###############

# Function to load required R packages
load_required_packages <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  invisible(lapply(packages, library, character.only = TRUE))
}

# Function to perform analysis
perform_analysis <- function(datasets,
                              projection_methods = "PCA",
                              clustering_methods = "none",
                              cluster_number_methods = "orig",
                              selected_cluster_metrics = c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                                                            "Rand_index", "DaviesBouldin_index", "dbcv_index",
                                                            "CalinskiHarabasz_index", "inertia", "adjusted_mutual_information"),
                              distance_metric = "euclidean",
                              highlight_best_clustering = FALSE,
                              method_for_hcpc = "ward",
                              label_points = TRUE,
                              highlight_misclassified = FALSE,
                              points_colored_for = "Target",
                              cells_colored_for = "Cluster",
                              palette_target = NULL,
                              palette_cluster = NULL,
                              seed = 42,
                              jitter_one_dimensional_projections = TRUE,
                              nProc = 12,
                              max_clusters = 5,
                              scaleX = TRUE,
                              remove_colinear_vars = FALSE,
                              correlation_threshold_for_colinearity = 0.9,
                              vif_threshold_for_colinearity = 10) {
  # Validate input parameters
  if (missing(datasets) || length(datasets) == 0) {
    stop("The `datasets` parameter must be a non-empty list of dataset names.")
  }

  if (!is.logical(highlight_best_clustering)) {
    stop("The argument `highlight_best_clustering` must be a logical (TRUE or FALSE).")
  }

  if (!is.logical(highlight_misclassified)) {
    stop("The argument `highlight_misclassified` must be a logical (TRUE or FALSE).")
  }

  if (!is.logical(jitter_one_dimensional_projections)) {
    stop("The argument `jitter_one_dimensional_projections` must be a logical (TRUE or FALSE).")
  }

  if (!is.logical(scaleX)) {
    stop("The argument `scaleX` must be a logical (TRUE or FALSE).")
  }

  # Call the validation function
  validated_params <- validate_and_filter_parameters(
    projection_methods = projection_methods,
    clustering_methods = clustering_methods,
    method_for_hcpc = method_for_hcpc,
    distance_metric = distance_metric,
    cluster_number_methods = cluster_number_methods,
    selected_cluster_metrics = selected_cluster_metrics
  )

  # Extract the validated and filtered parameters
  projection_methods <- validated_params$projection_methods
  clustering_methods <- validated_params$clustering_methods
  method_for_hcpc <- validated_params$method_for_hcpc
  distance_metric <- validated_params$distance_metric
  cluster_number_methods <- validated_params$cluster_number_methods
  selected_cluster_metrics <- validated_params$selected_cluster_metrics

  # Ensure random seed is set for reproducibility
  set.seed(seed)

  # Initialize output variables
  projections_plots <- list() # Store projection-related plots
  all_projection_results <- list() # Store projection results
  all_clustering_results <- list() # Store clustering results
  cluster_quality_results <- list() # Store cluster quality scores

  # Loop through all datasets
  for (dataset_name in datasets) {
    message("Processing dataset: ", dataset_name)

    # Fetch dataset
    dataset <- retrieve_dataset(dataset_name)
    if (is.null(dataset)) {
      warning("Dataset '", dataset_name, "' could not be retrieved. Skipping...")
      next
    }

    # Handle 'Target' with only one level and adjust projection methods
    local_projection_methods <- projection_methods
    if (length(unique(dataset$Target)) == 1 && "PLSDA" %in% local_projection_methods) {
      message("'Target' has only one level. Removing PLSDA from projection methods list.")
      local_projection_methods <- setdiff(local_projection_methods, "PLSDA")
      if (length(local_projection_methods) < 1) {
        warning("No available projection methods for dataset: ", dataset_name, ". Skipping...")
        next
      }
    }

    # Generate projections
    message("Generating projections...")
    on.exit({
      cat("Cleaning up processes...\n")
      parallel::mccollect(wait = FALSE)
    })
    tryCatch({
      projection_results <- generate_projections(
        data = dataset,
        projection_methods = local_projection_methods,
        seed = seed,
        jitter_one_dimensional_projections = jitter_one_dimensional_projections,
        nProc = min(nProc, length(local_projection_methods)),
        scaleX = scaleX,
        remove_colinear_vars = remove_colinear_vars,
        correlation_threshold_for_colinearity = correlation_threshold_for_colinearity,
        vif_threshold_for_colinearity = vif_threshold_for_colinearity
      )
      all_projection_results[[dataset_name]] <- projection_results
    }, error = function(e) {
      warning("Projection generation failed for dataset '", dataset_name, "' with error: ", e$message)
      next
    })
    # parallel::mccollect( wait = TRUE )
    # Sys.sleep( 3 )
    #
    # Apply clustering
    message("Applying clustering...")
    tryCatch({
      clustering_results <- apply_clustering(
        projection_results = projection_results,
        clustering_methods = clustering_methods,
        cluster_number_methods = cluster_number_methods,
        method_for_hcpc = method_for_hcpc,
        distance_metric = distance_metric,
        seed = seed,
        nProc = min(nProc, length(local_projection_methods)),
        max_clusters = max_clusters
      )
      all_clustering_results[[dataset_name]] <- clustering_results
    }, error = function(e) {
      warning("Clustering failed for dataset '", dataset_name, "' with error: ", e$message)
      next
    })
    # parallel::mccollect( wait = TRUE )
    # Calculate cluster quality indices
    message("Calculating cluster stability or quality indices...")
    tryCatch({
      cluster_indices <- calculate_cluster_scores_df(
        clustering_results = clustering_results,
        projection_methods = local_projection_methods,
        clustering_methods = clustering_methods,
        cluster_number_methods = cluster_number_methods,
        distance_metric = distance_metric,
        nProc = min(nProc, length(local_projection_methods))
      )
      cluster_quality_results[[dataset_name]] <- process_cluster_scores(
        cluster_scores_df = cluster_indices,
        selected_cluster_metrics = selected_cluster_metrics
      )
    }, error = function(e) {
      warning("Cluster quality calculation failed for dataset '", dataset_name, "' with error: ", e$message)
      next
    })
    # parallel::mccollect( wait = TRUE )
    # Create plots
    message("Creating plots...")
    tryCatch({
      plots <- create_dataset_plots(
        clustering_results = clustering_results,
        projection_methods = local_projection_methods,
        clustering_methods = clustering_methods,
        cluster_number_methods = cluster_number_methods,
        dataset_name = dataset_name,
        label_points = label_points,
        highlight_misclassified = highlight_misclassified,
        points_colored_for = points_colored_for,
        cells_colored_for = cells_colored_for,
        palette_target = palette_target,
        palette_cluster = palette_cluster,
        nProc = min(nProc, length(local_projection_methods))
      )

      # Highlight best clustering if required
      if (highlight_best_clustering) {
        plots <- do_highlight_best_clustering(
          cluster_quality_results = cluster_quality_results,
          dataset_name = dataset_name,
          plots = plots
        )
      }

      projections_plots[[dataset_name]] <- unlist(plots, recursive = FALSE)

    }, error = function(e) {
      warning("Plot creation failed for dataset '", dataset_name, "' with error: ", e$message)
      next
    })
  }
  # End of dataset loop

  # Final return
  return(list(
    projections_plots = projections_plots,
    projection_results = all_projection_results,
    clustering_results = all_clustering_results,
    cluster_quality_results = cluster_quality_results
  ))
}


do_highlight_best_clustering <- function(cluster_quality_results, dataset_name, plots) {
  # Retrieve the best combination rank
  best_combination_rank <- max(cluster_quality_results[[dataset_name]][["combined_rank"]])

  # Determine the second-best combination rank (if applicable)
  if (length(unique(cluster_quality_results[[dataset_name]][["combined_rank"]])) > 1) {
    second_best_combination_rank <- sort(unique(cluster_quality_results[[dataset_name]][["combined_rank"]]), decreasing = TRUE)[2]
  } else {
    second_best_combination_rank <- -1
  }

  # Extract best combination details
  best_combination <- cbind.data.frame(
    FirstSecond = 1,
    cluster_quality_results[[dataset_name]][cluster_quality_results[[dataset_name]]$combined_rank == best_combination_rank,
                                            c("projection_method", "clustering_method", "cluster_number_method")]
  )

  best_second_best_combination <- best_combination

  # If a second-best combination exists, add it to the result
  if (second_best_combination_rank > 0) {
    second_best_combination <- cbind.data.frame(
      FirstSecond = 2,
      cluster_quality_results[[dataset_name]][cluster_quality_results[[dataset_name]]$combined_rank == second_best_combination_rank,
                                              c("projection_method", "clustering_method", "cluster_number_method")]
    )

    best_second_best_combination <- rbind.data.frame(
      best_combination,
      second_best_combination
    )

    # Remove duplicates
    best_second_best_combination <- best_second_best_combination[
      !duplicated(best_second_best_combination[, c("projection_method", "clustering_method", "cluster_number_method")]),
    ]
  }

  # Update plots for the best and second-best combinations
  for (i in seq_len(nrow(best_second_best_combination))) {
    if (best_second_best_combination$FirstSecond[i] == 1) {
      plots[[best_second_best_combination[i, "projection_method"]]][[best_second_best_combination[i, "clustering_method"]]][[best_second_best_combination[i, "cluster_number_method"]]] <-
        suppressMessages(
          suppressMessages(
            plots[[best_second_best_combination[i, "projection_method"]]][[best_second_best_combination[i, "clustering_method"]]][[best_second_best_combination[i, "cluster_number_method"]]] +
              scale_color_hue() +
              scale_fill_hue()
          )
        )
    }

    if (best_second_best_combination$FirstSecond[i] == 2) {
      plots[[best_second_best_combination[i, "projection_method"]]][[best_second_best_combination[i, "clustering_method"]]][[best_second_best_combination[i, "cluster_number_method"]]] <-
        suppressMessages(
          suppressMessages(
            plots[[best_second_best_combination[i, "projection_method"]]][[best_second_best_combination[i, "clustering_method"]]][[best_second_best_combination[i, "cluster_number_method"]]] +
              scale_color_viridis_d() +
              scale_fill_viridis_d()
          )
        )
    }
  }

  return(plots)
}


############### Functions for data preparations ###############

# Function to prepare datasets
prepare_dataset <- function(input_X, Target = NULL, Label = NULL) {
  input_X <- if (is.matrix(input_X)) as.data.frame(input_X) else input_X
  if (!is.data.frame(input_X) && !is.matrix(input_X)) stop("Input must be a data frame or matrix.")

  # Use existing "Target" column if Target is NULL and "Target" exists in input_X
  output_Y <- if (is.null(Target)) {
    if (!"Target" %in% colnames(input_X)) {
      message("'Target' is missing, creating 'Target' = 1.")
      rep(1, nrow(input_X))
    } else {
      input_X$Target
    }
  } else {
    Target
  }

  # Use existing "Label" column if Labels is NULL and "Target" exists in input_X
  output_L <- if (is.null(Label)) {
    if (!"Label" %in% colnames(input_X)) {
      message("Taking row names as case labels.")
      rownames(input_X)
    } else {
      message("Taking 'Label' column as case labels.")
      input_X$Label
    }
  } else {
    if (length(Label) != nrow(input_X)) {
      message("Length of 'Label' does not match numbr of rows in input matrix\nTaking row names as case labels.")
      rownames(input_X)
    } else {
      Label
    }
  }


  data_frame <- as.data.frame(input_X)
  if (ncol(data_frame) < 3) stop("Matrix needs at least three columns including 'Target'.")

  # Ensure "Target" is numeric after checking availability
  data_frame$Target <- output_Y
  if (!is.numeric(data_frame$Target)) {
    data_frame$Target <- as.numeric(factor(data_frame$Target))
  }

  # Add "Label" is numeric after checking availability
  data_frame$Label <- output_L

  return(data_frame)
}


# Function to retrieve datasets
retrieve_dataset <- function(dataset_name) {
  if (!exists(dataset_name, envir = .GlobalEnv)) {
    warning("Dataset ", dataset_name, " not found. Skipping.")
    return(NULL)
  }

  data_set <- get(dataset_name, envir = .GlobalEnv)
  if (!is.data.frame(data_set)) {
    warning("Dataset ", dataset_name, " is not a data frame. Skipping.")
    return(NULL)
  }

  if (!"Target" %in% names(data_set)) {
    warning("No 'Target' column found in dataset ", dataset_name, ". Skipping.")
    return(NULL)
  }

  return(data_set)
}

# Combine results for plotting
create_combined_result <- function(projected_data, clusters, projection_result) {
  combined_result <- cbind.data.frame(
    projected_data,
    Cluster = clusters,
    Target = projection_result$UniqueData$Target,
    Label = projection_result$UniqueData$Label
  )
  combined_result$Misclassified <- ifelse(combined_result$Cluster == combined_result$Target, 0, 1)
  return(combined_result)
}

# Function for validation of the methods inputs
validate_and_filter_parameters <- function(projection_methods,
                                            clustering_methods,
                                            method_for_hcpc,
                                            distance_metric,
                                            cluster_number_methods,
                                            selected_cluster_metrics) {
  # Define valid methods
  valid_projection_methods <- c("none", "PCA", "ICA", "MDS", "tSNE", "isomap", "LDA",
                                 "PLSDA", "PLSLDA", "IPCA", "Umap", "NMF", "LLE", "autoencoder")
  valid_clustering_methods <- c("none", "kmeans", "kmedoids",
                                 "diana", "hcpc", "ward.D2", "single", "average", "median", "complete", "centroid")
  valid_hcpc_methods <- c("average", "single", "complete", "ward", "weighted", "flexible", "gaverage")
  valid_distances <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  valid_cluster_number_methods <- c("orig", "NbClust", "hcpc")
  valid_cluster_metrics <- c("cluster_accuracy", "Silhouette_index", "Dunn_index", "Rand_index", "DaviesBouldin_index", "dbcv_index",
                              "CalinskiHarabasz_index", "inertia", "adjusted_mutual_information")

  # Helper function to check input lists
  check_lists <- function(list1, list2) {
    vec1 <- unlist(list1)
    vec2 <- unlist(list2)
    not_in_list2 <- setdiff(vec1, vec2)
    if (length(not_in_list2) > 0) {
      return(list(result = TRUE, elements = not_in_list2))
    } else {
      return(list(result = FALSE, elements = NULL))
    }
  }

  # Check if all methods are required
  if ("all" %in% projection_methods) projection_methods <- valid_projection_methods
  if ("all" %in% clustering_methods) clustering_methods <- valid_clustering_methods
  if ("all" %in% cluster_number_methods) cluster_number_methods <- valid_cluster_number_methods
  if ("all" %in% projection_methods) projection_methods <- valid_projection_methods
  if ("all" %in% selected_cluster_metrics) selected_cluster_metrics <- valid_cluster_metrics

  # Validate input parameters BEFORE reassignment
  if (check_lists(projection_methods, valid_projection_methods)$result) {
    stop("Invalid `projection_methods`. Please choose from: ", paste(valid_projection_methods, collapse = ", "))
  }
  if (check_lists(clustering_methods, valid_clustering_methods)$result) {
    stop("Invalid `clustering_methods`. Please choose from: ", paste(valid_clustering_methods, collapse = ", "))
  }
  if (check_lists(method_for_hcpc, valid_hcpc_methods)$result) {
    stop("Invalid `method_for_hcpc`. Please choose from: ", paste(valid_hcpc_methods, collapse = ", "))
  }
  if (check_lists(distance_metric, valid_distances)$result) {
    stop("Invalid `distance_metric`. Please choose from: ", paste(valid_distances, collapse = ", "))
  }
  if (check_lists(cluster_number_methods, valid_cluster_number_methods)$result) {
    stop("Invalid `cluster_number_methods`. Please choose from: ", paste(valid_cluster_number_methods, collapse = ", "))
  }
  if (check_lists(selected_cluster_metrics, valid_cluster_metrics)$result) {
    stop("Invalid `selected_cluster_metrics`. Please choose from: ", paste(valid_cluster_metrics, collapse = ", "))
  }

  # AFTER validation, safely filter parameters
  projection_methods <- intersect(projection_methods, valid_projection_methods)
  clustering_methods <- intersect(clustering_methods, valid_clustering_methods)
  method_for_hcpc <- intersect(method_for_hcpc, valid_hcpc_methods)
  distance_metric <- intersect(distance_metric, valid_distances)
  cluster_number_methods <- intersect(cluster_number_methods, valid_cluster_number_methods)
  selected_cluster_metrics <- intersect(selected_cluster_metrics, valid_cluster_metrics)

  # Check if any method remains
  if (length(projection_methods) < 1) {
    stop(paste("Invalid projection_methods. Choose from:", paste(valid_projection_methods, collapse = ", ")))
  }
  if (length(clustering_methods) < 1) {
    stop(paste("Invalid clustering_methods. Choose from:", paste(valid_clustering_methods, collapse = ", ")))
  }
  if (length(method_for_hcpc) < 1) {
    stop(paste("Invalid method_for_hcpc method. Choose from:", paste(valid_hcpc_methods, collapse = ", ")))
  }
  if (length(distance_metric) < 1) {
    stop(paste("Invalid distance_metric. Choose from:", paste(valid_distances, collapse = ", ")))
  }
  if (length(cluster_number_methods) < 1) {
    stop(paste("Invalid cluster_number_methods. Choose from:", paste(valid_cluster_number_methods, collapse = ", ")))
  }
  if (length(selected_cluster_metrics) < 1) {
    stop(paste("Invalid selected_cluster_metrics. Choose from:", paste(valid_cluster_metrics, collapse = ", ")))
  }

  # Return filtered parameters as a list
  return(list(projection_methods = projection_methods,
                clustering_methods = clustering_methods,
                method_for_hcpc = method_for_hcpc,
                distance_metric = distance_metric,
                cluster_number_methods = cluster_number_methods,
                selected_cluster_metrics = selected_cluster_metrics))
}


############### Functions for data projection ###############

# Function to generate projections
performProjection <- function(X,
                               projection_method = "PCA",
                               scaleX = TRUE,
                               Target = NULL,
                               switchDims = TRUE,
                               labels = NULL,
                               seed = 42,
                               jitter_one_dimensional_projections = FALSE,
                               remove_colinear_vars = FALSE,
                               vif_threshold_for_colinearity = 10,
                               correlation_threshold_for_colinearity = 1,
                               newdata = NULL,
                               projection_object = NULL) {

  # Check required libraries
  required_packages <- c("FactoMineR", "fastICA", "MASS", "Rtsne", "RDRToolbox", "mixOmics",
                          "ABCanalysis", "umap", "NMF", "neuralnet", "pls", "car", "RANN")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it."))
    }
  }

  valid_projection_methods <- c("none", "PCA", "ICA", "MDS", "tSNE", "isomap", "LDA",
                                 "PLSDA", "PLSLDA", "IPCA", "Umap", "NMF", "LLE", "autoencoder")
  if (!(projection_method %in% valid_projection_methods)) {
    stop(paste("Invalid projection method. Choose from:", paste(valid_projection_methods, collapse = ", ")))
  }

  if (!is.data.frame(X) && !is.matrix(X)) stop("X must be a data frame or matrix.")
  if (is.null(X) || nrow(X) == 0) stop("X must not be NULL or empty.")

  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  # TRAINING MODE
  if (is.null(newdata)) {
    if (!is.null(Target) && length(Target) != nrow(X)) stop("Target length must match the number of rows in X.")
    if (!is.null(labels) && length(labels) != nrow(X)) stop("Labels length must match the number of rows in X.")

    dupl <- which(duplicated(X))
    if (!is.null(Target) && projection_method != "PCA" && !is.integer0(dupl)) {
      Target <- Target[-dupl]
    }
    if (!is.null(labels) && projection_method != "PCA" && !is.integer0(dupl)) {
      labels <- labels[-dupl]
    }
    if (projection_method != "PCA") X <- X[!duplicated(X),]

    if (scaleX && !projection_method %in% c("IPCA", "NMF")) {
      X <- scale(X)
    }

    if (is.null(labels)) {
      if (!is.null(rownames(X))) {
        labels <- rownames(X)
      } else {
        labels <- seq_len(nrow(X))
      }
    }

    if (!is.null(Target)) {
      df_unique <- data.frame(X, Target = Target, Label = labels)
    } else {
      df_unique <- data.frame(X, Label = labels)
    }

    X_unique <- df_unique[, setdiff(names(df_unique), c("Target", "Label")), drop = FALSE]
    Target_unique <- if (!is.null(Target)) df_unique$Target else NULL

    # Collinearity removal function (complete)
    remove_collinear_with_fallback <- function(df, target_col, vif_threshold = 10, correlation_threshold = 0.9) {
      if (!(target_col %in% colnames(df))) {
        stop("The specified target column is not present in the data frame.")
      }

      target <- df[[target_col]]
      predictors <- df[, colnames(df) != target_col, drop = FALSE]

      remove_collinear_vif <- function(target, predictors, vif_threshold) {
        cleaned_predictors <- predictors

        repeat {
          model <- try(stats::lm(target ~ ., data = data.frame(cleaned_predictors)), silent = TRUE)

          if (inherits(model, "try-error")) break

          vif_values <- try(car::vif(model), silent = TRUE)

          if (inherits(vif_values, "try-error")) break

          max_vif <- max(vif_values, na.rm = TRUE)
          if (max_vif > vif_threshold) {
            col_to_remove <- names(which.max(vif_values))
            message(paste("Removing variable due to high VIF:", col_to_remove))
            cleaned_predictors <- cleaned_predictors[, !colnames(cleaned_predictors) %in% col_to_remove, drop = FALSE]
          } else {
            break
          }
        }

        return(cleaned_predictors)
      }

      remove_collinear_correlation <- function(predictors, correlation_threshold) {
        cor_matrix <- stats::cor(predictors, use = "pairwise.complete.obs")
        high_cor <- which(abs(cor_matrix) > correlation_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)
        high_cor <- high_cor[high_cor[, 1] < high_cor[, 2],, drop = FALSE]
        cols_to_remove <- unique(high_cor[, 2])
        return(predictors[, - cols_to_remove, drop = FALSE])
      }

      try_vif <- try(remove_collinear_vif(target, predictors, vif_threshold), silent = TRUE)

      if (!inherits(try_vif, "try-error") && ncol(try_vif) > 0) {
        cleaned_predictors <- try_vif
      } else {
        message("VIF-based collinearity removal failed. Falling back to correlation-based method.")
        cleaned_predictors <- remove_collinear_correlation(predictors, correlation_threshold)
      }

      cleaned_df <- cbind(Target = target, cleaned_predictors)
      return(cleaned_df)
    }

    if (remove_colinear_vars && !is.null(Target_unique)) {
      df_unique_no_labels <- df_unique[, setdiff(names(df_unique), c("Label")), drop = FALSE]
      target_column <- "Target"

      cleaned_data <- remove_collinear_with_fallback(df_unique_no_labels,
                                                      target_col = target_column,
                                                      vif_threshold = vif_threshold_for_colinearity,
                                                      correlation_threshold = correlation_threshold_for_colinearity)

      X_unique <- cleaned_data[, setdiff(names(cleaned_data), c("Target")), drop = FALSE]
      df_unique <- cbind.data.frame(cleaned_data, Label = labels)
      Target_unique <- cleaned_data$Target
    }

    # COMPLETE get_recommended_dimensions function
    get_recommended_dimensions <- function(data, variance_threshold = 0.8, correlation_threshold = 0.7) {
      numeric_data <- data

      pca_result <- prcomp(numeric_data, scale. = TRUE)

      var_explained <- cumsum(pca_result$sdev ^ 2 / sum(pca_result$sdev ^ 2))
      var_dims <- which(var_explained >= variance_threshold)[1]

      eigenvalues <- pca_result$sdev ^ 2
      kaiser_dims <- sum(eigenvalues > 1)

      cor_matrix <- cor(numeric_data)
      high_cor_count <- sum(abs(cor_matrix[upper.tri(cor_matrix)]) > correlation_threshold)
      cor_dims <- ncol(numeric_data) - high_cor_count

      recommended_dims <- max(2, ceiling(mean(c(var_dims, kaiser_dims, cor_dims))))

      return(recommended_dims)
    }

    ncomp <- get_recommended_dimensions(X_unique)

  } else {
    # TRANSFORM MODE
    if (is.null(projection_object)) stop("projection_object required when newdata is provided")

    # Apply same preprocessing as training
    if (scaleX && !projection_method %in% c("IPCA", "NMF")) {
      if (is.null(attr(X, "scaled:scale"))) {
        X <- scale(X)
      }
      newdata <- scale(newdata, center = attr(X, "scaled:center"), scale = attr(X, "scaled:scale"))
    }

    X_unique <- newdata
    if (!is.null(labels) && length(labels) == nrow(newdata)) {
      df_unique <- data.frame(newdata, Label = labels)
    } else {
      df_unique <- data.frame(newdata, Label = seq_len(nrow(newdata)))
    }
    Target_unique <- Target
    ncomp <- ncol(projection_object$Projected)
  }

  # Validate Target for supervised methods
  if (projection_method %in% c("PLSDA", "LDA", "PLSLDA") && is.null(Target_unique) && is.null(newdata)) {
    stop(paste(projection_method, "requires Target to be provided."))
  }

  set.seed(seed)
  Proj <- NULL
  local_projection_object <- projection_object

  # Helper functions needed by specific methods
  sort_ica_components <- function(S, A) {
    var_explained <- apply(S, 2, stats::var)
    order_components <- order(var_explained, decreasing = TRUE)
    S_sorted <- S[, order_components]
    A_sorted <- A[, order_components]
    return(list(S = S_sorted, A = A_sorted))
  }

  # MAIN PROJECTION SWITCH
  switch(projection_method,
          "none" = {
    Proj <- X_unique
  },
          "PCA" = {
    if (is.null(newdata)) {
      pca_obj <- FactoMineR::PCA(X_unique, graph = FALSE, ncp = ncomp)
      Proj <- pca_obj$ind$coord
      local_projection_object <- pca_obj
    } else {
      Proj <- predict(local_projection_object, newdata = X_unique)$ind$coord[, 1:min(ncomp, ncol(Proj))]
    }
  },
          "ICA" = {
    if (is.null(newdata)) {
      ica_obj <- fastICA::fastICA(X_unique, n.comp = ncomp, alg.typ = "parallel",
                                           fun = "logcosh", alpha = 1, method = "C",
                                           row.norm = FALSE, maxit = 1000, tol = 1e-04, verbose = FALSE)
      sorted_ica <- sort_ica_components(S = ica_obj$S, A = ica_obj$A)
      Proj <- sorted_ica$S
      local_projection_object <- ica_obj
    } else {
      Proj <- as.matrix(X_unique) %*% local_projection_object$A[, 1:ncomp]
    }
  },
          "MDS" = {
    if (is.null(newdata)) {
      d <- dist(X_unique)
      mds_obj <- MASS::isoMDS(d, k = ncomp, maxit = 1000, tol = 1e-4)
      Proj <- mds_obj$points
      local_projection_object <- mds_obj
    } else {
      # Procrustes alignment to training MDS
      d_new <- dist(X_unique)
      mds_new <- MASS::isoMDS(d_new, init = local_projection_object$points[1:min(nrow(local_projection_object$points), nrow(X_unique)),], k = ncomp)
      Proj <- mds_new$points
    }
  },
          "tSNE" = {
    if (is.null(newdata)) {
      perpl <- ifelse(nrow(X_unique) / 3 - 1 / 3 < 30, nrow(X_unique) / 3 - 0.34, 30)
      tsne_obj <- Rtsne::Rtsne(X_unique, dims = min(ncomp, 3), perplexity = perpl)
      Proj <- tsne_obj$Y
      local_projection_object <- tsne_obj
    } else {
      perpl <- local_projection_object$perplexity
      combined_data <- rbind(local_projection_object$X, X_unique)
      tsne_full <- Rtsne::Rtsne(combined_data, dims = min(ncomp, 3), perplexity = perpl)
      Proj <- tsne_full$Y[(nrow(local_projection_object$X) + 1):nrow(combined_data),]
    }
  },
          "isomap" = {
    if (is.null(newdata)) {
      isomap_obj <- RDRToolbox::Isomap(as.matrix(X_unique), dims = ncomp, k = 5)
      Proj <- isomap_obj[[1]]
      local_projection_object <- isomap_obj
    } else {
      isomap_obj <- RDRToolbox::Isomap(as.matrix(X_unique), dims = ncomp, k = 5)
      Proj <- isomap_obj[[1]]
    }
  },
          "PLSDA" = {
    if (is.null(newdata)) {
      plsda_obj <- mixOmics::plsda(X_unique, factor(Target_unique))
      plot_df <- mixOmics::plotIndiv(plsda_obj)$df
      Proj <- plot_df[c("x", "y")]
      local_projection_object <- plsda_obj
    } else {
      pred <- predict(local_projection_object, newdata = X_unique)
      Proj <- data.frame(x = pred$variates$X[, 1], y = pred$variates$X[, 2])
    }
  },
          "LDA" = {
    if (is.null(newdata)) {
      dfLDA <- cbind.data.frame(Target_unique = Target_unique, X_unique)
      lda_obj <- MASS::lda(factor(Target_unique) ~ ., data = dfLDA, dimen = ncomp)
      lda_proj <- predict(lda_obj, X_unique)
      Proj <- lda_proj$x
      local_projection_object <- lda_obj
    } else {
      lda_proj <- predict(local_projection_object, newdata = X_unique)
      Proj <- lda_proj$x[, 1:min(ncomp, ncol(lda_proj$x))]
    }
  },
          "PLSLDA" = {
    pls_lda_projection <- function(X, y, ncomp) {
      dfPLSLDA <- cbind.data.frame(y = y, X)
      pls_model <- pls::plsr(y ~ ., data = dfPLSLDA, ncomp = ncomp, scale = TRUE)
      pls_scores <- pls::scores(pls_model)[, 1:ncomp]
      lda_model <- MASS::lda(pls_scores, grouping = y)
      lda_projection <- predict(lda_model, newdata = pls_scores)$x

      return(list(pls_model = pls_model, lda_model = lda_model,
                            pls_scores = pls_scores, lda_projection = lda_projection))
    }

    if (is.null(newdata)) {
      results <- pls_lda_projection(X = as.data.frame(X_unique), y = Target_unique, ncomp = ncomp)
      Proj <- results$lda_projection
      local_projection_object <- results
    } else {
      pls_scores <- predict(local_projection_object$pls_model, newdata = X_unique)[, 1:ncomp]
      lda_proj <- predict(local_projection_object$lda_model, newdata = pls_scores)$x
      Proj <- lda_proj
    }
  },
          "IPCA" = {
    if (is.null(newdata)) {
      ipca_obj <- mixOmics::ipca(as.matrix(X_unique), ncomp = ncomp, mode = "deflation",
                                          fun = "logcosh", scale = F, w.init = NULL,
                                          max.iter = 1000, tol = 1e-04)
      plot_df <- mixOmics::plotIndiv(ipca_obj)$df
      Proj <- plot_df[c("x", "y")]
      if (switchDims && ipca_obj$prop_expl_var$X[1] < ipca_obj$prop_expl_var$X[2]) {
        Proj <- plot_df[c("y", "x")]
      }
      local_projection_object <- ipca_obj
    } else {
      pred <- predict(local_projection_object, newdata = X_unique)
      Proj <- data.frame(x = pred$variates$X[, 1], y = pred$variates$X[, 2])
    }
  },
          "Umap" = {
    if (is.null(newdata)) {
      umap_obj <- uwot::umap(X_unique, ret_model = TRUE)
      Proj <- umap_obj$embedding
      local_projection_object <- umap_obj
    } else {
      Proj <- uwot::umap_transform(X_unique, local_projection_object)
    }
  },
          "NMF" = {
    perform_nmf <- function(data, ncomp) {
      if (any(data <= 0)) {
        message("Input matrix contains non-positive values. Shifting all values to be positive.")
        min_value <- min(data)
        data <- data - min_value + 1e-10
      }
      nmf_result <- NMF::nmf(data, ncomp, .options = list(p = FALSE))
      coef_matrix <- NMF::basis(nmf_result)
      return(coef_matrix)
    }

    if (is.null(newdata)) {
      nmf_result <- suppressMessages(perform_nmf(data = as.matrix(X_unique), ncomp = ncomp))
      Proj <- nmf_result
      local_projection_object <- nmf_result
    } else {
      Proj <- suppressMessages(perform_nmf(data = as.matrix(X_unique), ncomp = ncomp))
    }
  },
          "LLE" = {
    lle_optimized <- function(data, max_dim = 10, max_k = 30, perplexity = 30) {
      X <- as.matrix(data)
      X <- scale(X)

      lle <- function(X, k, d, reg = 1e-3) {
        n <- nrow(X)
        p <- ncol(X)

        nn <- RANN::nn2(X, k = k + 1)
        neighbors <- nn$nn.idx[, -1]

        W <- matrix(0, n, n)
        for (i in 1:n) {
          z <- X[neighbors[i,],, drop = FALSE] - matrix(X[i,], k, p, byrow = TRUE)
          C <- tcrossprod(z)
          C <- C + diag(k) * reg * sum(diag(C))
          w <- solve(C, rep(1, k))
          w <- w / sum(w)
          W[i, neighbors[i,]] <- w
        }

        M <- diag(n) - W - t(W) + t(W) %*% W
        eig <- eigen(M, symmetric = TRUE)
        Y <- eig$vectors[, (n - d):(n - 1), drop = FALSE]

        var_explained <- apply(Y, 2, stats::var)
        order <- order(var_explained, decreasing = TRUE)
        Y <- Y[, order, drop = FALSE]

        return(list(Y = Y, var_explained = var_explained[order]))
      }

      optimize_params <- function(X, max_k, max_d, perplexity) {
        scores <- matrix(Inf, max_k, max_d)
        ncomp_local <- min(ncol(X), 3)

        tsne_emb <- Rtsne::Rtsne(X, dims = ncomp_local, perplexity = perplexity, verbose = FALSE)$Y

        for (k in seq(5, max_k, by = 5)) {
          for (d in 1:min(max_d, k - 1)) {
            tryCatch({
              result <- lle(X, k, d)
              Y <- result$Y

              if (ncol(Y) < 2) {
                Y <- cbind(Y, rep(0, nrow(Y)))
              }

              kl_div <- KL_divergence(tsne_emb, Y[, 1:2, drop = FALSE])
              trust <- trustworthiness(X, Y, k = min(10, k - 1))
              scores[k, d] <- kl_div - trust

            }, error = function(e) {
              message("Error for k=", k, ", d=", d, ": ", e$message)
            })
          }
        }

        best <- which(scores == min(scores, na.rm = TRUE), arr.ind = TRUE)[1,]
        return(list(k = best[1], d = best[2]))
      }

      KL_divergence <- function(P, Q) {
        P <- as.matrix(dist(P))
        Q <- as.matrix(dist(Q))
        P <- P / sum(P)
        Q <- Q / sum(Q)
        return(sum(P * log((P + 1e-10) / (Q + 1e-10))))
      }

      trustworthiness <- function(X, Y, k) {
        n <- nrow(X)
        r_x <- apply(as.matrix(dist(X)), 2, rank)
        r_y <- apply(as.matrix(dist(Y)), 2, rank)
        sum_penalty <- sum(sapply(1:n, function(i) {
          sum(pmax(r_x[i, r_y[i,] <= k] - k, 0))
        }))
        return(1 - 2 / (n * k * (2 * n - 3 * k - 1)) * sum_penalty)
      }

      perpl <- ifelse(nrow(X) / 3 - 1 / 3 < 30, nrow(X) / 3 - 0.34, 30)

      params <- optimize_params(X, max_k, max_dim, perplexity = perpl)
      params$d <- max(2, params$d)
      result <- lle(X, params$k, params$d)

      return(list(projection = result$Y, k = params$k, d = params$d, var_explained = result$var_explained))
    }

    if (is.null(newdata)) {
      result_lle <- lle_optimized(X_unique)
      Proj <- result_lle$projection
      local_projection_object <- result_lle
    } else {
      result_lle <- lle_optimized(X_unique)
      Proj <- result_lle$projection
    }
  },
          "autoencoder" = {
    AutoEncode2 <- function(Data, Hidden = NULL, Learningrate = c(0.01, 0.1, 0.5),
                                     Threshold = c(0.01, 0.1, 0.5), n_times_n_features = 3,
                                     Epochs = 10000, ActivationFunction = c("logistic", "tanh"), ...) {
      if (!is.numeric(Data) || any(is.na(Data))) {
        stop("Data must be numeric and contain no missing values.")
      }

      if (!is.null(Hidden) && (length(Hidden) %% 2) != 1)
        stop("AutoEncode: Number of Layers of Hidden Neurons has to be odd")

      n_features <- ncol(Data)

      if (is.null(Hidden)) {
        bottleneck_size <- max(2, round(n_times_n_features * n_features))
        middle_layer <- bottleneck_size / 2
        Hidden <- c(bottleneck_size, middle_layer, bottleneck_size)
      }

      d <- as.data.frame(Data)
      formel <- paste0(paste(colnames(d), collapse = " + "), " ~ ",
                                paste(colnames(d), collapse = " + "))

      param_grid <- expand.grid(
                Learningrate = Learningrate,
                Threshold = Threshold
              )

      train_model <- function(Learningrate, Threshold, ActivationFunction) {
        tryCatch({
          net <- neuralnet::neuralnet(formel, d, hidden = Hidden, stepmax = Epochs,
                                               lifesign = "none", linear.output = FALSE,
                                               algorithm = "backprop", learningrate = Learningrate,
                                               threshold = Threshold)

          x <- neuralnet::compute(net, d)
          reproduction <- x$net.result
          mse <- mean((as.matrix(Data) - reproduction) ^ 2)
          return(list(net = net, mse = mse))
        }, error = function(e) {
          return(list(net = NULL, mse = Inf))
        })
      }

      results <- apply(param_grid, 1, function(params) {
        train_model(params[1], params[2], ActivationFunction[1])
      })

      best_model <- results[[which.min(sapply(results, function(x) x$mse))]]

      if (is.null(best_model$net)) {
        stop("Failed to train any valid model. Try adjusting the parameter ranges.")
      }

      net <- best_model$net
      x <- neuralnet::compute(net, d)
      ProjectionLayer <- ceiling((length(Hidden) + 2) / 2)
      nrOfDim <- ncol(x$neurons[[ProjectionLayer]]) - 1
      projection <- x$neurons[[ProjectionLayer]][, 2:(nrOfDim + 1), drop = FALSE]

      var_importance <- apply(projection, 2, stats::var)
      ordered_projection <- projection[, order(var_importance, decreasing = TRUE), drop = FALSE]

      return(ordered_projection)
    }

    if (is.null(newdata)) {
      Proj <- AutoEncode2(as.matrix(X_unique),
                                   Learningrate = c(0.01, 0.1, 0.5),
                                   Threshold = c(0.01, 0.1, 0.5),
                                   ActivationFunction = c("logistic", "tanh"),
                                   Epochs = 5000,
                                   n_times_n_features = 3)
      local_projection_object <- "trained_autoencoder"
    } else {
      Proj <- AutoEncode2(as.matrix(X_unique),
                                   Learningrate = c(0.01, 0.1, 0.5),
                                   Threshold = c(0.01, 0.1, 0.5),
                                   ActivationFunction = c("logistic", "tanh"),
                                   Epochs = 1000,
                                   n_times_n_features = 3)
    }
  }, {
    Proj <- X_unique
  }
  )

  # Ensure that at least two columns are returned
  if (ncol(Proj) == 1) {
    if (jitter_one_dimensional_projections) {
      Proj <- cbind.data.frame(Proj, jitter(Proj))
    } else {
      Proj <- cbind.data.frame(Proj, Proj)
    }
  }

  Proj <- data.frame(Proj)
  names(Proj) <- paste0("Dim", seq_len(ncol(Proj)))

  # Return the projection, unique data frame, and projection object
  return(list(Projected = Proj,
                UniqueData = df_unique,
                projection_object = local_projection_object))
}

############### Functions for clustering ###############

# Apply clustering algorithms
apply_clustering <- function(projection_results, clustering_methods, cluster_number_methods, method_for_hcpc,
                              distance_metric = "euclidean", seed = 42, nProc = 12,
                              max_clusters = max_clusters) {

  cluster_list <- pbmcapply::pbmclapply(names(projection_results), function(projection) {
    projection_result <- projection_results[[projection]]
    tryCatch({
      if (is.null(projection_result$UniqueData$Target)) {
        stop("Required data is missing for clustering.")
      }
      projected_data <- as.data.frame(projection_result$Projected)
      clusters_per_projection <- lapply(clustering_methods, function(cluster_alg) {
        clusters_results <- lapply(cluster_number_methods, function(cluster_number_method) {
          set.seed(seed)
          clusters <- performClustering(X = projected_data,
                                         Target = projection_result$UniqueData$Target,
                                         clustering_method = cluster_alg,
                                         cluster_number_method = cluster_number_method,
                                         method_for_hcpc = method_for_hcpc,
                                         distance_metric = distance_metric,
                                         max_clusters = max_clusters)
          if (cluster_alg != "none") clusters <- renameClusters(trueCls = projection_result$UniqueData$Target, clusters)
          combined_result <- create_combined_result(projected_data, clusters, projection_result)
          return(combined_result)
        })
        names(clusters_results) <- cluster_number_methods
        return(clusters_results)
      })
      names(clusters_per_projection) <- clustering_methods
      return(clusters_per_projection)

    }, error = function(err) {
      message("Error in clustering with projection ", projection, ": ", err$message)
      return(NULL)
    })
  }, mc.cores = min(length(projection_results), nProc))

  names(cluster_list) <- names(projection_results)
  return(cluster_list)
}


# Function to determine the number of clusters
findOptimalClusters <- function(X, clustering_method, cluster_number_method, Target = NULL,
                                 method_for_hcpc = "ward", distance_metric = "euclidean", max_clusters = max_clusters) {

  # Load necessary libraries
  if (!requireNamespace("NbClust", quietly = TRUE)) {
    stop("Package 'NbClust' is required but not installed.")
  }
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("Package 'FactoMineR' is required but not installed.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but not installed.")
  }

  # Function to calculate mode
  calculateMode <- function(v) {
    if (length(v) == 0) return(NA)
    freqTable <- table(v)
    modes <- as.numeric(names(freqTable)[freqTable == max(freqTable)])
    if (length(modes) > 1) modes else modes[1]
  }

  nClusters <- switch(cluster_number_method,
                       "NbClust" = {
    # Adjust clustering method for compatibility with NbClust
    clustering_method <- switch(clustering_method,
                                                      "kmedoids" = "kmeans",
                                                      "diana" = "ward.D2",
                                                      "hcpc" = "ward",
                                                      clustering_method
                         )

    # Try running NbCLust

    res <- try(suppressWarnings(NbClust::NbClust(data = X, diss = NULL, distance = distance_metric,
                                                                         min.nc = 2, max.nc = max_clusters, method = clustering_method)), silent = TRUE)

    if (!inherits(res, "try-error")) {
      max(unlist(res[4]))
    } else {
      # Try running NbClust using parallel processing
      nbclustIndices <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
                                                "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
                                                "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                                "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")

      # Determine number of clusters by running NbClust method by method
      nClustersTest <- lapply(nbclustIndices, function(i) {
        res <- try(suppressWarnings(NbClust::NbClust(data = X, diss = NULL, distance = distance_metric,
                                                                             min.nc = 2, max.nc = max_clusters, method = clustering_method, index = i)), silent = TRUE)
        if (!inherits(res, "try-error")) {
          max(unlist(res[4]))
        } else {
          warning(paste("Error with index", i))
          NA
        }
      })

      # Compute mode of results to determine optimal number of clusters
      nC <- na.omit(unlist(nClustersTest))
      result <- if (length(nC) > 0) calculateMode(nC) else NA
      result[length(result)]
    }
  },
                       "hcpc" = {
    res.hcpc <- FactoMineR::HCPC(X, consol = TRUE, method = method_for_hcpc, metric = distance_metric, nb.clust = -1,
                                                       iter.max = 100, graph = FALSE, nstart = 100)
    length(unique(res.hcpc$data.clust$clust))
  },
  # Default case if Target is provided or none matches
                       if (!is.null(Target)) {
                        length(unique(Target))
                       } else {
    1
  }
  )

  if (is.na(nClusters)) {
    warning("Unable to determine the number of clusters. Setting to 1 or the original classes if provided.")
    if (!is.null(Target)) {
      nClusters <- length(unique(Target))
    } else {
      nClusters <- 1
    }
  }

  return(nClusters)
}


# Function to perform clustering
performClustering <- function(X, Target = NULL, clustering_method = "kmeans", cluster_number_method = "orig", method_for_hcpc = "ward",
                               distance_metric = "euclidean", max_clusters = max_clusters) {

  # Check if input X is a non-empty data frame
  if (!is.data.frame(X) || nrow(X) == 0) stop("Input data must be a non-empty data frame.")

  # Define valid clustering methods
  valid_clustering_methods <- c("none", "kmeans", "kmedoids",
                                 "diana", "hcpc", "ward.D2", "single", "average", "median", "complete", "centroid")
  valid_hcpc_methods <- c("average", "single", "complete", "ward", "weighted", "flexible", "gaverage")
  valid_distances <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  valid_cluster_number_methods <- c("orig", "NbClust", "hcpc")

  # Validate the selected method
  if (!(clustering_method %in% valid_clustering_methods))
    stop(paste("Invalid clustering method. Choose from:", paste(valid_clustering_methods, collapse = ", ")))
  if (!(method_for_hcpc %in% valid_hcpc_methods))
    stop(paste("Invalid method for HCPC clustering. Choose from:", paste(valid_hcpc_methods, collapse = ", ")))
  if (!(distance_metric %in% valid_distances))
    stop(paste("Invalid distance metric. Choose from:", paste(valid_distances, collapse = ", ")))
  if (!(cluster_number_method %in% valid_cluster_number_methods))
    stop(paste("Invalid method for cluster number determination. Choose from:", paste(valid_cluster_number_methods, collapse = ", ")))

  # If Target is NULL, initialize it with default values
  if (is.null(Target)) Target <- rep(1, nrow(X))

  # Determine the number of clusters using the custom function

  if (clustering_method == "none") {
    nClusters <- length(unique(Target))
  } else {
    nClusters <- findOptimalClusters(X, clustering_method = clustering_method,
                                      cluster_number_method = cluster_number_method, Target = Target, max_clusters = max_clusters)
  }

  # Perform clustering based on the specified method
  cluster_result <- switch(clustering_method,
                            "none" = Target, # Return Target directly if method is "none"
                            "kmeans" = {
    # Perform K-means clustering
    kmeans(X, centers = nClusters, nstart = 100,)$cluster
  },
                            "kmedoids" = {
    # Perform K-medoids clustering using PAM
    cluster::pam(X, k = nClusters)$clustering
  },
                            "diana" = {
    # Perform divisive clustering
    clusterObject <- cluster::diana(X)
    as.numeric(cutree(clusterObject, k = nClusters))
  },
                            "hcpc" = {
    # Perform hierarchical clustering and principal component analysis
    res.pca <- FactoMineR::PCA(X, graph = FALSE)
    FactoMineR::HCPC(res.pca, nb.clust = nClusters, metric = distance_metric,
                                                graph = FALSE, nstart = 100, method = method_for_hcpc)$data.clust$clust
  }, {
    # Perform other hierarchical clustering methods
    clusterObject <- stats::hclust(dist(X, method = distance_metric), method = clustering_method)
    as.numeric(cutree(clusterObject, k = nClusters))
  })

  # Return the clustering result
  return(cluster_result)
}


# Function to align cluster names with likely class labels of the prior classification
renameClusters <- function(trueCls, currentCls, K = 12) {

  # Helper function to reduce clusters
  ReduceClsToK <- function(Cls, K = 12) {
    uniqueCls <- unique(Cls)
    while (length(uniqueCls) > K) {
      counts <- table(Cls)
      to_merge <- names(sort(counts)[1:2])
      Cls[Cls %in% to_merge] <- to_merge[1]
      uniqueCls <- unique(Cls)
    }
    return(Cls)
  }

  # Preprocess input
  trueCls[!is.finite(trueCls)] <- 9999
  currentCls[!is.finite(currentCls)] <- 9999

  # Warnings
  if (length(unique(trueCls)) > 9) {
    warning("Too many clusters in PriorCls. Consider using cloud computing.")
  }
  if (length(unique(currentCls)) > K) {
    warning("Too many clusters in CurrentCls. Combining smallest clusters.")
    currentCls <- ReduceClsToK(currentCls, K)
  }

  # Normalize clusters
  trueCls <- as.numeric(factor(trueCls))
  currentCls <- as.numeric(factor(currentCls))

  # Get unique labels
  uniqueLabels <- sort(unique(c(trueCls, currentCls)))
  nLabels <- length(uniqueLabels)

  # Generate permutations
  permutations <- combinat::permn(nLabels)

  # Find best permutation
  bestAccuracy <- 0
  bestPermutation <- NULL

  for (perm in permutations) {
    newLabels <- perm[match(currentCls, seq_along(perm))]
    accuracy <- sum(trueCls == newLabels) / length(trueCls)
    if (accuracy > bestAccuracy) {
      bestAccuracy <- accuracy
      bestPermutation <- perm
    }
  }

  # Rename clusters
  renamedCls <- bestPermutation[match(currentCls, seq_along(bestPermutation))]

  return(renamedCls)
}


# Function to calculate Adjusted Rand Index (ARI) without using mclust
calculate_adjusted_rand_index <- function(Target, Clusters) {
  # Helper function to compute a contingency table
  contingency_table <- function(x, y) {
    table(factor(x, levels = unique(c(x, y))), factor(y, levels = unique(c(x, y))))
  }

  # Compute the contingency table
  cont_table <- contingency_table(Target, Clusters)

  # Sum over rows and columns
  sum_comb_n <- function(n) sum(choose(n, 2))

  # Total number of pairs
  N <- sum(cont_table)

  # Sum of combinations for each cell
  sum_comb_c_ij <- sum(mapply(choose, as.numeric(cont_table), 2))

  # Sum of combinations for each row
  sum_comb_a_i <- sum_comb_n(rowSums(cont_table))

  # Sum of combinations for each column
  sum_comb_b_j <- sum_comb_n(colSums(cont_table))

  # Calculate ARI
  expected_index <- sum_comb_a_i * sum_comb_b_j / choose(N, 2)
  max_index <- (sum_comb_a_i + sum_comb_b_j) / 2
  ari <- (sum_comb_c_ij - expected_index) / (max_index - expected_index)

  return(ari)
}


# Function to calculate Dunn's Index using a distance matrix
calculate_dunn_index <- function(distance_matrix, clusters) {
  # Convert the distance vector to a complete distance matrix
  distance_matrix <- as.matrix(distance_matrix)

  # Number of clusters
  unique_clusters <- unique(clusters)
  k <- length(unique_clusters)

  # Initialize variables for inter-cluster and intra-cluster distances
  min_inter_cluster_dist <- Inf
  max_intra_cluster_dist <- 0

  # Calculate inter-cluster distances
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      cluster_i_indices <- which(clusters == unique_clusters[i])
      cluster_j_indices <- which(clusters == unique_clusters[j])

      for (p in cluster_i_indices) {
        for (q in cluster_j_indices) {
          dist <- distance_matrix[p, q]
          if (!is.na(dist) && dist < min_inter_cluster_dist) {
            min_inter_cluster_dist <- dist
          }
        }
      }
    }
  }

  # Calculate intra-cluster distances
  for (i in 1:k) {
    cluster_indices <- which(clusters == unique_clusters[i])
    if (length(cluster_indices) > 1) {
      # Ensure there's more than one point in the cluster
      for (p in 1:(length(cluster_indices) - 1)) {
        for (q in (p + 1):length(cluster_indices)) {
          dist <- distance_matrix[cluster_indices[p], cluster_indices[q]]
          if (!is.na(dist) && dist > max_intra_cluster_dist) {
            max_intra_cluster_dist <- dist
          }
        }
      }
    }
  }

  # Ensure the max_intra_cluster_dist is valid for division
  if (max_intra_cluster_dist == 0) {
    stop("Intra-cluster distance is zero, which might indicate that each cluster has a single point.")
  }

  # Calculate Dunn's Index
  dunn_index <- min_inter_cluster_dist / max_intra_cluster_dist

  return(dunn_index)
}


# Calculate the Davies-Bold cluster index
calculate_davies_bouldin_index <- function(clusters, data) {

  # Number of clusters
  k <- length(unique(clusters))

  # Ensure data is numeric
  numeric_data <- data[, sapply(data, is.numeric)]

  # Calculate centroids for each cluster
  centroids <- aggregate(numeric_data, by = list(cluster = clusters), FUN = mean)
  rownames(centroids) <- centroids$cluster
  centroids <- centroids[, -1] # Remove the cluster column

  # Calculate average distance within clusters
  avg_within_dist <- tapply(seq_len(nrow(data)), clusters, function(idx) {
    cluster_data <- data[idx,]
    centroid <- centroids[as.character(clusters[idx[1]]),]
    mean(apply(cluster_data, 1, function(x) sqrt(sum((x - centroid) ^ 2))))
  })

  emicheps <- .Machine$double.eps ^ 0.5

  # Calculate Davies-Bouldin index
  db_index <- mean(sapply(1:k, function(i) {
    max_ratio <- 0
    for (j in 1:k) {
      if (i != j) {
        dist_between <- sqrt(sum((centroids[i,] - centroids[j,]) ^ 2))
        ratio <- (avg_within_dist[i] + avg_within_dist[j]) / (dist_between + emicheps)
        max_ratio <- max(max_ratio, ratio)
      }
    }
    max_ratio
  }))

  return(db_index)
}


# Calculate the Density-Based Clustering Validation (DBCV) index
calculate_dbcv <- function(data, clusters, dist_method = "euclidean") {

  # Calculate the DBCV index (using a simplified heuristic approach)
  intra_density <- numeric(max(clusters))
  inter_density <- numeric(max(clusters))

  for (i in unique(clusters)) {
    cluster_points <- data[clusters == i,]
    if (nrow(cluster_points) > 1) {
      intra_cluster_dist_matrix <- dist(cluster_points, method = dist_method)
      intra_density[i] <- mean(intra_cluster_dist_matrix)
    }
  }

  for (i in unique(clusters)) {
    for (j in unique(clusters)) {
      if (i < j) {
        combined_points <- rbind(data[clusters == i,], data[clusters == j,])
        inter_cluster_dist_matrix <- dist(combined_points, method = dist_method)
        inter_density[i] <- mean(inter_cluster_dist_matrix)
      }
    }
  }

  dbcv_index <- sum(intra_density) / (sum(intra_density) + sum(inter_density))

  return(dbcv_index)
}


# Calinski-Harabasz Index
calculate_calinski_harabasz_index <- function(data, cluster_labels) {
  n <- nrow(data)
  k <- length(unique(cluster_labels))

  # Calculate overall mean
  overall_mean <- colMeans(data)

  # Calculate between-cluster sum of squares
  between_ss <- 0
  for (cluster in unique(cluster_labels)) {
    cluster_data <- data[cluster_labels == cluster,]
    cluster_mean <- colMeans(cluster_data)
    between_ss <- between_ss + nrow(cluster_data) * sum((cluster_mean - overall_mean) ^ 2)
  }

  # Calculate within-cluster sum of squares
  within_ss <- 0
  for (cluster in unique(cluster_labels)) {
    cluster_data <- data[cluster_labels == cluster,]
    cluster_mean <- colMeans(cluster_data)
    within_ss <- within_ss + sum(apply(cluster_data, 1, function(x) sum((x - cluster_mean) ^ 2)))
  }

  # Calculate Calinski-Harabasz Index
  ch_index <- ((n - k) / (k - 1)) * (between_ss / within_ss)
  return(ch_index)
}

# Inertia (Within-Cluster Sum of Squares)
calculate_inertia <- function(data, cluster_labels) {
  # Ensure data is a matrix
  data <- as.matrix(data)

  # Calculate cluster centers
  unique_labels <- unique(cluster_labels)
  centers <- t(sapply(unique_labels, function(k) {
    colMeans(data[cluster_labels == k,, drop = FALSE])
  }))

  # Calculate inertia
  total_inertia <- sum(sapply(seq_along(unique_labels), function(i) {
    cluster_data <- data[cluster_labels == unique_labels[i],, drop = FALSE]
    sum(apply(cluster_data, 1, function(row) sum((row - centers[i,]) ^ 2)))
  }))

  return(total_inertia)
}


# Normalized Mutual Information
# Helper functions
entropy <- function(labels) {
  p <- table(labels) / length(labels)
  - sum(p * log(p))
}

calculate_nmi <- function(true_labels, cluster_labels) {
  contingency <- table(true_labels, cluster_labels)

  # Calculate entropies
  h_true <- entropy(true_labels)
  h_cluster <- entropy(cluster_labels)

  # Calculate mutual information
  n <- sum(contingency)
  mi <- sum(contingency * log(contingency * n / (rowSums(contingency) %*% t(colSums(contingency)))))

  # Calculate NMI
  nmi <- 2 * mi / (h_true + h_cluster)
  return(nmi)
}

# Adjusted Mutual Information
mutual_information <- function(true_labels, cluster_labels) {
  contingency <- table(true_labels, cluster_labels)
  n <- sum(contingency)
  mi <- 0
  for (i in seq_len(nrow(contingency))) {
    for (j in seq_len(ncol(contingency))) {
      if (contingency[i, j] > 0) {
        mi <- mi + contingency[i, j] * log(contingency[i, j] * n / (sum(contingency[i,]) * sum(contingency[, j])))
      }
    }
  }
  mi / n
}

expected_mutual_information <- function(true_labels, cluster_labels) {
  contingency <- table(true_labels, cluster_labels)
  n <- sum(contingency)
  a <- rowSums(contingency)
  b <- colSums(contingency)
  emi <- 0
  for (i in seq_len(nrow(contingency))) {
    for (j in seq_len(ncol(contingency))) {
      if (a[i] > 0 && b[j] > 0) {
        mij <- max(1, a[i] + b[j] - n)
        xij <- min(a[i], b[j])
        for (nij in mij:xij) {
          v <- nij / n * log((nij * n) / (a[i] * b[j]))
          emi <- emi + v * exp(lchoose(a[i], nij) + lchoose(n - a[i], b[j] - nij) - lchoose(n, b[j]))
        }
      }
    }
  }
  emi
}

calculate_ami <- function(true_labels, cluster_labels) {
  mi <- mutual_information(true_labels, cluster_labels)
  emi <- expected_mutual_information(true_labels, cluster_labels)
  h_true <- entropy(true_labels)
  h_cluster <- entropy(cluster_labels)
  max_h <- max(h_true, h_cluster)

  if (mi == emi) {
    return(0)
  } else {
    return((mi - emi) / (max_h - emi))
  }
}

# Calculate various cluster indices and return as a data frame
calculate_cluster_scores_df <- function(clustering_results, projection_methods, clustering_methods, cluster_number_methods,
                                         distance_metric = "euclidean", nProc = 2) {

  # Validate input
  if (!distance_metric %in% c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")) {
    stop("Unsupported distance metric")
  }

  cluster_scores <- pbmcapply::pbmclapply(projection_methods, function(projection_method) {
    if (is.null(clustering_results[[projection_method]])) {
      warning("Projection for method ", projection_method, " returned NULL. Skipping.")
      return(NULL)
    }

    cluster_scores_per_alg <- lapply(clustering_methods, function(cluster_alg) {
      cluster_scores_per_number_method <- lapply(cluster_number_methods, function(cluster_number_method) {

        combined_result <- clustering_results[[projection_method]][[cluster_alg]][[cluster_number_method]]

        projected_data <- combined_result[, 1:2]
        Target <- combined_result$Target
        Clusters <- combined_result$Cluster


        distance_matrix <- stats::dist(projected_data, method = distance_metric)

        cluster_accuracy <- sum(Clusters == Target) / length(Clusters)
        Silhouette_index <- mean(cluster::silhouette(Clusters, distance_matrix)[, "sil_width"])
        Dunn_index <- calculate_dunn_index(distance_matrix = distance_matrix, clusters = Clusters)
        Rand_index <- calculate_adjusted_rand_index(Clusters, Target)
        DaviesBouldin_index <- calculate_davies_bouldin_index(clusters = Clusters, data = projected_data)
        dbcv_index <- calculate_dbcv(data = projected_data, clusters = Clusters, dist_method = distance_metric)
        CalinskiHarabasz_index <- calculate_calinski_harabasz_index(projected_data, Clusters)
        inertia_value <- calculate_inertia(projected_data, Clusters)
        ami_value <- calculate_ami(Target, Clusters)


        return(data.frame(
          projection_method = projection_method,
          clustering_method = cluster_alg,
          cluster_number_method = cluster_number_method,
          cluster_accuracy = cluster_accuracy,
          Silhouette_index = Silhouette_index,
          Dunn_index = Dunn_index,
          Rand_index = Rand_index,
          DaviesBouldin_index = DaviesBouldin_index,
          dbcv_index = dbcv_index,
          CalinskiHarabasz_index = CalinskiHarabasz_index,
          inertia = inertia_value,
          adjusted_mutual_information = ami_value,
          stringsAsFactors = FALSE
        ))
      })
      do.call(rbind, cluster_scores_per_number_method)
    })
    do.call(rbind, cluster_scores_per_alg)
  }, mc.cores = min(length(projection_methods), nProc))

  cluster_scores_df <- do.call(rbind, cluster_scores)
  rownames(cluster_scores_df) <- apply(cluster_scores_df[, 1:3], 1, function(x) paste0(x, collapse = "_"))

  return(cluster_scores_df)
}


# Rank the metrics and find the best combination of methods
process_cluster_scores <-
  function(cluster_scores_df, selected_cluster_metrics = c("cluster_accuracy", "Silhouette_index", "Dunn_index",
                                                             "Rand_index", "DaviesBouldin_index", "dbcv_index",
                                                             "CalinskiHarabasz_index", "inertia", "adjusted_mutual_information")) {

    # Check if there is something to rank
    if (nrow(cluster_scores_df) > 1) {
      # Validate selected metrics
      valid_cluster_metrics <- c("cluster_accuracy", "Silhouette_index", "Dunn_index", "Rand_index", "DaviesBouldin_index", "dbcv_index",
                                  "CalinskiHarabasz_index", "inertia", "adjusted_mutual_information")

      cluster_metrics_without_orig_classes <- c("Silhouette_index", "Dunn_index", "DaviesBouldin_index",
                                                 "dbcv_index", "inertia", "CalinskiHarabasz_index")
      cluster_metrics_with_orig_classes <- c("cluster_accuracy", "Rand_index", "adjusted_mutual_information")


      selected_cluster_metrics <- intersect(selected_cluster_metrics, valid_cluster_metrics)
      selected_cluster_metrics_without_orig_classes <- intersect(selected_cluster_metrics, cluster_metrics_without_orig_classes)
      selected_cluster_metrics_with_orig_classes <- intersect(selected_cluster_metrics, cluster_metrics_with_orig_classes)

      if (length(selected_cluster_metrics) == 0) {
        message("No valid cluster metrics selected for ranking.\n Reverting to all impleneted metrics: ")
        message(valid_cluster_metrics)
        selected_cluster_metrics <- valid_cluster_metrics
      }

      # Ranking scores: higher is better for most, lower is better for Davies-Bouldin
      rank_columns <- as.data.frame(sapply(valid_cluster_metrics, function(metric) {
        if (metric %in% c("DaviesBouldin_index", "inertia")) {
          -cluster_scores_df[[metric]] # Lower is better, so negate values to rank
        } else {
          cluster_scores_df[[metric]] # Higher is better
        }
      }))

      # Calculate ranks for each metric

      ranked_columns <- as.data.frame(apply(rank_columns, 2, rank, ties.method = "min"))
      colnames(ranked_columns) <- paste0(names(ranked_columns), "_rank")

      # Calculate combined ranks

      if (length(selected_cluster_metrics_without_orig_classes) > 1) {
        ranked_columns$combined_rank_metrics_without_orig_classes <-
          apply(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics_without_orig_classes, "_rank")], 1, prod)
      } else {
        if (length(selected_cluster_metrics_without_orig_classes) > 0) {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics_without_orig_classes, "_rank")]))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank_metrics_without_orig_classes"
        } else {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(rep(1, nrow(ranked_columns))))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank_metrics_without_orig_classes"
        }
      }

      if (length(selected_cluster_metrics_with_orig_classes) > 1) {
        ranked_columns$combined_rank_metrics_with_orig_classes <-
          apply(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics_with_orig_classes, "_rank")], 1, prod)
      } else {
        if (length(selected_cluster_metrics_with_orig_classes) > 0) {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics_with_orig_classes, "_rank")]))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank_metrics_with_orig_classes"
        } else {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(rep(1, nrow(ranked_columns))))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank_metrics_with_orig_classes"
        }
      }

      if (length(selected_cluster_metrics) > 1) {
        ranked_columns$combined_rank <-
          apply(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics, "_rank")], 1, prod)
        # rowMedians( ranked_columns[names( ranked_columns ) %in% paste0( selected_cluster_metrics, "_rank" )] )
      } else {
        if (length(selected_cluster_metrics) > 0) {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(ranked_columns[names(ranked_columns) %in% paste0(selected_cluster_metrics, "_rank")]))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank"
        } else {
          ranked_columns <- cbind.data.frame(ranked_columns,
                                              as.vector(rep(1, nrow(ranked_columns))))
          names(ranked_columns)[ncol(ranked_columns)] <- "dummy"
          names(ranked_columns)[ncol(ranked_columns)] <- "combined_rank"
        }
      }

    } else {

      ranked_columns <- cbind.data.frame(combined_rank_metrics_without_orig_classes = 1,
                                          combined_rank_metrics_with_orig_classes = 1,
                                          combined_rank = 1
      )
    }

    # Add columns for individual ranks to the data frame
    cluster_scores_df <- cbind.data.frame(cluster_scores_df, ranked_columns)

    # Return results as a list
    return(
      ranked_data = cluster_scores_df)
  }



extract_and_safe_cluster_data <- function(cluster_list, selected_columns = c("Label", "Cluster", "Target", "Misclassified")) {
  if (is.null(cluster_list)) stop("Error: `cluster_list` is NULL.")

  # Wrapper to safely perform extraction
  result_list <- list() # To store extracted data
  tryCatch({
    # Traverse the nested list hierarchy
    for (data_set in names(cluster_list)) {
      current_data_set <- cluster_list[[data_set]]

      for (projection_method in names(current_data_set)) {
        current_projection <- current_data_set[[projection_method]]

        for (clustering_method in names(current_projection)) {
          current_clustering <- current_projection[[clustering_method]]

          for (cluster_number_method in names(current_clustering)) {
            # Extract the matrix (data.frame) at the deepest level
            current_matrix <- current_clustering[[cluster_number_method]]

            if (is.data.frame(current_matrix) && all(selected_columns %in% colnames(current_matrix))) {
              # Generate a unique column prefix
              hierarchical_name <- paste(data_set, projection_method, clustering_method, cluster_number_method, sep = "_")

              # Extract required columns as a data.frame
              filtered_data <- current_matrix[, selected_columns, drop = FALSE]

              # Add row number for alignment later (only within the function)
              filtered_data$row_number <- seq_len(nrow(filtered_data))

              # Rename the columns with hierarchical prefixes
              colnames(filtered_data)[colnames(filtered_data) != "row_number"] <-
                paste0(hierarchical_name, "_", selected_columns)

              # Add this data.frame to the result list
              result_list[[hierarchical_name]] <- filtered_data
            }
          }
        }
      }
    }
  }, error = function(e) {
    stop("Error while extracting cluster data: ", conditionMessage(e))
  })

  # Combine all individual data.frames by row_number (using full outer join)
  combined_result <- Reduce(function(x, y) merge(x, y, by = "row_number", all = TRUE), result_list)

  # Remove the temporary row_number column
  combined_result$row_number <- NULL

  return(combined_result)
}

#Function to extract the cluster assignments for each case
extract_cluster_data <- function(cluster_list, selected_columns = c("Label", "Cluster", "Target", "Misclassified")) {
  result_list <- list() # To store extracted data

  # Traverse the nested list hierarchy
  for (data_set in names(cluster_list)) {
    current_data_set <- cluster_list[[data_set]]

    for (projection_method in names(current_data_set)) {
      current_projection <- current_data_set[[projection_method]]

      for (clustering_method in names(current_projection)) {
        current_clustering <- current_projection[[clustering_method]]

        for (cluster_number_method in names(current_clustering)) {
          # Extract the matrix (data.frame) at the deepest level
          current_matrix <- current_clustering[[cluster_number_method]]

          if (is.data.frame(current_matrix) && all(selected_columns %in% colnames(current_matrix))) {
            # Generate a unique column prefix
            hierarchical_name <- paste(data_set, projection_method, clustering_method, cluster_number_method, sep = "_")

            # Extract required columns as a data.frame
            filtered_data <- current_matrix[, selected_columns, drop = FALSE]

            # Add row number for alignment later (only within the function)
            filtered_data$row_number <- seq_len(nrow(filtered_data))

            # Rename the columns with hierarchical prefixes
            colnames(filtered_data)[colnames(filtered_data) != "row_number"] <-
              paste0(hierarchical_name, "_", selected_columns)

            # Add this data.frame to the result list
            result_list[[hierarchical_name]] <- filtered_data
          }
        }
      }
    }
  }

  # Extract the row numbers directly from the 5th column
  row_numbers <- unique(unlist(lapply(result_list, function(df) df[, "row_number"])))

  # Process each extracted data.frame to ensure row alignment
  combined_data <- lapply(result_list, function(df) {
    # Identify missing rows
    missing_rows <- setdiff(row_numbers, df[, "row_number"])
    if (length(missing_rows) > 0) {
      # Create a placeholder for missing rows with NA
      missing_df <- as.data.frame(matrix(NA, ncol = ncol(df), nrow = length(missing_rows)))
      colnames(missing_df) <- colnames(df)
      missing_df[, "row_number"] <- missing_rows

      # Add the missing rows
      df <- rbind(df, missing_df)
    }
    # Order the data.frame by row_number
    return(df[order(df[, "row_number"]),])
  })

  # Combine all data.frames column-wise
  final_data <- do.call(cbind, combined_data)
  names(final_data) <- gsub(".D2", "_D2", names(final_data))

  # Simplify column names: take the part AFTER the dot
  colnames(final_data) <- sapply(colnames(final_data), function(name) {
    if (grepl("\\.", name)) {
      sub(".*\\.", "", name) # Extract only the part after the dot
    } else {
      name # Keep names without dots as is
    }
  })

  # Remove the temporary `row_number` column
  final_data <- final_data[, !colnames(final_data) %in% "row_number", drop = FALSE]

  return(final_data)
}

# Function to create an ABC plot
generate_ABC_plot <- function(misclassified_data, plot_title) {
  p_ABC <- ABCplotGG(misclassified_data) +
    theme_light() +
    theme(
      legend.position = c(0.9, 0.5),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.background = element_rect(color = "transparent", fill = alpha("white", 0.2))
    ) +
    labs(title = plot_title) +
    ggthemes::scale_color_colorblind() + ggthemes::scale_fill_colorblind()
  return(p_ABC)
}


# Function to analyze and plot misclassifications
analyze_and_plot_misclassifications <-
  function(projections_and_plots, palette_ABC, cluster_number_methods_list, dataset_names, alternative_classes) {
    library(ABCanalysis)
    library(ggplot2)

    # Ensure important utilities are defined
    if (missing(palette_ABC)) {
      palette_ABC <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    }

    # Extract cluster assignments and handle "Misclassified" data
    df_cluster_assignments <- safe_extract_cluster_data(cluster_list = projections_and_plots$clustering_results)
    df_cluster_assignments_Misclassified <- df_cluster_assignments[, grep("Misclassified", names(df_cluster_assignments))]
    if (any(grepl("_none", names(df_cluster_assignments_Misclassified)))) {
      df_cluster_assignments_Misclassified <- df_cluster_assignments_Misclassified[, - grep("_none", names(df_cluster_assignments_Misclassified))]
    }
    names(df_cluster_assignments_Misclassified) <- gsub("_orig.*", "", names(df_cluster_assignments_Misclassified))

    # Calculate TimesMisclassified
    df_cluster_assignments_Misclassified_cases <- data.frame(TimesMisclassified = rowSums(df_cluster_assignments_Misclassified, na.rm = TRUE))
    df_cluster_assignments_Misclassified_cases$Case <- rownames(df_cluster_assignments_Misclassified_cases)

    # Ensure `alternative_classes` exists and initialize if not provided
    if (missing(alternative_classes)) {
      alternative_classes <- rep("no", nrow(df_cluster_assignments_Misclassified_cases))
    }

    # Check length of `alternative_classes`
    if (length(alternative_classes) != nrow(df_cluster_assignments_Misclassified_cases)) {
      stop("Error: `alternative_classes` must have the same length as the number of rows in `df_cluster_assignments_Misclassified_cases`.")
    }

    df_cluster_assignments_Misclassified_cases$Suspect <- alternative_classes

    # Define auxiliary functions
    generate_uniform_distribution <- function(values) {
      runif(
        n = length(values),
        min = min(values),
        max = max(values)
      )
    }

    check_uniformity <- function(values, uniform_distribution) {
      p_value <- twosamples::ad_test(values, uniform_distribution)[["P-Value"]]
      return(p_value < 0.05) # Return TRUE if not uniform
    }

    # Perform ABC Analysis

    ABC_plotlist <- list()
    LettersA_for_name <- ""
    ABC_misclassified <- list()
    ABC_misclassified$Aind <- setNames(df_cluster_assignments_Misclassified_cases$TimesMisclassified,
                                        df_cluster_assignments_Misclassified_cases$Case)

    for (i in 1:5) {
      # Generate a uniform distribution for comparison
      TimesMisclassified_uniform <- generate_uniform_distribution(
        values = df_cluster_assignments_Misclassified_cases$TimesMisclassified[
          get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
        ]
      )

      # Check if the data is uniform
      is_not_uniform <- check_uniformity(
        df_cluster_assignments_Misclassified_cases$TimesMisclassified[
          get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
        ],
        TimesMisclassified_uniform
      )

      # Stop the loop if data becomes uniform
      if (!is_not_uniform) break

      # Generate an ABC plot for current iteration
      p_ABC_misclassified <- generate_ABC_plot(
        misclassified_data = setNames(
          df_cluster_assignments_Misclassified_cases$TimesMisclassified[
            get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
          ],
          df_cluster_assignments_Misclassified_cases$Case[
            get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
          ]
        ),
        plot_title = paste0("ABC Analysis #", i)
      ) + scale_color_manual(values = palette_ABC)

      # Perform ABC analysis
      assign(
        paste0(LettersA_for_name, "ABC_misclassified"),
        ABCanalysis::ABCanalysis(
          setNames(
            df_cluster_assignments_Misclassified_cases$TimesMisclassified[
              get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
            ],
            df_cluster_assignments_Misclassified_cases$Case[
              get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
            ]
          ),
          PlotIt = FALSE
        )
      )

      # Store the plot in the list
      ABC_plotlist[[length(ABC_plotlist) + 1]] <- p_ABC_misclassified

      # Recheck uniformity after ABC analysis
      TimesMisclassified_uniform <- generate_uniform_distribution(
        values = df_cluster_assignments_Misclassified_cases$TimesMisclassified[
          get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
        ]
      )
      is_not_uniform <- check_uniformity(
        df_cluster_assignments_Misclassified_cases$TimesMisclassified[
          get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind
        ],
        TimesMisclassified_uniform
      )

      # Stop the loop if data becomes uniform
      if (!is_not_uniform) break

      # Update prefix for the next iteration
      if (i > 1) LettersA_for_name <- paste0(LettersA_for_name, "A")
    }

    if (length(ABC_plotlist) < 1) {
      p_ABC_misclassified_empty <- ggplot() + theme_void() + geom_text(label = "Uniform distribution\nof n missclassified", show.legend = FALSE) +

        ABC_plotlist[[length(ABC_plotlist) + 1]] <- p_ABC_misclassified_empty
    }


    A_set <- df_cluster_assignments_Misclassified_cases$TimesMisclassified[order(df_cluster_assignments_Misclassified_cases$TimesMisclassified, decreasing = TRUE)[1:length(get(paste0(LettersA_for_name, "ABC_misclassified"))$Aind)]]
    df_cluster_assignments_Misclassified_cases$Aset <- ifelse(df_cluster_assignments_Misclassified_cases$TimesMisclassified >= min(A_set), "A", "BC")

    # Plot Misclassified Cases
    p_Misclassified_cases <-
      ggplot(df_cluster_assignments_Misclassified_cases, aes(x = reorder(Case, - TimesMisclassified), y = TimesMisclassified, color = Aset, fill = Suspect)) +
      geom_bar(stat = "identity", alpha = 0.5) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 9),
             legend.position = c(0.7, 0.6),
             legend.direction = "vertical",
             legend.box = "vertical",
             legend.background = element_rect(color = "transparent", fill = ggplot2::alpha("white", 0.2))) +
      scale_color_manual(values = palette_ABC) +
      scale_fill_manual(values = palette_ABC) +
      labs(x = "Case number", title = "Times misclassified per case", color = "ABC set 'A'", fill = "Suspect cases")

    # Combined Plot ABC
    ABC_plots <-
      cowplot::plot_grid(plotlist = ABC_plotlist,
                          nrow = 1)

    # combined_plot_ABC <-
    #   cowplot::plot_grid(p_Misclassified_cases,
    #                      p_ABC_misclassified,
    #                      p_AABC_misclassified,
    #                      labels = LETTERS[1:3],
    #                      nrow = 1, rel_widths = c(2, 1, 1), align = "h", axis = "tb")

    # Return the combined plot
    return(list(
      barplot_Misclassified_cases = p_Misclassified_cases,
      combined_plot_ABCanaylsis = ABC_plots)
    )
  }

# analyze_and_plot_misclassifications <- function(projections_and_plots, dfUmx_Cls, df_Misclassified,
#                                                 final_plot, cb_palette, cluster_number_methods_list, dataset_names) {
#   library(ABCanalysis)
#   library(ggplot2)
#   library(cowplot)
#   library(scales) # Required for `alpha`.
#
#   # Ensure important utilities are defined
#   if (missing(cb_palette)) {
#     cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   }
#
#   safe_extract_cluster_data <- function(cluster_list) {
#     if (is.null(cluster_list)) stop("Error: `cluster_list` is NULL.")
#     extracted_data <- tryCatch(
#       extract_cluster_data(cluster_list),
#       error = function(e) stop("Error while extracting cluster data:", conditionMessage(e))
#     )
#     return(extracted_data)
#   }
#
#   # Extract cluster assignments and handle "Misclassified" data
#   df_cluster_assignments <- safe_extract_cluster_data(cluster_list = projections_and_plots$clustering_results)
#   df_cluster_assignments_Misclassified <- df_cluster_assignments[, grep("Misclassified", names(df_cluster_assignments))]
#   if (any(grepl("_none", names(df_cluster_assignments_Misclassified)))) {
#     df_cluster_assignments_Misclassified <- df_cluster_assignments_Misclassified[, -grep("_none", names(df_cluster_assignments_Misclassified))]
#   }
#   names(df_cluster_assignments_Misclassified) <- gsub("_orig.*", "", names(df_cluster_assignments_Misclassified))
#
#   # Calculate TimesMisclassified
#   df_cluster_assignments_Misclassified_cases <- data.frame(TimesMisclassified = rowSums(df_cluster_assignments_Misclassified, na.rm = TRUE))
#   df_cluster_assignments_Misclassified_cases$Case <- rownames(df_cluster_assignments_Misclassified_cases)
#   df_cluster_assignments_Misclassified_cases$Suspect <-
#     ifelse(df_cluster_assignments_Misclassified_cases$Case %in% rownames(dfUmx_Cls[dfUmx_Cls$ControLikePsA == 1,]), "Yes", "No")
#
#   # Perform ABC Analysis
#   ABC_misclassified <- ABCanalysis::ABCanalysis(
#     setNames(df_cluster_assignments_Misclassified_cases$TimesMisclassified,
#              df_cluster_assignments_Misclassified_cases$Case),
#     PlotIt = TRUE
#   )
#
#   # Generate the first ABC plot using the reusable function
#   p_ABC_misclassified <- generate_ABC_plot(
#     misclassified_data = setNames(df_cluster_assignments_Misclassified_cases$TimesMisclassified,
#                                   df_cluster_assignments_Misclassified_cases$Case),
#     plot_title = "1st ABC analysis"
#   ) + scale_color_manual(values = cb_palette)
#
#   # Recursive ABC Analysis
#   AABC_misclassified <- ABCanalysis::ABCanalysis(
#     setNames(df_cluster_assignments_Misclassified_cases$TimesMisclassified[ABC_misclassified$Aind],
#              names(ABC_misclassified$Aind)),
#     PlotIt = TRUE
#   )
#
#   df_cluster_assignments_Misclassified_cases$Aset <-
#     ifelse(df_cluster_assignments_Misclassified_cases$Case %in% names(AABC_misclassified$Aind), "A", "BC")
#
#   # Generate the second (recursive) ABC plot using the reusable function
#   p_AABC_misclassified <- generate_ABC_plot(
#     misclassified_data = setNames(df_cluster_assignments_Misclassified_cases$TimesMisclassified[ABC_misclassified$Aind],
#                                   names(ABC_misclassified$Aind)),
#     plot_title = "2nd (recursive) ABC analysis"
#   ) + scale_color_manual(values = cb_palette)
#
#   # Plot Misclassified Cases
#   p_Misclassified_cases <-
#     ggplot(df_cluster_assignments_Misclassified_cases, aes(x = reorder(Case, -TimesMisclassified), y = TimesMisclassified, color = Aset, fill = Suspect)) +
#       geom_bar(stat = "identity", alpha = 0.5) +
#       theme_light() +
#       theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 9),
#             legend.position = c(0.7, 0.6),
#             legend.direction = "vertical",
#             legend.box = "vertical",
#             legend.background = element_rect(color = "transparent", fill = alpha("white", 0.2))) +
#       scale_color_manual(values = cb_palette[c(2, 3)]) +
#       scale_fill_manual(values = cb_palette[c(3, 2)]) +
#       labs(x = "Case number", title = "Times misclassified per case", color = "ABC set 'A'", fill = "Controls in PsA patients cluster")
#
#   # Barplot for Misclassification Summary
#   if (nrow(df_Misclassified) == 0) warning("Empty `df_Misclassified`. Barplot might not render correctly.")
#   missclassified_barplot <-
#     ggplot(df_Misclassified, aes(x = reorder(Method, -n_Misclassified), y = n_Misclassified, color = PriorClass, fill = PriorClass)) +
#       geom_bar(stat = "identity", alpha = 0.2) +
#       theme_light() +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#             legend.position = c(0.8, 0.2),
#             legend.direction = "horizontal",
#             legend.box = "horizontal",
#             legend.background = element_rect(color = "transparent", fill = alpha("white", 0.2))) +
#       scale_color_manual(values = cb_palette) +
#       scale_fill_manual(values = cb_palette) +
#       labs(x = "Dataset name_projection method_clustering method", title = "Clustering results related to different prior classifications")
#
#   # Final Combined Plot
#   combined_plot_misclassifieds <-
#     cowplot::plot_grid(cowplot::plot_grid(final_plot, labels = LETTERS[1]),
#                        cowplot::plot_grid(p_Misclassified_cases,
#                                           p_ABC_misclassified,
#                                           p_AABC_misclassified,
#                                           labels = LETTERS[2:4],
#                                           nrow = 1, rel_widths = c(2, 1, 1), align = "h", axis = "tb"),
#                        ncol = 1, rel_heights = c(2, 1))
#
#   # Save the plot
#   ggsave(
#     filename = paste0("Missclassifieds_plot_", cluster_number_methods_list, "_clusters_", dataset_names[1], ".svg"),
#     plot = combined_plot_misclassifieds, width = 26, height = 20, limitsize = FALSE
#   )
#
#   return(combined_plot_misclassifieds)
# }

############### Functions for plotting ###############

# Function to plot the projected data and prior classes or clusters
create_dataset_plots <- function(clustering_results, projection_methods, clustering_methods, cluster_number_methods,
                                  dataset_name, label_points, highlight_misclassified, points_colored_for, cells_colored_for,
                                  palette_target = NULL, palette_cluster = NULL, nProc = 2) {

  projection_plots <- pbmcapply::pbmclapply(projection_methods, function(projection_method) {

    cluster_algs_plots <- lapply(clustering_methods, function(cluster_alg) {
      cluster_number_method_plots <- lapply(cluster_number_methods, function(cluster_number_method) {

        combined_result <- clustering_results[[projection_method]][[cluster_alg]][[cluster_number_method]]

        if (is.null(combined_result)) {
          dfPlot <- data.frame(x = 1, y = 1)
          text_to_show <- "Projection\nfailed"
          plot <- ggplot(data = dfPlot, aes(x = x, y = y, color = x, fill = y)) +
            geom_text(label = text_to_show, show.legend = FALSE) +
            theme_light() +
            labs(x = "Dim 1", y = "Dim 2")
        } else {

          # Initialize variables
          max_attempts <- 10
          plot <- NULL

          # Loop for attempting to plot with increasing jitter
          for (attempt in 1:max_attempts) {
            plot1 <- tryCatch({
              plotVoronoiTargetProjection(
                X = combined_result[, 1:2],
                targets = combined_result$Target,
                clusters = combined_result$Cluster,
                labels = combined_result$Label,
                misclassified = combined_result$Misclassified,
                LabelPoints = label_points,
                ColorMisClassified = highlight_misclassified,
                points_colored_for = points_colored_for,
                cells_colored_for = cells_colored_for,
                palette_target = palette_target,
                palette_cluster = palette_cluster,
                cluster_alg = cluster_alg
              )
            }, error = function(e) {
              return(structure(list(message = e$message), class = "try-error"))
            })

            # Check if the plot was successful
            if ((!inherits(plot1, "try-error") && !is.null(plot1)) || inherits(plot1, "gg")) {
              plot <- plot1
              break # Exit the loop if successful
            } else {
              # Prepare columns for jittering
              columns_to_jitter <- intersect(setdiff(names(combined_result),
                                                       c("Cluster", "Target", "Label", "Misclassified")),
                                              names(combined_result))
              numeric_columns_to_jitter <- sapply(combined_result[, columns_to_jitter], is.numeric)

              # Apply jitter with increasing factor based on the attempt number
              jitter_factor <- attempt * 10 # Increase jitter factor with each attempt
              combined_result[, columns_to_jitter[numeric_columns_to_jitter]] <-
                lapply(combined_result[, columns_to_jitter[numeric_columns_to_jitter], drop = FALSE], function(x) jitter(x, factor = jitter_factor))
            }
          }

          # If all attempts failed, create a fallback plot
          if (is.null(plot)) {
            dfPlot <- data.frame(x = 1, y = 1)
            text_to_show <- "Projection\nfailed"
            plot <- ggplot(data = dfPlot, aes(x = x, y = y, color = x, fill = y)) +
              geom_text(label = text_to_show, show.legend = FALSE) +
              theme_light() +
              labs(x = "Dim 1", y = "Dim 2")
          }
        }

        plot_title <- paste(dataset_name, ": ", projection_method)
        if (cluster_alg != "none") {
          plot_title <- paste(plot_title, "- ", cluster_alg, "- ", cluster_number_method)
        }
        plot <- plot + labs(title = plot_title)
        if (cluster_alg == "none")
          plot <- plot + labs(fill = "Target")

        return(plot) # Ensuring plot is returned
      })
      names(cluster_number_method_plots) <- cluster_number_methods
      return(cluster_number_method_plots)
    })
    names(cluster_algs_plots) <- clustering_methods
    return(cluster_algs_plots)
  }, mc.cores = min(length(projection_methods), nProc))
  names(projection_plots) <- projection_methods
  return(projection_plots)
}


# Function to plot Voronoi cells with projection and clustering results
plotVoronoiTargetProjection <- function(X, targets, clusters = NULL,
                                         labels = NULL, LabelPoints = FALSE,
                                         misclassified = NULL, ColorMisClassified = FALSE,
                                         points_colored_for = "Target", cells_colored_for = "Cluster",
                                         palette_target = NULL, palette_cluster = NULL,
                                         cluster_alg = NULL) {
  # Extended colorblind palette
  cb_palette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  shape_values <- c(17, 16, 15, 18, 2, 1, 0, 5, 7, 8, 9, 19, 11, 12)

  # Initial checks of arguments
  if (ncol(X) < 2) stop("The input data `X` must have at least two columns.")
  if (nrow(X) != length(targets)) stop("The length of `targets` must match the number of rows in `X`.")

  if (!is.null(clusters) && length(clusters) != nrow(X))
    stop("The length of `clusters` must match the number of rows in `X`.")
  if (!is.null(misclassified) && length(misclassified) != length(labels))
    stop("The length of `misclassified` must match the number of labels.")

  if (!is.null(palette_target) && length(palette_target) < length(unique(targets))) {
    stop("The `palette_target` must provide enough color values for the number of target groups.")
  }
  if (!is.null(palette_cluster) && length(palette_cluster) < length(unique(clusters))) {
    stop("The `palette_cluster` must provide enough color values for the number of clusters.")
  }

  # Define default palettes if none are provided
  if (is.null(palette_target)) palette_target <- rep(cb_palette, length.out = length(unique(targets)))
  if (is.null(palette_cluster)) palette_cluster <- rep(cb_palette, length.out = length(unique(clusters)))

  # Assign labels if not present
  if (is.null(labels)) {
    if (!is.null(rownames(X))) {
      labels <- rownames(X)
    } else {
      labels <- seq_len(nrow(X))
    }
  }

  # Default values for clusters and misclassified if NULL
  if (is.null(clusters)) clusters <- targets
  if (is.null(misclassified)) misclassified <- rep(0, nrow(X))

  # Data preparation
  plotData <- data.frame(Proj1 = X[, 1], Proj2 = X[, 2], Target = targets,
                          Clusters = clusters, Label = labels, Misclassified = misclassified)

  plotData$DotsInformation <- if (points_colored_for == "Target") plotData$Target else plotData$Clusters
  plotData$CellsInformation <- if (cells_colored_for == "Cluster") plotData$Clusters else plotData$Target

  # Voronoi diagram computation
  voronoiData <- deldir::deldir(plotData$Proj1, plotData$Proj2)

  # Convert Voronoi tessellation to a data frame for plotting
  vor_polys <- deldir::tile.list(voronoiData)
  voronoi_df <- do.call(rbind, lapply(seq_along(vor_polys), function(i) {
    data.frame(x = vor_polys[[i]]$x, y = vor_polys[[i]]$y, id = i)
  }))

  # Create plot with ggplot2
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = voronoi_df, ggplot2::aes(x = x, y = y, group = id, fill = as.factor(plotData$CellsInformation[id])),
                           alpha = 0.3, color = NA) +
    ggplot2::geom_point(data = plotData, ggplot2::aes(x = Proj1, y = Proj2,
                                                        color = as.factor(DotsInformation), shape = as.factor(DotsInformation))) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.5, 0.08),
                    legend.box = "horizontal", legend.direction = "horizontal",
                    legend.background = element_rect(color = "transparent", fill = ggplot2::alpha("white", 0.2))) +
  #ggplot2::labs( x = "Dim 1", y = "Dim 2", color = "Target", fill = "Cluster", shape = "Target" ) +
  scale_shape_manual(values = shape_values) +
    ggplot2::scale_fill_manual(values = palette_cluster) +
    ggplot2::scale_color_manual(values = palette_target) +
    guides(shape = guide_legend(override.aes = list(label = "")))

  if (cluster_alg == "none") {
    plot <- plot + ggplot2::labs(x = "Dim 1", y = "Dim 2", color = "Target", fill = "Target", shape = "Target")
  } else {
    plot <- plot + ggplot2::labs(x = "Dim 1", y = "Dim 2", color = "Target", fill = "Cluster", shape = "Target")
  }
  # Conditional labeling of points
  if (LabelPoints) {
    plot <- plot +
      ggrepel::geom_text_repel(data = plotData,
                                ggplot2::aes(x = Proj1, y = Proj2, color = as.factor(DotsInformation),
                                              label = Label, fontface = 2), vjust = -1, size = 3, max.overlaps = Inf)
  }

  # Color misclassified cases in red
  if (sum(plotData$Misclassified) > 0 & ColorMisClassified) {
    suppressMessages({
      plot <- plot +
        geom_text(
          aes(x = Inf, y = Inf,
               label = paste0(round(100 * sum(plotData$Misclassified) / nrow(X), 1), "% misclassified")),
          hjust = 1,
          vjust = 1
        )
    })
    q <- ggplot_build(plot)
    q$data[[2]]$colour[plotData$Misclassified == 1] <- "red"
    gtable <- ggplot2::ggplot_gtable(q)
    plot <- ggplotify::as.ggplot(function() grid::grid.draw(gtable))
  }

  return(plot)
}


# Function to extract all plots from the generated lists and combine them into a figure
combine_all_plots <- function(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods, nProc = 2) {
  if (any(sapply(list(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods), is.null))) {
    stop("All input lists must be non-null and contain elements.")
  }

  if (any(sapply(list(datasets, projection_plots, projection_methods, clustering_methods, cluster_number_methods), length) == 0)) {
    stop("All input lists must contain non-empty elements.")
  }

  calculate_plots_per_dataset <- function()
    length(projection_methods) *
    length(clustering_methods) *
    length(cluster_number_methods)

  plots_per_dataset <- calculate_plots_per_dataset()

  if (plots_per_dataset == 0) {
    stop("No valid combinations found in input lists.")
  }

  all_plots <- unlist(lapply(projection_plots, function(projection) {
    if (!is.null(projection)) {
      return(unlist(projection, recursive = FALSE, use.names = FALSE))
    }
  }), recursive = FALSE, use.names = FALSE)

  if (length(all_plots) == 0) {
    stop("No plots available to combine.")
  }

  figure_count <- length(all_plots) %/% plots_per_dataset

  if (figure_count == 0) {
    stop("Insufficient plots to create figures based on input methods.")
  }

  if (length(clustering_methods) == 1 && clustering_methods == "none") {
    columns <- calculate_golden_matrix_dims(n_items = length(projection_methods))$ncol
    rows <- calculate_golden_matrix_dims(n_items = length(projection_methods))$nrow
  } else {
    columns <- length(clustering_methods) * length(cluster_number_methods)
    rows <- length(projection_methods)
  }

  create_combined_plot <- function(idx) {
    start_idx <- 1 + (idx - 1) * plots_per_dataset
    end_idx <- plots_per_dataset + (idx - 1) * plots_per_dataset

    if (end_idx > length(all_plots)) {
      stop(sprintf("Index out of range when creating plot: %s", idx))
    }

    if (length(clustering_methods) == 1 && clustering_methods == "none") {
      fig <- cowplot::plot_grid(plotlist = all_plots[start_idx:end_idx], ncol = columns, nrow = rows, labels = "AUTO")
    } else {
      fig <- cowplot::plot_grid(plotlist = all_plots[start_idx:end_idx], ncol = columns, nrow = rows,
                                 labels =
                                  paste0(rep(LETTERS[1:rows], each = columns), rep(1:columns, rows))
                                )
    }
    return(fig)
  }


  if (figure_count > 1) {
    figures <- lapply(
      seq_along(datasets),
      create_combined_plot
    )
  } else {
    figures <- lapply(seq_along(datasets), create_combined_plot)
  }

  names(figures) <- datasets
  return(figures)
}


# Function to calculate a compact matrix for a combined results figure
calculate_compact_matrix_dims <- function(n_items) {
  if (n_items < 1) {
    stop("Number of items must be at least 1.")
  }

  # Calculate the number of rows (ceiling of square root)
  nrow <- ceiling(sqrt(n_items))

  # Calculate the number of columns (ceiling of items divided by rows)
  ncol <- ceiling(n_items / nrow)

  return(list(nrow = nrow, ncol = ncol))
}


calculate_full_matrix_dims <- function(n_items, vertical = TRUE) {
  # Find all possible factors of n
  factors <- which(n_items %% 1:n_items == 0)

  # Find the pair of factors that minimizes the difference between rows and columns
  best_pair <- which.min(abs(factors - n_items / factors))


  # Determine nrow and ncol
  if (vertical) {
    ncol <- min(factors[best_pair], n_items / factors[best_pair])
    nrow <- max(factors[best_pair], n_items / factors[best_pair])
  } else {
    nrow <- min(factors[best_pair], n_items / factors[best_pair])
    ncol <- max(factors[best_pair], n_items / factors[best_pair])
  }
  # Return the result as a named list
  return(list(nrow = nrow, ncol = ncol))
}

calculate_golden_matrix_dims <- function(n_items, vertical = TRUE, golden_ratio = 1.618, max_ratio = 4) {
  if (n_items <= 0) {
    stop("n_items must be a positive integer")
  }

  # Find all possible factors of n_items
  factors <- which(n_items %% 1:n_items == 0)

  # Calculate potential ratios for each pair of factors
  ratios <- factors / (n_items / factors)
  inverted_ratios <- (n_items / factors) / factors

  # Filter out invalid ratios
  valid_indices <- which(ratios <= max_ratio & inverted_ratios <= max_ratio)

  if (length(valid_indices) == 0 || max(factors[valid_indices]) == n_items) {
    # Revert to compact matrix calculation if no valid indices or single row
    nrow <- ceiling(sqrt(n_items))
    ncol <- ceiling(n_items / nrow)

    if (!vertical) {
      temp <- nrow
      nrow <- ncol
      ncol <- temp
    }
  } else {
    # Find the pair of factors closest to the golden ratio
    best_index <- valid_indices[which.min(abs(ratios[valid_indices] - golden_ratio))]
    best_factor <- factors[best_index]

    if (vertical) {
      ncol <- min(best_factor, n_items / best_factor)
      nrow <- max(best_factor, n_items / best_factor)
    } else {
      nrow <- min(best_factor, n_items / best_factor)
      ncol <- max(best_factor, n_items / best_factor)
    }
  }

  # Calculate how full the matrix is
  fullness <- n_items / (nrow * ncol)

  # Return the result as a named list
  return(list(nrow = nrow, ncol = ncol, fullness = fullness, ratio = max(nrow / ncol, ncol / nrow)))
}


# Update plot titles if the are too long
update_plot_titles <- function(projection_plots, title_size = 7, wrap_chars = 25) {
  lapply(projection_plots, function(proj_plots) {
    lapply(proj_plots, function(cluster_plots) {
      lapply(cluster_plots, function(single_plot) {
        if (!is.null(single_plot) && inherits(single_plot, "ggplot")) {
          # Get current title and wrap it
          current_title <- single_plot$labels$title
          if (!is.null(current_title)) {
            wrapped_title <- str_wrap(current_title, width = wrap_chars)
            single_plot <- single_plot + labs(title = wrapped_title)
          }
          single_plot + theme(
            plot.title = element_text(
              size = title_size,
              face = "plain",
              lineheight = 0.8,    # Controls spacing between wrapped lines
              margin = margin(b = 4)
            )
          )
        } else {
          single_plot
        }
      })
    })
  })
}
