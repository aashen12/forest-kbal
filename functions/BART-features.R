# =============================================================================
# BART Kernel Matrix Extraction and PCA
# =============================================================================
#
# Trains a BART model on a pilot sample of controls, then extracts the
# implied kernel matrix K^{bart} for the analysis sample based on leaf-node
# co-occurrences across posterior draws. See Section 3.1 of Shen et al. (2025).
#
# Functions:
#   bart_kernel_matrix()  - Train BART and compute kernel matrix
#   pca_bart()            - PCA on the BART kernel for feature extraction
#
# Depends on: sim-utils.R (find_elbow)
# =============================================================================


#' Train BART on pilot controls and extract kernel matrix for analysis sample.
#'
#' The kernel is K^{bart}_{ij} = (1/BT) sum_b sum_t I(i,j share leaf in tree t, draw b).
#' We thin to nkeep=50 posterior draws to reduce computation.
#'
#' @param train Pilot sample (controls only, used for training BART)
#' @param test Analysis sample (treated + control, kernel computed here)
#' @param seed Random seed for BART
#' @param verbose Print progress messages
#' @param dataset Dataset identifier (affects covariate selection)
#' @param covs Character vector of covariate names
#' @return List with kernel (expected), kernel_post, kernel_post_f (posterior with prediction variance)
bart_kernel_matrix <- function(train, test, seed = 1022, verbose = FALSE, dataset = "soldiering", covs) {
  train_con <- train %>% dplyr::filter(Z == 0)
  target <- test

  # Select covariates based on dataset type
  if (dataset == "simulation") {
    X_train  <- train_con %>% dplyr::select(starts_with("X"))
    X_target <- target %>% dplyr::select(starts_with("X"))
  } else {
    X_train  <- train_con %>% dplyr::select(all_of(covs))
    X_target <- target %>% dplyr::select(all_of(covs))
  }

  y_train <- train_con$Y
  X_train  <- as.matrix(X_train)
  X_target <- as.matrix(X_target)

  # BART model: 100 trees, sigma_leaf = 1/100, 200 burnin, 1000 MCMC draws
  num_trees <- 100
  sigma_leaf <- 1 / num_trees
  mean_forest_params <- list(num_trees = num_trees, sigma2_leaf_init = sigma_leaf)

  bart_model <- bart(X_train = X_train, y_train = y_train,
                     mean_forest_params = mean_forest_params,
                     num_burnin = 200, num_mcmc = 1000,
                     general_params = list(random_seed = seed))

  if (verbose) print("BART model trained")

  nsamp <- bart_model$model_params$num_samples
  nkeep <- 50  # Thin MCMC output to save time

  n_test <- nrow(test)
  kernel_post <- matrix(0, nrow = n_test, ncol = n_test)
  kernel <- kernel_post
  muhat <- matrix(0, nrow = n_test, ncol = nkeep)

  # Iterate over thinned posterior draws
  ix <- round(seq(1, nsamp, length.out = nkeep))
  i <- 1
  for (ii in ix) {
    if (verbose & (i %% 5 == 0)) print(paste("On sample", i, "of", nkeep))

    leaf_mat_train  <- computeForestLeafIndices(bart_model, X_train,
                                                forest_type = "mean",
                                                forest_inds = ii - 1)
    leaf_mat_target <- computeForestLeafIndices(bart_model, X_target,
                                                forest_type = "mean",
                                                forest_inds = ii - 1)

    n_features <- bart_model$mean_forests$num_forest_leaves(ii - 1)
    var_leaf <- bart_model$sigma2_leaf_samples[ii]
    sigma2   <- bart_model$sigma2_global_samples[ii]

    # Sparse leaf-encoding matrices (scaled by sqrt of leaf variance)
    W_train  <- sparseMatrix(i = rep(1:nrow(X_train), num_trees),
                             j = leaf_mat_train + 1, x = 1,
                             dims = c(nrow(X_train), n_features)) * sqrt(var_leaf)
    W_target <- sparseMatrix(i = rep(1:nrow(X_target), num_trees),
                             j = leaf_mat_target + 1, x = 1,
                             dims = c(nrow(X_target), n_features)) * sqrt(var_leaf)
    W <- W_target

    cur_kernel <- W %*% t(W)
    kernel <- kernel + cur_kernel / nkeep

    A <- W %*% t(W_train)
    B <- solve(W_train %*% t(W_train) + diag(sigma2, nrow(W_train)), t(A))
    muhat[, i] <- matrix(A %*% solve(W_train %*% t(W_train) + diag(sigma2, nrow(W_train)), y_train))
    kernel_post <- kernel_post + (cur_kernel - A %*% B) / nkeep

    i <- i + 1
  }

  if (verbose) print("Finished all samples")
  kernel_post_f <- kernel_post + cov(t(muhat))
  if (verbose) print("Finished all BART kernels")

  list(kernel = kernel, kernel_post_f = kernel_post_f, kernel_post = kernel_post)
}


#' PCA on the BART kernel matrix for feature extraction.
#'
#' Applies truncated SVD (via irlba) to the BART kernel matrix and returns
#' principal components scaled by sqrt of singular values: sigma_k * U_k.
#'
#' @param kernel BART kernel matrix (from bart_kernel_matrix)
#' @param data Analysis sample data frame
#' @param X Covariate matrix (unused, kept for interface consistency)
#' @param n_components Number of principal components to extract
#' @param verbose Print progress
#' @return List with data_bart, features, explained_variance, scree plot, elbow index, K data frame
pca_bart <- function(kernel, data, X, n_components = 10, verbose = FALSE) {
  if (verbose) print("Starting BART PCA")
  K <- as.matrix(kernel)
  svd_result <- irlba(K, nv = n_components, maxit = 2000, verbose = FALSE)
  if (verbose) print("BART PCA finished")

  # Scale PCs by sqrt of singular values (Section 2.3)
  features <- svd_result$u %*% diag(sqrt(svd_result$d)) %>% data.frame()
  names(features) <- paste0("PC", 1:n_components)
  data_bart <- cbind(data, features)

  # Scree plot and elbow detection
  eigs_raw <- svd_result$d^2
  elbow <- find_elbow(eigs_raw)

  pc_df <- data.frame(PC = seq_along(eigs_raw), Eigenvalue = eigs_raw)
  gg <- ggplot2::ggplot(pc_df, ggplot2::aes(PC, Eigenvalue)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_point(data = pc_df[elbow, ], color = "red", size = 3) +
    ggplot2::geom_vline(xintercept = elbow, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = elbow, y = pc_df$Eigenvalue[elbow],
                      label = paste0("Elbow = ", elbow), vjust = -1, color = "red") +
    ggplot2::labs(title = "BART Scree Plot", x = "Principal Component", y = "Eigenvalue") +
    ggplot2::theme_minimal()

  K_df <- as.data.frame(as.matrix(K))
  names(K_df) <- paste0("PC", 1:ncol(K))

  list(data_bart = data_bart, features = features,
       explained_variance = svd_result$d^2 / sum(svd_result$d^2),
       scree = gg, elbow = elbow, K = K_df)
}
