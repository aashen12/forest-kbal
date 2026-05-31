# =============================================================================
# Random Forest Kernel Matrix Extraction and PCA
# =============================================================================
#
# Extracts the implied kernel matrix K^{rf} from a fitted random forest model
# based on leaf-node co-occurrence: K^{rf}_{ij} = (1/T) sum_t I(i,j share leaf).
# See Section 3.1 of Shen et al. (2025).
#
# Functions:
#   rf_kernel_matrix()  - Build RF kernel and extract PCA features
#
# Depends on: sim-utils.R (find_elbow)
# =============================================================================


#' Build RF kernel matrix and extract PCA features.
#'
#' Given a fitted randomForest model, computes the kernel matrix where
#' K[i,j] = proportion of trees in which i and j land in the same leaf.
#' Then applies truncated SVD to extract principal components.
#'
#' @param model A fitted randomForest object (trained on pilot sample)
#' @param data Analysis sample data frame (with Z, Y, and covariates)
#' @param X Covariate matrix for the analysis sample
#' @param n_components Number of principal components to extract
#' @param verbose Print progress messages
#' @return List with data_rf, features, explained_variance, scree plot, elbow index, K data frame
rf_kernel_matrix <- function(model, data, X, n_components, verbose = FALSE) {
  # Extract terminal-node IDs: N x n_trees matrix
  preds <- predict(model, newdata = as.data.frame(data), nodes = TRUE, proximity = TRUE)
  terminal_nodes <- attributes(preds)$nodes
  if (verbose) print("RF terminal nodes extracted")

  N <- nrow(terminal_nodes)
  n_trees <- ncol(terminal_nodes)

  # Build kernel by counting leaf-agreements
  K <- matrix(0, N, N)
  for (tree_idx in seq_len(n_trees)) {
    K <- K + (outer(terminal_nodes[, tree_idx], terminal_nodes[, tree_idx], "=="))
  }
  K <- K / n_trees

  if (!is.null(rownames(data))) {
    rownames(K) <- colnames(K) <- rownames(data)
  }

  if (verbose) print("Forest kernel matrix created. Running SVD.")

  # Truncated SVD for PCA
  svd_result <- irlba(K, nv = n_components, maxit = 2000, verbose = FALSE)
  if (verbose) print("RF PCA finished")

  # Scale PCs by sqrt of singular values (Section 2.3)
  features <- svd_result$u %*% diag(sqrt(svd_result$d)) %>% data.frame()
  names(features) <- paste0("PC", 1:n_components)
  data_rf <- cbind(data, features)

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
    ggplot2::labs(title = "RF Scree Plot", x = "Principal Component", y = "Eigenvalue") +
    ggplot2::theme_minimal()

  K_df <- as.data.frame(K)
  names(K_df) <- paste0("PC", 1:ncol(K))

  if (verbose) print("Finished RF kernel")

  list(data_rf = data_rf, features = features,
       explained_variance = svd_result$d^2 / sum(svd_result$d^2),
       scree = gg, elbow = elbow, K = K_df)
}
