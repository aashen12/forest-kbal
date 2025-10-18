library(randomForest)
library(tidyverse)
library(irlba)

rf_kernel_matrix <- function(model, data, X, n_components, verbose = FALSE) {
  # model: a fitted randomForest object
  # data: data frame or matrix of predictors (no response column)
  
  # extract terminal‐node IDs: an N × n_trees matrix
  preds <- predict(model, newdata = as.data.frame(data), nodes = TRUE, proximity = TRUE)
  terminal_nodes <- attributes(preds)$nodes
  if (verbose) print("RF Terminal nodes extracted")
  
  K1 <- preds$proximity
  
  N <- nrow(terminal_nodes)
  n_trees <- ncol(terminal_nodes)

  # build the kernel by counting leaf‐agreements
  K <- matrix(0, N, N)
  for (tree_idx in seq_len(n_trees)) {
    K <- K + (outer(terminal_nodes[, tree_idx], terminal_nodes[, tree_idx], "=="))
  }

  # normalize
  K <- K / n_trees
  
  check <- K1 == K
  
  # print output of check to tell user if matrices are equal
  if (all(check)) {
    if (verbose) print("GOOD!. Proximity matrix matches leaf-agreement kernel matrix.")
  } else {
    if (verbose) print("Warning: Proximity matrix does not match leaf-agreement kernel matrix.")
  }
  
  # preserve names if present
  if (!is.null(rownames(data))) {
    rownames(K) <- colnames(K) <- rownames(data)
  }
  
  if (verbose) print(paste("Finished creating forest kernel matrix. Running SVD."))
  
  svd_result <- irlba(K, nv = n_components, maxit = 2000, verbose = F)
  if (verbose) print("First RF PCA finished")
  # Scale by square root of singular values
  features <- svd_result$u %*% diag(sqrt(svd_result$d)) %>% data.frame()
  names(features) <- paste0("PC", 1:n_components)
  data_rf <- cbind(data, features)
  
  # Xs <- scale(X)
  # K_mixed <- K + Xs %*% t(Xs)
  # svd_result_mixed <- irlba(K_mixed, nv = n_components, maxit = 2000, verbose = F)
  # if (verbose) print("Second RF PCA finished")
  # # Scale by square root of singular values
  # features_mixed <- svd_result_mixed$u %*% diag(svd_result_mixed$d) %>% data.frame()
  # names(features_mixed) <- paste0("PC", 1:n_components)
  
  # Find elbo point
  # 1. raw eigenvalues
  eigs_raw <- svd_result$d^2
  n <- nrow(K)
  
  # 2. helper to find elbow index
  find_elbow <- function(y) {
    n <- length(y)
    x <- seq_len(n)
    
    # endpoints
    x1 <- x[1]; y1 <- y[1]
    x2 <- x[n]; y2 <- y[n]
    
    # distances from each point to the line (x1,y1)-(x2,y2)
    num <- abs((y2 - y1) * (x - x1) - (x2 - x1) * (y - y1))
    denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    dists <- num / denom
    
    which.max(dists)
  }
  
  elbow <- find_elbow(y = eigs_raw)
  
  # 3. build data.frame for ggplot
  pc_df <- data.frame(
    PC = seq_along(eigs_raw),
    Eigenvalue = eigs_raw
  )
  
  # 4. plot
  library(ggplot2)
  gg <- ggplot(pc_df, aes(PC, Eigenvalue)) +
    geom_line() +
    geom_point() +
    # highlight elbow
    geom_point(
      data = pc_df[elbow, ],
      aes(PC, Eigenvalue),
      color = "red", size = 3
    ) +
    geom_vline(xintercept = elbow, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = elbow, y = pc_df$Eigenvalue[elbow],
      label = paste0("Elbow = ", elbow),
      vjust = -1, color = "red"
    ) +
    labs(
      title = "Scree Plot with Elbow",
      x = "Principal Component",
      y = "Eigenvalue"
    ) +
    theme_minimal()
  
  K_df <- as.data.frame(K)
  names(K_df) <- paste0("PC", 1:ncol(K))
  
  if (verbose) print("Finished RF Kernel")
  
  return(list(
    # leaf_encoding = leaf_encoding,
    data_rf = data_rf,
    features = features,
    #features_mixed = features_mixed,
    explained_variance = svd_result$d^2 / sum(svd_result$d^2),
    scree = gg,
    elbow = elbow,
    K = K_df
  ))
}

randomForestRegression <- function(data, true_att, feat_rep, ncomp = 5, verbose = FALSE) {
  #print(paste("True ATT:", round(true_att, 3)))
  data <- data %>% 
    dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  covs <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(covs), data))
  trt <- data$treat
  n <- nrow(data)
  
  data0 <- data %>% dplyr::filter(treat == 0)
  data1 <- data %>% dplyr::filter(treat == 1)
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data0, ntree = 100)
  rfpreds_1 <- predict(rfmod, newdata = data1, type = "response")
  
  # Next, extract feat reps
  rf_obj <- extract_rf_features(rfmod, data1, n_components = ncomp)
  rf_feats <- rf_obj$features # U matrix
  data_rf_feats <- data.frame(Y = data1$Y, rf_feats)
  #return(data_rf_feats)
  lmmod <- lm(Y ~ ., data = data_rf_feats)
  lm_rf_preds <- lmmod$fitted.values #predict(lmmod, newdata = data1, type = "response")
  
  att.raw.rf <- mean(data1$Y) - mean(rfpreds_1)
  att.ols.rf <- mean(data1$Y) - mean(lm_rf_preds)
  
  bias.raw.rf <- true_att - att.raw.rf
  bias.ols.rf <- true_att - att.ols.rf
  
  df_raw_rf <- data.frame(est = "raw.rf", "dgp" = dgp, "feat_rep" = feat_rep, 
                          "bias" = bias.raw.rf, "cvg" = NA, pbr = NA, ess = NA,
                          "true.att" = true_att, "est.att" = att.raw.rf)
  
  df_ols_rf <- data.frame(est = "ols.rf", "dgp" = dgp, "feat_rep" = feat_rep, 
                          "bias" = bias.ols.rf, "cvg" = NA, pbr = NA, ess = NA,
                          "true.att" = true_att, "est.att" = att.ols.rf)
  data.frame(bind_rows(df_raw_rf, df_ols_rf))
}

