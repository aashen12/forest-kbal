library(Matrix)
library(stochtree)
library(tidyverse)
library(rsample)
library(irlba)
library(kernlab)

bart_kernel_matrix <- function(train, test, seed = 1022, verbose = FALSE, simulation = TRUE) {
  train_con <- train %>% filter(Z == 0) 
  test_con <- test %>% filter(Z == 0)
  test_trt <- test %>% filter(Z == 1)
  target = test

  if (simulation) {
    X_train  <- train_con %>% dplyr::select(starts_with('X'))
    #X_source <- test_con  %>% dplyr::select(starts_with('X'))
    X_target <- target %>% dplyr::select(starts_with('X'))
  } else {
    covs <- c(
      "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
      "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
      "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
      "mathPre"
    )
    X_train  <- train_con %>% dplyr::select(all_of(covs))
    #X_source <- test_con  %>% dplyr::select(all_of(covs))
    X_target <- target %>% dplyr::select(all_of(covs))
  }
  
  y_train = train_con$Y
  n_train = length(y_train)
  
  # Run BART on the training control sample
  num_trees <- 100
  sigma_leaf <- 1/num_trees
  
  mean_forest_params <- list(num_trees=num_trees, sigma2_leaf_init=sigma_leaf)
  bart_model <- bart(X_train=X_train, y_train=y_train, 
                     mean_forest_params = mean_forest_params,
                     num_burnin=200, num_mcmc=1000, 
                     general_params=list(random_seed=seed))
  
  if (verbose) print("BART model trained")
  
  nsamp = bart_model$model_params$num_samples
  nkeep = 50 # Thinning MCMC output to save time

  n_test = nrow(test)
  kernel_post = matrix(0, nrow = n_test, ncol=n_test)
  kernel = kernel_post
  
  muhat = matrix(0, nrow=n_test, ncol=nkeep)
  
  ix = round(seq(1, bart_model$model_params$num_samples, length.out=nkeep))
  i = 1
  for(ii in ix) {
    if (verbose & (i %% 5 == 0)) print(paste("On sample", i, "of", nkeep))
    
    leaf_mat_train  <- computeForestLeafIndices(bart_model, X_train, 
                                                forest_type = "mean", 
                                                forest_inds = ii - 1)
    leaf_mat_target <- computeForestLeafIndices(bart_model, X_target, 
                                                forest_type = "mean", 
                                                forest_inds = ii - 1)
    
    n_features = bart_model$mean_forests$num_forest_leaves(ii-1)
    
    var_leaf = bart_model$sigma2_leaf_samples[ii]
    sigma2   = bart_model$sigma2_global_samples[ii]
    W_train   <- sparseMatrix(i=rep(1:nrow(X_train),num_trees), 
                              j=leaf_mat_train + 1, x=1,
                              dims=c(nrow(X_train), n_features))*sqrt(var_leaf)

    W_target  <- sparseMatrix(i=rep(1:nrow(X_target),num_trees), 
                              j=leaf_mat_target + 1, x=1, 
                              dims=c(nrow(X_target), n_features))*sqrt(var_leaf)
    W = W_target

    cur_kernel = W%*%t(W)
    kernel      = kernel + cur_kernel/nkeep
    
    A = W %*% t(W_train)
    B = solve(W_train %*% t(W_train) + diag(sigma2, nrow(W_train)), t(A))
    muhat[,i] = matrix(A%*%solve(W_train %*% t(W_train) + diag(sigma2, nrow(W_train)), y_train))
    kernel_post = kernel_post + (cur_kernel -  A%*% B)/nkeep

    i = i+1
  }
  if (verbose) print("Finished all samples")
  kernel_post_f = kernel_post + cov(t(muhat))
  return(list(kernel = kernel, kernel_post_f = kernel_post_f, kernel_post=kernel_post))
}

pca_bart <- function(kernel, data, n_components = 10) {
  K <- as.matrix(kernel)
  svd_result <- irlba(K, nv = n_components, maxit = 2000, verbose = F)
  
  # Scale by square root of singular values
  features <- svd_result$u %*% diag(svd_result$d) %>% data.frame()
  names(features) <- paste0("PC", 1:n_components)
  data_bart <- cbind(data, features)
  
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
  
  elbow <- find_elbow(eigs_raw)
  
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
  
  K_df <- as.data.frame(as.matrix(K))
  names(K_df) <- paste0("PC", 1:ncol(K))
  
  return(list(
    # leaf_encoding = leaf_encoding,
    data_bart = data_bart,
    features = features,
    explained_variance = svd_result$d^2 / sum(svd_result$d^2),
    scree = gg,
    elbow = elbow,
    K = K_df
  ))
}





