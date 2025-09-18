
bal_weights_varsweep <- function(data, true_att, feat_rep, raw_covs, kernel_covs, 
                                 varsweep = 1, verbose = FALSE) {
  # l2 or inf
  #print(paste("True ATT:", round(true_att, 3)))
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  if (is.null(raw_covs)) {
    # scenario where we use kernel only
    covs <- kernel_covs
    basis <- kernel_basis <- c(kernel_covs, "-1")
    X <- model.matrix(reformulate(kernel_basis), data)
    X <- X / sqrt(1/varsweep * sum(apply(X, 2, var)))
  } else if (is.null(kernel_covs)) {
    # raw covs only
    covs <- raw_covs
    basis <- raw_basis <- c(raw_covs, "-1")
    X <- scale(model.matrix(reformulate(raw_basis), data))
    # X <- X / sqrt(1/varsweep * sum(apply(X, 2, var)))
  } else {
    # Raw + Kernel
    covs <- c(raw_covs, kernel_covs)
    raw_basis <- c(raw_covs, "-1")
    kernel_basis <- c(kernel_covs, "-1")
    basis <- c(raw_covs, kernel_covs, "-1")
    X_raw <- scale(model.matrix(reformulate(raw_basis), data))
    X_raw <- X_raw / sqrt(sum(apply(X_raw, 2, var)))
    X_kernel <- model.matrix(reformulate(kernel_basis), data)
    # varsweep is the intended target variance
    X_kernel <- X_kernel / sqrt(1/varsweep * sum(apply(X_kernel, 2, var)))
    X <- cbind(X_raw, X_kernel)
  }
  
  # balancing weights
  l <- estimate_regularization(data, covs = covs)
  
  bal_weights <- multilevel_qp(X, data$treat, Z = rep(1, nrow(data)), lambda = l, verbose= FALSE,
                               exact_global = FALSE, scale_sample_size = FALSE)
  data$wts <- pmax(bal_weights$weights, 0)
  data$wts[data$treat == 1] <- 1
  
  if (verbose) {
    print("balancing weights estimated")
  }
  
  bal.data <- as.data.frame(X)
  bal.names <- names(bal.data)
  n_covs <- length(bal.names)
  bal.data$treat <- data$treat
  bal.data$wts <- data$wts
  
  # ESS
  data0 <- data %>% dplyr::filter(treat==0)
  n.0 <- (sum(data0$wts)^2) / (sum(data0$wts^2))
  data1 <- data %>% dplyr::filter(treat==1)
  n.1 <- (sum(data1$wts)^2) / (sum(data1$wts^2))
  
  if (verbose) {
    print("ESS calculated")
  }
  
  ## Balance diagnostics
  
  data.var <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~var(.x))) %>% as.data.frame()
  
  c.var <- as.numeric(data.var[1,])
  t.var <- as.numeric(data.var[2,])
  c.var <- c.var[-1]
  t.var <- t.var[-1]
  pooled.var <- sqrt((t.var + c.var)/2)
  
  ## Balance
  um.wt <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~mean(.x))) %>% as.data.frame()
  
  bal.st <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~ weighted.mean(.x, wts))) %>% as.data.frame()
  
  um.wt.tab <- matrix(NA, length(bal.names), 3)
  um.wt.tab[,1] <- unlist(um.wt[1,-1]) 
  um.wt.tab[,2] <- unlist(um.wt[2,-1])                        
  um.wt.tab[,3] <- (unlist(um.wt[2,-1]) - unlist(um.wt[1,-1]))/pooled.var
  
  bal.st.tab <- matrix(NA, length(bal.names), 3)
  bal.st.tab[,1] <- unlist(bal.st[1,-1]) 
  bal.st.tab[,2] <- unlist(bal.st[2,-1])                        
  bal.st.tab[,3] <- (unlist(bal.st[2,-1]) - unlist(bal.st[1,-1]))/pooled.var # SMD
  
  um.wt.bias <- um.wt.tab[,3]
  bal.bias <- bal.st.tab[,3] 
  
  ## Bias Reduction
  pbr.bal.wt <- (1 - (mean(abs(bal.bias))/mean(abs(um.wt.bias))))*100
  if (verbose) {
    print("PBR calculated")
  }
  
  # Outcomes for weighting
  bal.wt <- geeglm(Y ~ treat, data = data, std.err = 'san.se', 
                   weights = wts, id=1:nrow(data),
                   corstr="independence")
  att.bw <- msm.out(bal.wt)
  att.bw
  sr.att.bw <- sum((data$treat -(1 - data$treat)*data$wts)*data$Y)/sum(data$treat)
  sr.att.bw # singly robust
  
  with(data, {
    mean(Y[treat==1]) - sum(wts[treat==0] * Y[treat==0]) / sum(wts[treat==0])
  })
  
  sum(data$wts[data$treat == 0])
  sum(data$treat)
  
  bias.bw <- sr.att.bw - true_att
  bias.bw
  # ATT of ACIC-17 around 0.118
  rel.bias <- compute_relative_bias(sr.att.bw, true_att)
  
  coverage.bw <- (true_att >= att.bw["lcl.treat"]) & (true_att <= att.bw["ucl.treat"])
  names(coverage.bw) <- NULL
  
  if (verbose) {
    print("simulation finished")
  }
  
  if (bias.bw > 3) {
    warning("Bias is large for balancing weights. Caution advised.")
  }
  
  data.frame(est = "bal.wgt", "feat_rep" = feat_rep, 
             "bias" = bias.bw, rel.bias=rel.bias, "cvg" = coverage.bw, 
             "pbr" = pbr.bal.wt, "ess" = n.0, "est.att" = sr.att.bw, se = att.bw["SE"],
             varsweep = varsweep)
}


eval_data_varsweep <- function(dat, pilot.dat, treat.true = 5, verbose = FALSE, simulation = TRUE) {
  
  data.c <- dat %>% filter(Z==0)
  data.t <- dat %>% filter(Z==1)
  
  # set.seed(1031)
  # train <- sample(nrow(data.c), round(.1*nrow(data.c)))
  # length(train)
  # data.c.train <- data.c[train, ]
  # data.c.test <- data.c[-train, ]
  data.c.train <- pilot.dat
  data <- rbind(data.t, data.c)
  # gets rid of the training sample used to select interactions & nonlinear covariates
  
  if (simulation) {
    covs <- names(dat)[!names(dat) %in% c("Z", "Y")]
  } else {
    covs <- c(
      "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
      "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
      "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
      "mathPre"
    )
  }
  
  #### Generate Random Forest features #####
  nc_rf <- nc_bart <- 10
  
  
  
  if (verbose) print("Starting RF Kernel")
  ## First, run random forest
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  rf_obj <- rf_kernel_matrix(model = rfmod, data = data, X = as.matrix(data %>% dplyr::select(-Y)),
                             n_components = max(nc_rf), verbose = TRUE)
  elbo_rf <- rf_obj$elbow
  
  data_rf_only_whole <- rf_obj$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_plus_whole <- rf_obj$data_rf
  # data_rf_mixed_whole <- rf_obj$features_mixed %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_K <- rf_obj$K %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  orig_varsweep <- sum(apply(rf_obj$features, 2, var))
  ##############################################################################
  
  ############################# Extract BART Kernel Features #################################
  if (verbose) print("Starting BART Kernel")
  
  # bart_obj <- bart_kernel_matrix(train = data.c.train, test = data, seed = 1022, verbose = TRUE, simulation = simulation)
  # K_generic <- bart_obj$kernel
  # bart.pca <- pca_bart(kernel = K_generic, data = data, X = as.matrix(data %>% dplyr::select(-Y)), 
  #                      n_components = max(nc_bart))
  # elbo_bart <- bart.pca$elbow
  # data_bart_only_whole <- bart.pca$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  # data_bart_plus_whole <- bart.pca$data_bart
  # data_bart_mixed_whole <- bart.pca$features_mixed %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  
  
  varsweeps <- c(10, 50, 100, 1000, 2500, 5000, 7500, 10000, 100000)#, orig_varsweep)
  varsweep_vec <- c(1 / varsweeps, 1, varsweeps)
  
  raw.scenarios <- expand.grid(fr = "raw", varsweep = NA)
  rf.scenarios <- rbind(expand.grid(fr = c("rf_plus"), varsweep = varsweep_vec), 
                        expand.grid(fr = c("rf_only"), varsweep = varsweep_vec))
  bart.scenarios <- rbind(expand.grid(fr = c("bart_plus"), varsweep = varsweep_vec), 
                          expand.grid(fr = c("bart_only"), varsweep = varsweep_vec))
  
  ### Run the actual simulation ###
  scenarios <- rbind(raw.scenarios, rf.scenarios)#, bart.scenarios)#, kbal.scenarios)
  if (verbose) print("Starting ATT estimation")
  out <- lapply(1:nrow(scenarios), function(i){
    fr <- scenarios[i,1]
    varsweep <- scenarios[i,2]
    nc <- nc_rf
    feat_rep <- paste(fr, varsweep, sep = "_")
    if (fr == "raw") {
      dataset <- data
      raw_covs <- covs
      kernel_covs <- NULL
      input_covs <- raw_covs
    } else if (fr == "rf_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_rf_plus_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "rf_only") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- kernel_covs
      dataset <- data_rf_only_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "rf_mixed") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- kernel_covs
      dataset <- data_rf_mixed_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "bart_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_bart_plus_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "bart_only") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- kernel_covs
      dataset <- data_bart_only_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "bart_mixed") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- kernel_covs
      dataset <- data_bart_mixed_whole %>% dplyr::select(Z, Y, all_of(input_covs))
    } 
    
    # compute ATT, balance metrics, error metrics, and so forth for each estimator 
    
    bw <- bal_weights_varsweep(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                           raw_covs = raw_covs, kernel_covs = kernel_covs, varsweep = varsweep,
                           verbose = FALSE) # simplex
    
    list(bw = bw, elbo_rf = elbo_rf, elbo_bart = NA)
  })
  out
}