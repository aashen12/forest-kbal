eval_data <- function(dat, pilot.dat, treat.true = 5, verbose = FALSE, simulation = TRUE) {
  
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
  if (simulation) {
    nc_rf <- c(5, 10, 25, 50, 100)
  } else {
    nc_rf <- 2:50
  }
  rf.scenarios.partial <- expand.grid(fr = c("rf_only", "rf_plus"), ncomp = nc_rf)
  rf.scenarios.ker <- expand.grid(fr = c("rf_K"), ncomp = nrow(data))
  rf.scenarios <- rbind(rf.scenarios.partial, rf.scenarios.ker)
  ## First, run random forest
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  rf_obj <- rf_kernel_matrix(model = rfmod, data = data, X = as.matrix(data %>% dplyr::select(-Y)),
                             n_components = max(nc_rf), verbose = TRUE)
  
  data_rf_only_whole <- rf_obj$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_plus_whole <- rf_obj$data_rf
  data_rf_K <- rf_obj$K %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  
  
  ############################# Extract Gaussian Kernel Features #################################
  kbal.scenarios <- expand.grid(fr = c("kbal_only", "kbal_plus"), ncomp = nc_rf)
  X_unscaled <- model.matrix(reformulate(c(covs, "-1")), data)
  # scale to have variance 1 but not mean 0
  X <- scale(X_unscaled)
  kbal_objs <- lapply(1:length(nc_rf), function(i) {
    nc <- nc_rf[i]
    kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE)
    data_kbal_only <-  data.frame(
      kbal_obj$svdK$u[, 1:nc] %*% diag(sqrt(kbal_obj$svdK$d[1:nc]))
    )
    names(data_kbal_only) <- paste0("PC", 1:nc)
    var_exp <- kbal_obj$explained_variance
    data_kbal_plus <- cbind(data, data_kbal_only)
    ol <- list(
      data_kbal_only = data_kbal_only %>% dplyr::mutate(treat = data$Z, Y = data$Y),
      data_kbal_plus = data_kbal_plus
    )
    ol
  }
  )
  names(kbal_objs) <- paste0("kbal_", nc_rf)
  if (verbose) print("KBal Finished!")
  ###########################################################################
  
  raw.scenarios <- expand.grid(fr = "raw", ncomp = 0)
  ### Run the actual simulation ###
  scenarios <- rbind(raw.scenarios, rf.scenarios, kbal.scenarios)
  
  out <- lapply(1:nrow(scenarios), function(i){
    fr <- scenarios[i,1]
    nc <- scenarios[i,2]
    feat_rep <- paste(fr, nc, sep = "_")
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
    } else if (fr == "rf_K") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_rf_K %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "kbal_only") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_only
    } else if (fr == "kbal_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_plus
    }
    
    # compute ATT, balance metrics, error metrics, and so forth for each estimator 
    
    bw <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                           raw_covs = raw_covs, kernel_covs = kernel_covs,
                           verbose = FALSE) # simplex
    
    ipw <- logisticIPW(data = dataset, true_att = treat.true, 
                       feat_rep = feat_rep, verbose = FALSE)
    
    rf <- outcomeRegression(data = dataset, true_att = treat.true, 
                            feat_rep = feat_rep, type = "rf", verbose = FALSE)
    
    aug.l2 <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                                  raw_covs = raw_covs, kernel_covs = kernel_covs,
                                  out.mod = "ols", ps.mod = "l2", verbose = FALSE)
    aug.vanilla <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                                       raw_covs = raw_covs, kernel_covs = kernel_covs,
                                       out.mod = "ols", ps.mod = "ipw", verbose = FALSE)
    
    
    list(bw = bw, ipw = ipw, 
         rf = rf,
         aug.l2 = aug.l2, aug.vanilla = aug.vanilla)
  })
  out
}