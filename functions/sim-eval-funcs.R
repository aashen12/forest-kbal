eval_data <- function(dat, pilot.dat, treat.true = 5, verbose = FALSE, covs = NULL, dataset = "soldiering") {
  
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
  
  if (dataset == "simulation") {
    covs <- names(dat)[!names(dat) %in% c("Z", "Y")]
  } else if (dataset == "soldiering") {
    # FILL IN FROM CONSEQUENCES PAPER
    AC_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$",
                    names(raw.data), value = TRUE)
    ACG_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$|^G[1-8]_",
                     names(raw.data), value = TRUE)

    ctrls <- c("fthr_ed0", "fthr_ed4", "fthr_ed7", "mthr_ed0", "mthr_ed4", "mthr_ed7",
               "no_fthr96", "no_mthr96", "orphan96", "hh_fthr_frm", "hh_size96",
               "hh_size96_2", "hh_size96_3", "hh_size96_4", "hh_land", "hh_land_2",
               "hh_land_3", "hh_land_4", "hh_cattle", "hh_cattle_2", "hh_cattle_3",
               "hh_cattle_4", "hh_stock", "hh_stock_2", "hh_stock_3", "hh_stock_4", "hh_plow")

    ctrl_hh <- c("age", "hh_fthr_frm", "hh_size96", "hh_land", "landrich",
                 "hh_cattle", "hh_stock", "hh_plow")

    ctrl_i <- c("fthr_ed0", "fthr_ed", "mthr_ed0", "mthr_ed", "no_fthr96",
                "no_mthr96", "orphan96")
    # 
    # ## Covariates include geographical region, age in 1996, 
    # ### father’s education, mother’s education, whether the parents had died during 
    # ### or before 1996, whether the father is a farmer, and household size in 1996.
    # 
    # ## AC: age, geog region
    # ## ctrl_i: parental education, orphan, no parent, father farmer
    # ## orphan96: orphan in 1996
    # ## hh_fthr_frm: father farmer
    # ## hh_size96: household size in 1996
    # 
    # covs <- c(AC_vars, "fthr_ed", "mthr_ed", "orphan96", "hh_fthr_frm", "hh_size96")
  } else if (dataset == "wsc") {
    covs <- c(
      "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
      "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
      "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
      "mathPre"
    )
  } else if (dataset == "lalonde") {
    covs <- names(dat)[!names(dat) %in% c("Z", "Y")]
  } else {
    print("Dataset not recognized, defaulting to simulation covariates")
    covs <- names(dat)[!names(dat) %in% c("Z", "Y")]
  }
  
  if (dataset == "simulation") {
    nc_rf <- nc_bart <- c(2, 5, 10, 15, 25, 50, 100)
  } else {
    nc_rf <- nc_bart <- 2:15
  }
  
  ############################# Extract BART Kernel Features #################################
  print("Starting BART Kernel")
  # bart.scenarios.partial <- expand.grid(fr = c("bart_only", "bart_plus", "bart_mixed"), ncomp = nc_bart)
  bart.scenarios.partial <- expand.grid(fr = c("bart_only", "bart_plus"), ncomp = nc_bart)
  bart.scenarios.ker <- expand.grid(fr = c("bart_K", "bart_K_plus"), ncomp = nrow(data))
  bart.scenarios <- rbind(bart.scenarios.partial, bart.scenarios.ker)
  bart_obj <- bart_kernel_matrix(train = data.c.train, test = data, seed = 1022, verbose = TRUE, dataset = "soldiering", covs = covs)
  K_generic <- bart_obj$kernel
  bart.pca <- pca_bart(kernel = K_generic, data = data, X = as.matrix(data %>% dplyr::select(-Y)), 
                       n_components = max(nc_bart))
  bart_expl_var <- bart.pca$explained_variance
  elbo_bart <- bart.pca$elbow
  data_bart_only_whole <- bart.pca$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_bart_plus_whole <- bart.pca$data_bart
  # data_bart_mixed_whole <- bart.pca$features_mixed %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_bart_K <- data.frame(as.matrix(K_generic))
  names(data_bart_K) <- paste0("PC", 1:ncol(data_bart_K))
  data_bart_K <- data_bart_K %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  
  data_bart_K_plus <- cbind(data_bart_K, data %>% dplyr::select(all_of(covs)))
  ##############################################################################
  
  ############################# Extract Gaussian Kernel Features #################################
  if (verbose) print("Starting Kbal Kernel")
  kbal.scenarios.partial <- expand.grid(fr = c("kbal_only", "kbal_plus"), ncomp = nc_rf)
  kbal.scenarios.ker <- expand.grid(fr = c("kbal_K", "kbal_K_plus"), ncomp = nrow(data))
  kbal.scenarios <- rbind(kbal.scenarios.partial, kbal.scenarios.ker)
  X_unscaled <- model.matrix(reformulate(c(covs, "-1")), data)
  # scale to have variance 1 but not mean 0
  X <- scale(X_unscaled)
  
  
  
  kbal_objs <- lapply(1:length(nc_rf), function(i) {
    nc <- nc_rf[i]
    if (dataset == "simulation") {
      sim2 <- length(unique(X[,6])) == 2 # X6 in sim2 is bernoulli, we just need a quick check for Sim1 or Sim2
      print(paste("Sim2:", sim2))
      kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,
                             mixed_data = sim2, cat_columns = ifelse(sim2, "X6", NULL))
    } else if (dataset == "lalonde") {
      kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,
                             mixed_data = TRUE,
                             cat_columns = c("black", "hisp", "married", "u74", "u75", "nodegr", "educ")
                             )
    } else if (dataset == "wsc") {
      kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,
                             mixed_data = TRUE,
                             cat_columns = c("female", "white", "black", "asian", "hisp", "married",
                                             "collegeS", "collegeM", "collegeD", "calc",
                                             "mathLike", "income"))
    } else if (dataset == "soldiering") {
      num_vals <- apply(data %>% dplyr::select(all_of(covs)), 2, function(x) length(unique(x)))
      discrete_covs <- names(num_vals[num_vals <= 10])
      kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,
                             mixed_data = TRUE,
                             cat_columns = discrete_covs)
    } else {
      kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE)
    }
    kbal_only_weights <- kbal_obj$w
    data_kbal_only <-  data.frame(
      kbal_obj$svdK$u[, 1:nc] %*% diag(sqrt(kbal_obj$svdK$d[1:nc]))
    )
    names(data_kbal_only) <- paste0("PC", 1:nc)
    data_kbal_only$kbal_weights <- kbal_only_weights
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
  
  if (verbose) print("KBal PC Finished. Moving onto KBal Full K")
  
  ### EXTRACT kbal_K ###
  if (dataset == "simulation") {
    sim2 <- length(unique(X[,6])) == 2 # X6 in sim2 is bernoulli, we just need a quick check for Sim1 or Sim2
    print(paste("Sim2:", sim2))
    kbal_K_obj <- kbal::kbal(X, treatment = data$Z, printprogress = FALSE,
                           mixed_data = sim2, cat_columns = ifelse(sim2, "X6", NULL))
  } else if (dataset == "lalonde") {
    kbal_K_obj <- kbal::kbal(X, treatment = data$Z, printprogress = FALSE,
                           mixed_data = TRUE,
                           cat_columns = c("black", "hisp", "married", "u74", "u75", "nodegr", "educ")
    )
  } else if (dataset == "wsc") {
    kbal_K_obj <- kbal::kbal(X, treatment = data$Z, printprogress = FALSE,
                           mixed_data = TRUE,
                           cat_columns = c("female", "white", "black", "asian", "hisp", "married",
                                           "collegeS", "collegeM", "collegeD", "calc",
                                           "mathLike", "income"))
  } else if (dataset == "soldiering") {
    num_vals <- apply(data %>% dplyr::select(all_of(covs)), 2, function(x) length(unique(x)))
    discrete_covs <- names(num_vals[num_vals <= 10])
    kbal_K_obj <- kbal::kbal(X, treatment = data$Z, printprogress = FALSE,
                           mixed_data = TRUE,
                           cat_columns = discrete_covs)
  } else {
    kbal_K_obj <- kbal::kbal(X, treatment = data$Z, printprogress = FALSE)
  }
  kbal_expl_var <- NULL
  data_kbal_K <- data.frame(kbal_K_obj$K) 
  names(data_kbal_K) <- paste0("PC", 1:ncol(data_kbal_K))
  data_kbal_K <- data_kbal_K %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_kbal_K_plus <- cbind(data_kbal_K, data %>% dplyr::select(all_of(covs)))
  
  if (verbose) print("KBal Finished!")
  ##############################################################################
  
  ##################### Generate Random Forest features ############################
  
  if (verbose) print("Starting RF Kernel")
  # rf.scenarios.partial <- expand.grid(fr = c("rf_only", "rf_plus", "rf_mixed"), ncomp = nc_rf)
  rf.scenarios.partial <- expand.grid(fr = c("rf_only", "rf_plus"), ncomp = nc_rf)
  rf.scenarios.ker <- expand.grid(fr = c("rf_K", "rf_K_plus"), ncomp = nrow(data))
  rf.scenarios <- rbind(rf.scenarios.partial, rf.scenarios.ker)
  #rf.scenarios <- rf.scenarios.partial
  ## First, run random forest
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  rf_obj <- rf_kernel_matrix(model = rfmod, data = data, X = as.matrix(data %>% dplyr::select(-Y)),
                             n_components = max(nc_rf), verbose = TRUE)
  elbo_rf <- rf_obj$elbow
  rf_expl_var <- rf_obj$explained_variance
  data_rf_only_whole <- rf_obj$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_plus_whole <- rf_obj$data_rf
  # data_rf_mixed_whole <- rf_obj$features_mixed %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_K <- rf_obj$K %>% dplyr::mutate(Y = data$Y, Z = data$Z)
  # append raw covs with data_rf_K
  data_rf_k_plus <- cbind(data_rf_K, data %>% dplyr::select(all_of(covs)))
  ##############################################################################
  
  raw.scenarios <- expand.grid(fr = "raw", ncomp = 0)
  ### Run the actual simulation ###
  scenarios <- rbind(raw.scenarios, rf.scenarios, bart.scenarios, kbal.scenarios)
  if (verbose) print("Starting ATT estimation")
  out <- lapply(1:nrow(scenarios), function(i){
    fr <- scenarios[i,1]
    nc <- scenarios[i,2]
    feat_rep <- paste(fr, nc, sep = "_")
    print(feat_rep)
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
    } else if (fr == "rf_K") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_rf_K %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "bart_K") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_bart_K %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "kbal_K") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_kbal_K %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "rf_K_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_rf_k_plus %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "bart_K_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_bart_K_plus %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "kbal_K_plus") {
      raw_covs <- covs
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- data_kbal_K_plus %>% 
        dplyr::select(Z, Y, all_of(input_covs))
    } else if (fr == "kbal_only") {
      raw_covs <- NULL
      kernel_covs <- c(paste0("PC", 1:nc))
      input_covs <- c(raw_covs, kernel_covs)
      dataset <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_only
      att_kbal <- with(dataset, mean(Y[treat==1]) - weighted.mean(Y[treat==0], w = kbal_weights[treat==0]))
      bal.wt <- geeglm(Y ~ treat, data = dataset, std.err = 'san.se', 
                       weights = kbal_weights, id=1:nrow(dataset),
                       corstr="independence")
      att.bw <- msm.out(bal.wt)
      att.bw
      sr.att.bw <- sum((dataset$treat -(1 - dataset$treat)*dataset$kbal_weights)*dataset$Y)/sum(dataset$treat)
      sr.att.bw # singly robust
      bias.bw <- sr.att.bw - treat.true
      bias.bw
      # ATT of ACIC-17 around 0.118
      rel.bias <- compute_relative_bias(sr.att.bw, treat.true)
      
      coverage.bw <- (treat.true >= att.bw["lcl.treat"]) & (treat.true <= att.bw["ucl.treat"])
      names(coverage.bw) <- NULL
      
      kbal.bw.df <- data.frame(est = "kbal.bw", "feat_rep" = feat_rep, 
                 "bias" = bias.bw, rel.bias=rel.bias, "cvg" = coverage.bw, 
                 "pbr" = NA, "ess" = NA, "est.att" = sr.att.bw, se = att.bw["SE"])
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
    
    if(fr %in% c("kbal_only")) {
      kbal_bw <- kbal.bw.df
    } else {
      kbal_bw <- data.frame(est = "kbal.bw", "feat_rep" = feat_rep, 
                            "bias" = NA, rel.bias=NA, "cvg" = NA, 
                            "pbr" = NA, "ess" = NA, "est.att" = NA, se = NA)
    }
      
    
    ipw <- logisticIPW(data = dataset, true_att = treat.true, 
                       feat_rep = feat_rep, covs = input_covs, verbose = FALSE)
    
    # rf <- outcomeRegression(data = dataset, true_att = treat.true, 
    #                         feat_rep = feat_rep, type = "rf", covs = input_covs, verbose = FALSE)
    
    aug.l2 <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                                  raw_covs = raw_covs, kernel_covs = kernel_covs,
                                  out.mod = "ols", ps.mod = "l2", verbose = FALSE)
    aug.vanilla <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, 
                                       raw_covs = raw_covs, kernel_covs = kernel_covs,
                                       out.mod = "ols", ps.mod = "ipw", verbose = FALSE)
    
    
    list(bw = bw, ipw = ipw, 
         # rf = rf, 
         kbal_bw = kbal_bw,
         aug.l2 = aug.l2, aug.vanilla = aug.vanilla, 
         elbo_rf = elbo_rf, 
         elbo_bart = elbo_bart,
         rf_expl_var = rf_expl_var, bart_expl_var = bart_expl_var, kbal_expl_var = kbal_expl_var)
  })
  out
}