library(foreign)
library(balancer)
library(dplyr)
library(ggplot2)
library(sandwich)
library(WeightIt)
library(splines)
library(Matrix)
library(glmnet)
library(stochtree)

rm(list=ls())

load("data/wsc-obs-study-data.RData")

source("R/estimationFunctions.R")
source("R/randomForestFeatures.R")
source("R/fitSplines.R")
source("R/BARTFeatures.R")
source("R/crossfit-funcs.R")

data.pre <- data.obs
data.pre$std_math <- scale(data.pre$mathPost)
load("data/covs.Rdata")

data.pre <- data.pre %>% dplyr::rename(treat = m.treat)

table(data.pre$treat)


data.c <- data.pre %>% filter(treat==0)
data.t <- data.pre %>% filter(treat==1)


n_repeat <- 6

set.seed(12)

num_comps <- c(5, 10, 25, 50, 100, 0)
num_comps_bart <- num_comps_rf <- num_comps[num_comps != 0]

single.fit <- function(pilot.dat, est.dat) {
  data <- rbind(data.t, est.dat)
  rf_fit_covs <- covs[covs != "-1"]
  rfmod <- randomForest::randomForest(mathPost ~ ., 
                                      data = pilot.dat %>% dplyr::select(all_of(rf_fit_covs), mathPost), 
                                      ntree = 500)
  
  ####### Extract the random forest features #######
  all_rf <- lapply(num_comps_rf, function(nc) {
    rf_obj <- extract_rf_features(rf_model = rfmod, data = data, n_components = nc, verbose = TRUE)
    data_rf_only <- rf_obj$features %>% 
      dplyr::mutate(mathPost = data$mathPost, treat = data$treat)
    var_exp <- rf_obj$explained_variance
    print(paste("Explained variance for first 5 components:", paste0(round(var_exp[1:5], 3), collapse = ", ")))
    data_rf_plus <- rf_obj$data_rf
    ol <- list(
      data_rf_only = data_rf_only,
      data_rf_plus = data_rf_plus
    )
    ol
  })
  names(all_rf) <- paste0("rf_", num_comps_rf)
  ##############################################################################
  
  ############################## Extract Gaussian Kernel Features #################################
  X_unscaled <- model.matrix(reformulate(c(rf_fit_covs, "-1")), data)
  # scale to have variance 1 but not mean 0
  X <- scale(X_unscaled)
  
  kbal_objs <- lapply(1:length(num_comps_rf), function(i) {
    print(i)
    nc <- num_comps[i]
    kbal_obj <- kbal::kbal(X, treatment = data$treat, numdims = nc)
    data_kbal_only <-  data.frame(
      kbal_obj$svdK$u[, 1:nc] %*% diag(sqrt(kbal_obj$svdK$d[1:nc]))
    )
    names(data_kbal_only) <- paste0("PC", 1:nc)
    var_exp <- kbal_obj$explained_variance
    data_kbal_plus <- cbind(data, data_kbal_only)
    ol <- list(
      data_kbal_only = data_kbal_only %>% dplyr::mutate(treat = data$treat, mathPost = data$mathPost),
      data_kbal_plus = data_kbal_plus
    )
    ol
  }
  )
  names(kbal_objs) <- paste0("kbal_", num_comps_rf)
  
  data_raw <- data #data.pre
  
  tree.fr <- c("bart_only", "bart_plus", "rf_only", "rf_plus", "kbal_only", "kbal_plus")
  tree.fr <- c("rf_only", "rf_plus", "kbal_only", "kbal_plus")
  
  scen.init <- expand.grid(fr = c(tree.fr, "raw", "regression"), 
                           nc = num_comps)
  
  # delete the rows where the first column is "raw" or "regression" and the second column is NOT 0
  scenarios <- scen.init %>% 
    dplyr::filter(!(fr %in% c("raw", "regression") & nc != 0)) %>% 
    dplyr::filter(!(fr %in% tree.fr & nc == 0)) 
  scenarios
  
  out <- map(
    1:nrow(scenarios),
    function(i) {
      fr <- scenarios[i, 1] 
      nc <- scenarios[i, 2]
      outcome <- "mathPost"
      print(paste("Running scenario", i, "with ncomp", nc, "and fr", fr))
      
      if (fr != "regression") {
        if (fr == "raw") {
          input_data <- data_raw
          input_covs <- covs
        } else if (fr == "rf_plus") {
          input_data <- all_rf[[paste0("rf_", nc)]]$data_rf_plus
          input_covs <- c(covs, paste0("PC", 1:nc))
        } else if (fr == "rf_only") {
          input_data <- all_rf[[paste0("rf_", nc)]]$data_rf_only
          input_covs <- c(paste0("PC", 1:nc), "-1")
        } else if(fr == "bart_plus") {
          input_data <- all_bart[[paste0("bart_", nc)]]$data_bart_plus
          input_covs <- c(covs, paste0("BART_", 1:nc))
        } else if (fr == "bart_only") {
          input_data <- all_bart[[paste0("bart_", nc)]]$data_bart_only
          input_covs <- c(paste0("BART_", 1:nc), "-1")
        } else if (fr == "kbal_only") {
          input_data <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_only
          input_covs <- paste0("PC", 1:nc, "-1")
        } else if (fr == "kbal_plus") {
          input_data <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_plus
          input_covs <- c(covs, paste0("PC", 1:nc))
        }
        data.c <- input_data %>% filter(treat==0) %>% dplyr::mutate(std_math = scale(mathPost))
        lambda.reg <- lm(reformulate(input_covs, response="std_math"), data=data.c)
        l <- var(lambda.reg$resid)
        print(paste("Regularization parameter =", round(l, 2), "for", fr, "with ncomp", nc))
        bw.l2 <- balancingWeights(data = input_data, outcome_name = outcome, 
                                  covs = input_covs, feat_rep = fr, l = l, type = "l2", 
                                  verbose = T)
        
        # bw.inf <- balancingWeights(data = input_data, outcome_name = outcome, 
        #                            covs = input_covs, feat_rep = fr, l = l, type = "inf", 
        #                            verbose = F)
        
        ipw <- logisticIPW(data = input_data, outcome_name = outcome, covs = input_covs,
                           feat_rep = fr, verbose = F) 
        
        ols <- outcomeRegression(data = input_data, outcome_name = outcome, 
                                 covs = input_covs, feat_rep = fr, 
                                 type = "ols", verbose = F) %>% unlist()
        
        ridge <- outcomeRegression(data = input_data, outcome_name = outcome, 
                                   covs = input_covs, feat_rep = fr, 
                                   type = "ridge", verbose = F) %>% unlist()
        
        aug.l2 <- augmentedBalWeights(data = input_data, outcome_name = outcome, covs = input_covs, 
                                      feat_rep = fr, out.mod = "ols", ps.mod = "l2", verbose = F) %>% unlist()
        
        aug.vanilla <- augmentedBalWeights(data = input_data, outcome_name = outcome, covs = input_covs, 
                                           feat_rep = fr, out.mod = "ols", ps.mod = "ipw", verbose = FALSE) %>% unlist()
        
        # combine the results into a single dataframe
        results <- rbind(bw.l2$res %>% unlist(), 
                         #bw.inf$res %>% unlist(), 
                         ipw$res %>% unlist(), 
                         ols, ridge, aug.l2, aug.vanilla) %>% 
          data.frame() %>% 
          tibble() %>% 
          mutate(ncomp = nc, outcome = outcome, feat_rep = fr)
        names(results)[1:4] <- c("estimate", "lcl", "hcl", "se")
        ll <- list(res = results, p_bw = bw.l2$p, p_ipw = ipw$p)#, p_inf = bw.inf$p)
        print("DONE")
        return(ll)
      } else {
        input_covs <- rf_fit_covs
        input_data <- data_raw %>% dplyr::select(all_of(input_covs), treat, all_of(outcome))
        # run vanilla random forest to estimate the ATT
        data.0 <- input_data %>% filter(treat == 0)
        data.1 <- input_data %>% filter(treat == 1)
        rfmod <- randomForest::randomForest(as.formula(paste(outcome, "~ .")), 
                                            data = data.0, ntree = 100)
        rf_pred <- predict(rfmod, newdata = data.1)
        est.att <- mean(data.1[[outcome]]) - mean(rf_pred)
        lcl <- hcl <- se <- pbr <- ess <- NA
        out.df <- data.frame(
          estimate = as.character(est.att),
          lcl = lcl,
          hcl = hcl,
          se = se,
          pbr = pbr,
          ess = ess,
          est = "rf.reg",
          ncomp = nc,
          outcome = outcome,
          feat_rep = "rf.reg"
        ) %>% tibble()
        return(list(res = out.df))
      }
    }
  )
  res_df <- map(out, function(x) x$res) %>% bind_rows()
  res_df
  
  res_df <- res_df %>% 
    mutate_at(vars(estimate, lcl, hcl, se, pbr, ess), as.numeric) %>% 
    mutate(method = recode(est, 
                           bw.l2 = "Balancing Weights - L2", 
                           bw.inf = "Balancing Weights - L-Inf",
                           ipw = "IP Weights",
                           ols = "OLS Outcome Regression",
                           ols_ipw = "Augmented IP Weights with OLS",
                           ols_l2 = "Augmented L2 BW with Ridge",
                           rf.reg = "Random Forest Regression"))
  res_df
}

out <- lapply(1:n_repeat, function(i) {
  sample_split <- sample(nrow(data.c), round(.1*nrow(data.c)))
  length(sample_split)
  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]
  
  fit.2.1 <- single.fit(pilot.dat = data.c2, est.dat = data.c1)
  fit.1.2 <- single.fit(pilot.dat = data.c1, est.dat = data.c2)
  
  # average each entry across the two fits. take entry-wise everage between fit.1.2[i,j] and fit.2.1[i,j]. 
  # in other words, give me (fit.1.2[i,j] + fit.2.1[i,j]) / 2
  # Assume fit.1.2 and fit.2.1 are already loaded and have the same structure
  cfit.df <- fit.1.2  # initialize with one of the data frames
  
  # Identify numeric columns
  num_cols <- sapply(fit.1.2, is.numeric)
  
  # Average numeric columns
  cfit.df[, num_cols] <- (fit.1.2[, num_cols] + fit.2.1[, num_cols]) / 2
  
  cfit.df$se <- sqrt(fit.1.2$se^2 / 2 + fit.2.1$se^2 / 2)
  
  cfit.df$lcl <- cfit.df$estimate - 1.96 * cfit.df$se
  cfit.df$hcl <- cfit.df$estimate + 1.96 * cfit.df$se
  cfit.df
})


cfit.df <- average_cross_fit(out)

save(cfit.df, file = "results/wsc-math-cfit.RData")
