library(foreign)
library(balancer)
library(dplyr)
library(ggplot2)
library(sandwich)
library(splines)
library(GenericML)
library(Matrix)
library(glmnet)
library(stochtree)
library(foreach)
library(doParallel)
library(parallel)
library(future.apply)
library(randomForest)
library(kbal)
library(geepack)
library(irlba)

rm(list=ls())

setwd("~/Desktop/BalWeights/forest-kbal/lalonde")
setwd("~/Desktop/BalWeights/forest-kbal/lalonde")


source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")


data(lalonde, package = "kbal")
data.obs <- lalonde 
data.obs <- data.obs %>% dplyr::rename(Z = nsw, Y = re78) %>% dplyr::select(-race_ethnicity, -u78)

data.obs$Y <- log1p(data.obs$Y)  # log(1 + Y) to handle zeros

table(data.obs$Z)

head(data.obs)

covs <- names(data.obs)[!names(data.obs) %in% c("Z", "Y")]
length(covs)
covs

data.obs.log <- data.obs

naive.dim <- mean(data.obs$Y[data.obs$Z == 1]) - mean(data.obs$Y[data.obs$Z == 0])
naive.dim

log_trans <- TRUE

if (log_trans == TRUE) {
  # log transform
  data.obs.log$re74 <- log1p(data.obs.log$re74)
  data.obs.log$re75 <- log1p(data.obs.log$re75)
  data.obs.log$age <- log1p(data.obs.log$age)
  # data.obs.log$Y <- log1p(data.obs.log$Y)
} else {
  # exp transform
  data.obs.log <- data.obs.log %>%
    mutate(across(
      .cols = where(is.numeric) & !any_of(c("Z", "Y", "re74", "re75")),
      .fns = ~ if (n_distinct(.x) >= 2 && min(.x, na.rm = TRUE) >= 0) {
        exp(.x)
      } else {
        .x
      }
    ))
}

data.c <- data.obs %>% filter(Z == 0)
data.t <- data.obs %>% filter(Z == 1)

data.c.log <- data.obs.log %>% filter(Z == 0)
data.t.log <- data.obs.log %>% filter(Z == 1)

n_repeat <- 5

set.seed(12)

full.dat <- bind_rows(data.c, data.t)
full.dat.log <- bind_rows(data.c.log, data.t.log)

# check if any columns of full.dat.log contain Inf
any(sapply(full.dat.log, function(x) any(is.infinite(x))))

# which columns contain Inf?
which(sapply(full.dat.log, function(x) any(is.infinite(x))))



raw.all <- balancingWeights(data = full.dat, true_att = 0, feat_rep = "raw", 
                            raw_covs = covs, kernel_covs = NULL, verbose = TRUE)

raw.all

out_filename <- paste0("logs/lalonde-xfit-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

write("", out_filename, append = FALSE)   ##### ADDED (overwrite existing file)




single.fit <- function(pilot.dat, est.dat, treat.dat) {
  analysis.dat <- bind_rows(treat.dat, est.dat)
  sing.fit.out <- eval_data(dat = analysis.dat, 
                            pilot.dat = pilot.dat, 
                            treat.true = 0, verbose = TRUE, simulation = FALSE)
  sing.fit.out
}

process_finished_sim <- function(edat, id = 999) {
  out <- lapply(1:length(edat), function(i) {
    resi <- edat[[i]]
    elbo_rf <- resi$elbo_rf
    elbo_bart <- resi$elbo_bart
    resi_rest <- resi[!names(resi) %in% c("elbo_rf", "elbo_bart")]
    dplyr::bind_rows(resi_rest) %>% dplyr::mutate(elbo_rf = elbo_rf, elbo_bart = elbo_bart)
  })
  dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
}

run_cross_fit <- function(data.c1, data.c2, treat.dat, trans, id = 999) {
  fit.2.1.raw <- single.fit(pilot.dat = data.c2, est.dat = data.c1, treat.dat = treat.dat) 
  print("Finished first fit")
  fit.1.2.raw <- single.fit(pilot.dat = data.c1, est.dat = data.c2, treat.dat = treat.dat)
  print("Finished second fit") 

  
  # average each entry across the two fits. take entry-wise everage between fit.1.2[i,j] and fit.2.1[i,j]. 
  # in other words, give me (fit.1.2[i,j] + fit.2.1[i,j]) / 2
  # Assume fit.1.2 and fit.2.1 are already loaded and have the same structure
  
  fit.2.1 <- process_finished_sim(fit.2.1.raw, id = id)
  fit.1.2 <- process_finished_sim(fit.1.2.raw, id = id)
  
  
  cfit.df <- fit.1.2  # initialize with one of the data frames
  
  # Identify numeric columns
  num_cols <- sapply(fit.1.2, is.numeric)
  
  # Average numeric columns
  cfit.df[, num_cols] <- (fit.1.2[, num_cols] + fit.2.1[, num_cols]) / 2
  
  cfit.df$se <- sqrt(fit.1.2$se^2 / 2 + fit.2.1$se^2 / 2)
  
  cfit.df$lcl <- cfit.df$est.att - 1.96 * cfit.df$se
  cfit.df$hcl <- cfit.df$est.att + 1.96 * cfit.df$se
  
  cfit.df$trans <- trans
  
  return(cfit.df)
}


numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

log_message <- paste("Starting lalonde XFIT with log_trans", ifelse(log_trans == TRUE, "LOG", "EXP"), "\n")
cat(log_message, file = out_filename, append = TRUE)

out <- parallel::mclapply(1:n_repeat, function(i) {
  log_message <- paste("Starting repeat number", i, "at", Sys.time(), "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))
  
  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]
  
  data.c1.log <- data.c.log[sample_split, ]
  data.c2.log <- data.c.log[-sample_split, ]
  
  out.notrans <- run_cross_fit(data.c1 = data.c1, data.c2 = data.c2, treat.dat = data.t, trans = "none", id = i)
  out.log <- run_cross_fit(data.c1 = data.c1.log, data.c2 = data.c2.log, treat.dat = data.t.log, trans = "log", id = i)
  
  log_message <- paste("Finished repeat number", i, "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  out_df <- bind_rows(out.notrans, out.log)
  out_df
}, mc.set.seed = TRUE, mc.cores = numCores-1)

save(out, file = paste0("results/lalonde-xfit-kbal-logged.RData"))
