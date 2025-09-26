library(foreign)
library(balancer)
library(dplyr)
library(ggplot2)
library(sandwich)
library(WeightIt)
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

setwd("~/Desktop/BalWeights/forest-kbal/wsc")

load("data/wsc-obs-study-data.RData")

setwd("~/Desktop/BalWeights/forest-kbal/wsc")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")


outcome <- "math"

if (outcome == "math") {
  data.obs$Y <- data.obs$mathPost
  data.pre <- data.obs
  data.pre$std_math <- scale(data.pre$mathPost)
  data.pre <- data.pre %>% dplyr::rename(Z = m.treat)
} else if (outcome == "vocab") {
  data.obs$Y <- data.obs$vocabPost
  data.pre <- data.obs
  data.pre$std_math <- scale(data.pre$vocabPost)
  data.pre <- data.pre %>% dplyr::rename(Z = v.treat)
}



covs <- c(
  "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
  "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
  "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
  "mathPre"
)

data.pre.log <- data.pre

naive.dim <- mean(data.pre$Y[data.pre$Z == 1]) - mean(data.pre$Y[data.pre$Z == 0])
naive.dim

log_trans <- TRUE

if (log_trans == TRUE) {
  # log transform
  data.pre.log <- data.pre.log %>%
    mutate(across(
      .cols = where(is.numeric) & !any_of(c("Z", "Y")),
      .fns = ~ if (n_distinct(.x) >= 3 && min(.x, na.rm = TRUE) >= 0) {
        log(.x + 1)
      } else {
        .x
      }
    ))
} else {
  # exp transform
  data.pre.log <- data.pre.log %>%
    mutate(across(
      .cols = where(is.numeric) & !any_of(c("Z", "Y")),
      .fns = ~ if (n_distinct(.x) >= 2 && min(.x, na.rm = TRUE) >= 0) {
        exp(.x)
      } else {
        .x
      }
    ))
}

data.c <- data.pre %>% filter(Z == 0)
data.t <- data.pre %>% filter(Z == 1)

data.c.log <- data.pre.log %>% filter(Z == 0)
data.t.log <- data.pre.log %>% filter(Z == 1)

n_repeat <- 5

set.seed(12)

full.dat <- bind_rows(data.c.log, data.t.log)

raw.all <- balancingWeights(data = full.dat, true_att = 0, feat_rep = "raw", 
                            raw_covs = covs, kernel_covs = NULL, verbose = TRUE)

raw.all

out_filename <- paste0("logs/wsc-xfit-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

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

log_message <- paste("Starting WSC XFIT with log_trans", ifelse(log_trans == TRUE, "LOG", "EXP"), "\n")
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
}, mc.set.seed = TRUE, mc.cores = numCores - 1)


if (outcome == "math") {
  if (log_trans == TRUE) {
    save(out, file = "results/wsc-math-xfit.RData")
  } else {
    save(out, file = "results/wsc-math-xfit-exp.RData")
  }
} else if (outcome == "vocab") {
  if (log_trans == TRUE) {
    save(out, file = "results/wsc-vocab-xfit.RData")
  } else {
    save(out, file = "results/wsc-vocab-xfit-exp.RData")
  }
}
