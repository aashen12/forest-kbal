library(foreign)
library(balancer)
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
library(tidyverse)

rm(list=ls())

setwd("~/Desktop/BalWeights/forest-kbal/soldiering")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

setwd("~/Desktop/BalWeights/forest-kbal/soldiering")

# raw.data <- read_csv("blattman.csv")

raw.data <- read_csv("blattman_claude.csv") %>% dplyr::rename(log.wage = lwage_mo)


## For the wage variable, 504 observed values

hist(raw.data$log.wage)
hist(exp(raw.data$log.wage) - 1)

mean(raw.data$log.wage[raw.data$abd == 1], na.rm = TRUE) - mean(raw.data$log.wage[raw.data$abd == 0], na.rm = TRUE)
mean(raw.data$educ[raw.data$abd == 1], na.rm = TRUE) - mean(raw.data$educ[raw.data$abd == 0], na.rm = TRUE)


df <- raw.data %>% 
  dplyr::mutate(wage_mo = exp(log.wage) - 1) %>%
  dplyr::rename(Z = abd, Y = distress)
table(df$Z)
dim(df)
# names(df)


# df$Y <- log1p(df$Y)


AC_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$", 
                names(raw.data), value = TRUE)
ACG_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$|^G[1-8]_", 
                 names(raw.data), value = TRUE)

length(ACG_vars) - length(AC_vars)

length(ACG_vars) - (2 + 15 + 31)

ctrls <- c("fthr_ed0", "fthr_ed4", "fthr_ed7", "mthr_ed0", "mthr_ed4", "mthr_ed7",
           "no_fthr96", "no_mthr96", "orphan96", "hh_fthr_frm", "hh_size96", 
           "hh_size96_2", "hh_size96_3", "hh_size96_4", "hh_land", "hh_land_2", 
           "hh_land_3", "hh_land_4", "hh_cattle", "hh_cattle_2", "hh_cattle_3", 
           "hh_cattle_4", "hh_stock", "hh_stock_2", "hh_stock_3", "hh_stock_4", "hh_plow")

ctrl_hh <- c("age", "hh_fthr_frm", "hh_size96", "hh_land", "landrich", 
             "hh_cattle", "hh_stock", "hh_plow")

ctrl_i <- c("fthr_ed0", "fthr_ed", "mthr_ed0", "mthr_ed", "no_fthr96", 
            "no_mthr96", "orphan96")

full_covs <- c(AC_vars, ACG_vars, ctrls, ctrl_hh, ctrl_i)

## Covariates include geographical region, age in 1996, 
### father’s education, mother’s education, whether the parents had died during 
### or before 1996, whether the father is a farmer, and household size in 1996.

## AC: age, geog region
## ctrl_i: parental education, orphan, no parent, father farmer
## orphan96: orphan in 1996
## hh_fthr_frm: father farmer
## hh_size96: household size in 1996

covs <- full_covs #c(AC_vars, "fthr_ed", "mthr_ed", "orphan96", "hh_fthr_frm", "hh_size96")
length(covs)
# covs <- c(AC_vars, "fthr_ed", "mthr_ed", "no_fthr96", "no_mthr96", "orphan96", "hh_fthr_frm", "hh_size96")

dim(df)

# remove rows with missingness in Z or Y
df <- df %>% filter(!is.na(Z) & !is.na(Y))

dim(df)

table(df$Z)

############################################################################################################
####################################### REGRESSION ON RESIDUALS ############################################
############################################################################################################
# regress Y on all AC_vars and take residuals as outcome
# lm_formula <- as.formula(
#   paste("Y ~", paste(AC_vars, collapse = " + "))
# )
# 
# fit <- lm(lm_formula, data = df)
# 
# resids <- resid(fit)
# length(resids)
# 
# df$Y <- resids
# covs <- covs[!covs %in% AC_vars]
############################################################################################################
############################################################################################################
############################################################################################################


naive.dim <- mean(df$Y[df$Z == 1]) - mean(df$Y[df$Z == 0])
naive.dim


data.c <- df %>% filter(Z == 0)
data.t <- df %>% filter(Z == 1)


n_repeat <- 10

full.dat <- bind_rows(data.c, data.t)

raw.all <- balancingWeights(data = full.dat, true_att = 0, feat_rep = "raw", 
                            raw_covs = covs, kernel_covs = NULL, verbose = TRUE)

raw.all

out_filename <- paste0("logs/soldiering-xfit-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

write("", out_filename, append = FALSE)   ##### ADDED (overwrite existing file)




single.fit <- function(pilot.dat, est.dat, treat.dat, covs) {
  analysis.dat <- bind_rows(treat.dat, est.dat)
  sing.fit.out <- eval_data(dat = analysis.dat, 
                            pilot.dat = pilot.dat, 
                            treat.true = 0, verbose = TRUE, covs = covs, dataset = "soldiering")
  sing.fit.out
}

process_finished_sim <- function(edat, id = 999) {
  out <- lapply(1:length(edat), function(i) {
    resi <- edat[[i]]
    elbo_rf <- resi$elbo_rf
    elbo_bart <- resi$elbo_bart
    rf_expl_var <- resi$rf_expl_var
    bart_expl_var <- resi$bart_expl_var
    expl_var_df <- data.frame(rf_expl_var = rf_expl_var,
                                 bart_expl_var = bart_expl_var)
    resi_rest <- resi[!names(resi) %in% c("elbo_rf", "elbo_bart", "rf_expl_var", "bart_expl_var", "kbal_expl_var")]
    results_df <- dplyr::bind_rows(resi_rest) %>% dplyr::mutate(elbo_rf = elbo_rf, elbo_bart = elbo_bart)
    l <- list(results_df = results_df, expl_var_df = expl_var_df)
    l
  })
  results_df <- dplyr::bind_rows(lapply(out, function(x) x$results_df)) %>% dplyr::mutate(id = id)
  expl_var_df <- dplyr::bind_rows(lapply(out, function(x) x$expl_var_df)) %>% dplyr::mutate(id = id)
  list(results_df = results_df, expl_var_df = expl_var_df)
  ###results_df <- dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
}

run_cross_fit <- function(data.c1, data.c2, treat.dat, covs, trans, id = 999) {
  fit.2.1.raw <- single.fit(pilot.dat = data.c2, est.dat = data.c1, treat.dat = treat.dat, covs = covs) 
  print("Finished first fit")
  fit.1.2.raw <- single.fit(pilot.dat = data.c1, est.dat = data.c2, treat.dat = treat.dat, covs = covs)
  print("Finished second fit") 
  
  
  # average each entry across the two fits. take entry-wise everage between fit.1.2[i,j] and fit.2.1[i,j]. 
  # in other words, give me (fit.1.2[i,j] + fit.2.1[i,j]) / 2
  # Assume fit.1.2 and fit.2.1 are already loaded and have the same structure
  
  fit.2.1 <- process_finished_sim(edat = fit.2.1.raw, id = id)$results_df
  fit.1.2 <- process_finished_sim(fit.1.2.raw, id = id)$results_df
  
  fit.2.1.expl_var <- process_finished_sim(fit.2.1.raw, id = id)$expl_var_df
  fit.1.2.expl_var <- process_finished_sim(fit.1.2.raw, id = id)$expl_var_df
  
  # Average the two data frames for expl_var
  expl_var.df <- fit.1.2.expl_var  # initialize with one of the data frames
  expl_var.df$rf_expl_var <- (fit.1.2.expl_var$rf_expl_var + fit.2.1.expl_var$rf_expl_var) / 2
  expl_var.df$bart_expl_var <- (fit.1.2.expl_var$bart_expl_var + fit.2.1.expl_var$bart_expl_var) / 2
  
  
  
  cfit.df <- fit.1.2  # initialize with one of the data frames
  
  # Identify numeric columns
  num_cols <- sapply(fit.1.2, is.numeric)
  
  # Average numeric columns
  cfit.df[, num_cols] <- (fit.1.2[, num_cols] + fit.2.1[, num_cols]) / 2
  
  cfit.df$se <- sqrt(fit.1.2$se^2 / 2 + fit.2.1$se^2 / 2)
  
  cfit.df$lcl <- cfit.df$est.att - 1.96 * cfit.df$se
  cfit.df$hcl <- cfit.df$est.att + 1.96 * cfit.df$se
  
  cfit.df$trans <- trans
  
  return(list(cfit.df = cfit.df, expl_var.df = expl_var.df))
}


numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

set.seed(12)

log_message <- paste("Starting soldiering XFIT \n")
cat(log_message, file = out_filename, append = TRUE)

out <- parallel::mclapply(1:n_repeat, function(i) {
  log_message <- paste("Starting repeat number", i, "at", Sys.time(), "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))
  
  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]

  out.notrans <- run_cross_fit(data.c1 = data.c1, data.c2 = data.c2, treat.dat = data.t, covs = covs, trans = "none", id = i)

  log_message <- paste("Finished repeat number", i, "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  out_df <- out.notrans
  out_df
}, mc.set.seed = TRUE, mc.cores = numCores - 1)


# save(out, file = "results/soldiering-xfit-educ.RData")
save(out, file = "results/soldiering-xfit-distress.RData")


