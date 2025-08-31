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

setwd("~/Desktop/BalWeights/forest-kbal/wsc")

load("data/wsc-obs-study-data.RData")

setwd("~/Desktop/BalWeights/forest-kbal/wsc")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

data.obs$Y <- data.obs$mathPost

data.pre <- data.obs
data.pre$std_math <- scale(data.pre$mathPost)

covs <- c(
  "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
  "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
  "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
  "mathPre"
)


data.pre <- data.pre %>% dplyr::rename(Z = m.treat)
data.pre.log <- data.pre


data.pre.log <- data.pre.log %>%
  dplyr::mutate(across(
    where(is.numeric) & !any_of(c("Z","Y")),
    ~ if (n_distinct(.x, na.rm = TRUE) >= 3 && min(.x, na.rm = TRUE) >= 0) log1p(.x) else .x))

data.c <- data.pre %>% filter(Z == 0)
data.t <- data.pre %>% filter(Z == 1)

data.c.log <- data.pre.log %>% filter(Z == 0)
data.t.log <- data.pre.log %>% filter(Z == 1)

n_repeat <- 6

set.seed(12)


single.fit <- function(pilot.dat, est.dat, treat.dat) {
  analysis.dat <- bind_rows(treat.dat, est.dat)
  sing.fit.out <- eval_data(dat = analysis.dat, 
                            pilot.dat = pilot.dat, 
                            treat.true = 0, verbose = TRUE, simulation = FALSE)
}

process_finished_sim <- function(edat) {
  out <- lapply(1:length(edat), function(i) {
    resi <- edat[[i]]
    dplyr::bind_rows(resi)
  })
  dplyr::bind_rows(out) %>% dplyr::mutate(id = i)
}

run_cfit <- function(data.c1, data.c2, treat.dat, trans) {
  fit.2.1.raw <- single.fit(pilot.dat = data.c2, est.dat = data.c1, treat.dat = treat.dat) 
  print("Finished first fit")
  fit.1.2.raw <- single.fit(pilot.dat = data.c1, est.dat = data.c2, treat.dat = treat.dat)
  print("Finished second fit") 
  
  # average each entry across the two fits. take entry-wise everage between fit.1.2[i,j] and fit.2.1[i,j]. 
  # in other words, give me (fit.1.2[i,j] + fit.2.1[i,j]) / 2
  # Assume fit.1.2 and fit.2.1 are already loaded and have the same structure
  
  fit.2.1 <- process_finished_sim(fit.2.1.raw)
  fit.1.2 <- process_finished_sim(fit.1.2.raw)
  
  
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

out <- parallel::mclapply(1:n_repeat, function(i) {
  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))
  
  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]
  
  data.c1.log <- data.c.log[sample_split, ]
  data.c2.log <- data.c.log[-sample_split, ]
  
  out.notrans <- run_cfit(data.c1 = data.c1, data.c2 = data.c2, treat.dat = data.t, trans = "none")
  out.log <- run_cfit(data.c1 = data.c1.log, data.c2 = data.c2.log, treat.dat = data.t.log, trans = "log")
  
  
}, mc.set.seed = TRUE, mc.cores = numCores - 1)




cfit.df <- average_cross_fit(out)

save(cfit.df, file = "results/wsc-math-cfit.RData")
