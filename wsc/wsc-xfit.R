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

table(data.pre$Z)


data.c <- data.pre %>% filter(Z == 0)
data.t <- data.pre %>% filter(Z == 1)


n_repeat <- 6

set.seed(12)


single.fit <- function(pilot.dat, est.dat) {
  analysis.dat <- bind_rows(data.t, est.dat)
  sing.fit.out <- eval_data(dat = analysis.dat, 
                            pilot.dat = pilot.dat, 
                            treat.true = 0, 
                            verbose = TRUE, simulation = FALSE)
}

numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

out <- parallel::mclapply(1:n_repeat, function(i) {
  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))
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
},mc.set.seed = TRUE, mc.cores = numCores - 1)


cfit.df <- average_cross_fit(out)

save(cfit.df, file = "results/wsc-math-cfit.RData")
