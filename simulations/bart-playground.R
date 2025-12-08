rm(list=ls())
# generate Luke's data from Simulation 1 from Bal-Wgt-Educ
library(MASS)
library(tidyverse)
options(list(dplyr.summarise.inform = FALSE))
library(balancer)
library(geepack)
library(GenericML)
library(WeightIt)
library(glmnet)
library(kbal)
library(Matrix)
library(stochtree)
library(rsample)
library(kernlab)
library(randomForest)

setwd("~/Desktop/BalWeights/forest-kbal/simulations")
source("../functions/randomForestFeatures.R")
source("../functions/BART-features.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

library(foreach)
library(doParallel)
library(parallel)
library(furrr)

seed = 123


make_data <- function(n) {
  # Generate ZZ variables from standard normal distribution
  ZZ <- MASS::mvrnorm(n, mu = rep(0, 10), Sigma = diag(10))
  colnames(ZZ) <- paste0("ZZ", 1:10)
  
  # Define X variables using non-linear transformations of ZZ
  X1 <- exp(ZZ[, 1] / 2)
  X2 <- ZZ[, 2] / (1 + exp(ZZ[, 1]))
  X3 <- (ZZ[, 1] * ZZ[, 3] / 25 + 0.6)^3
  X4 <- (ZZ[, 2] + ZZ[, 4] + 20)^2
  X5 <- ZZ[, 5]
  X6 <- ZZ[, 6]
  X7 <- ZZ[, 7]
  X8 <- ZZ[, 8]
  X9 <- ZZ[, 9]
  X10 <- ZZ[, 10]
  
  # Assemble data frame for covariates
  X <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  
  # Define eta_3 using Weierstrass function form
  eta3 <- function(ZZ, a = 2, b = 13, N = 20) {
    Z_tilde <- (ZZ[,2] + ZZ[,4] + ZZ[,6] + ZZ[,8] + ZZ[,10]) / 5
    sapply(Z_tilde, function(x) {
      sum(sapply(0:N, function(n) a^(-n) * cos(b^n * pi * x)))
    })
  }
  
  
  # Compute treatment assignment probabilities
  logits <- -ZZ[, 1] - 0.1 * ZZ[, 4] #+ eta3(ZZ)
  prob_Z <- exp(logits) / (1 + exp(logits))
  Z <- rbinom(n, 1, prob_Z)
  
  
  # Compute outcome variable Y
  # Assumes ZZ is an n x 10 matrix, and T is a binary vector of length n
  
  interaction_term <- 27.4 * ZZ[,1] + 13.7 * ZZ[,2] + 13.7 * ZZ[,3] + 13.7 * ZZ[,4]
  
  Y0 <- 200 - 0.5 * interaction_term + rnorm(n)
  Y1 <- Y0 + 10 + 1.5 * interaction_term
  Y <- Z*Y1 + (1 - Z)*Y0
  
  
  # Return the full dataset
  out.df <- data.frame(X, Z = Z, Y = Y)
  true.att <- mean(Y1[Z == 1]) - mean(Y0[Z == 1])
  
  return(list(out.df=out.df, true.att = true.att))
}


bdat.obj <- make_data(1000)
bdat <- bdat.obj$out.df

pilot.dat.obj <- make_data(1000) # 1000 typically
pilot.dat <- pilot.dat.obj$out.df

bart_kernel <- bart_kernel_matrix(train = pilot.dat, test = bdat, verbose = TRUE, simulation = TRUE)

K_generic <- bart_kernel$kernel

bart.obj <- pca_bart(kernel = K_generic, data = bdat, n_components = 16)

