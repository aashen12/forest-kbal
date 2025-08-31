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


setwd("~/Desktop/BalWeights/forest-kbal/simulations")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

# library(doFuture)
library(foreach)
library(doParallel)
library(parallel)
library(furrr)

args <- commandArgs(trailingOnly = TRUE)

overlap <- as.numeric(args[1])
sim_reps <- as.numeric(args[2])


out_filename <- paste0("logs/sim2-", overlap, "-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

write("", out_filename, append = FALSE)   ##### ADDED (overwrite existing file)

log_message <- paste("Number of reps:", sim_reps, "\n")
cat(log_message, file = out_filename, append = TRUE)

make_data <- function(n, sig.ep) {
  # Generate ZZ variables from standard normal distribution
  #print(sig.ep)
  sig.123 <- diag(c(2,1,1))
  sig.123[1,2] <- 1; sig.123[1,3] <- -1; sig.123[2,3] <- -0.5;
  sig.123 <- Matrix::forceSymmetric(sig.123)
  beta_coef <- c(1,2,-2,-1,-0.5,1)
  
  X.123 <- as.matrix(mvrnorm(n, mu = rep(0,3), Sigma = sig.123))
  colnames(X.123) <- paste0("X", 1:3)
  X4 <- runif(n,-3,3)
  X5 <- rchisq(n,1)
  X6 <- rbinom(n,1,0.5)
  X <- cbind(X.123, X4, X5, X6)
  X1 <- X[,1]
  X2 <- X[,2]
  X3 <- X[,3]
  expression_value <- X1^2 + 2*X2^2 - 2*X3^2 - (X4 + 1)^3 - 0.5*log(X5 + 10) + X6 - 1.5
  Z <- ifelse(expression_value + rnorm(n,0,sig.ep) > 0, 1, 0)
  Y <- (X.123[,1] + X.123[,2] + X5)^2 + rnorm(n,0,1)
  
  out.df <- data.frame(X, Z = Z, Y = Y)
  
  return(list(out.df=out.df, true.att = 0))
}

numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

### Sim Run

set.seed(23967)


run_scenario = function() {
  
  log_message <- paste("Starting Simulation 2 at", Sys.time(), "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  
  # Run the Simulation              
  reps_qs0 = mclapply( 1:sim_reps, function( id ) {
    if (id < 25) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    if (id > 975) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE) 
    if (id %% 25 == 0) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    
    bdat.obj <- make_data(n = 1000, sig.ep = overlap)
    bdat <- bdat.obj$out.df
    
    pilot.dat.obj <- make_data(n = 1000, sig.ep = overlap) # 1000 typically
    pilot.dat <- pilot.dat.obj$out.df %>% dplyr::filter(Z == 0)
    
    true.att.bdat <- bdat.obj$true.att
    
    edat <- tryCatch({
      eval_data(dat = bdat, 
                pilot.dat = pilot.dat, 
                treat.true = true.att.bdat, verbose = FALSE)
    },
    error = function(e) {
      print("Eval data failed")
      print(e)
      l <- list(dat = bdat, pilot.dat = pilot.dat, true.att.bdat = true.att.bdat)
      save(l, file = paste0("eval-data-fail-", overlap, ".RData"))
    })
    
    #edat$id = id
    #edat
    out <- lapply(1:length(edat), function(i) {
      resi <- edat[[i]]
      dplyr::bind_rows(resi)
    })
    dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
  }, mc.set.seed = TRUE, mc.cores = numCores - 1) 
  #cat("Sim Done")
  dplyr::bind_rows(reps_qs0)
}


scenarios <- run_scenario()

log_message <- "temp file saved \n"
cat(log_message, file = out_filename, append = TRUE)

save(scenarios, file="simulations2-temp.RData")

filename <- paste0("results/simulation2-overlap-", overlap, "-", "numsim-", sim_reps, "-", format(Sys.time(), "%b-%d-%Y"), ".RData")

save(scenarios, file=filename)
print("saved file")

log_message <- "saved file \n"
cat(log_message, file = out_filename, append = TRUE)


