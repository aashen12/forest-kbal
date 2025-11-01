rm(list=ls())
if (!require(rsample)) install.packages("rsample"); library(rsample)
if (!require(kernlab)) install.packages("kernlab"); library(kernlab)
if (!require(stochtree)) install.packages("stochtree"); library(stochtree)
library(MASS)
library(tidyverse)
options(list(dplyr.summarise.inform = FALSE))
library(balancer)
library(geepack)
library(GenericML)
library(WeightIt)
library(glmnet)
library(kbal)
library(randomForest)

setwd("~/Desktop/BalWeights/forest-kbal/simulations")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")


library(parallel)
library(future)
library(future.apply)

verbose = TRUE


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


n <- 1000

bdat.obj <- make_data(n)
dat <- bdat.obj$out.df

pilot.dat.obj <- make_data(n) # 1000 typically

# keep the first 500 rows
pilot.dat <- pilot.dat.obj$out.df %>% dplyr::filter(Z == 0)

true.att.bdat <- bdat.obj$true.att


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

covs <- names(dat)[!names(dat) %in% c("Z", "Y")]

#### Generate Random Forest features #####
nc_rf <- nc_bart <- c(2, 5, 10, 15, 25, 50, 100)

if (verbose) print("Starting RF Kernel")
# rf.scenarios.partial <- expand.grid(fr = c("rf_only", "rf_plus", "rf_mixed"), ncomp = nc_rf)
rf.scenarios.partial <- expand.grid(fr = c("rf_only", "rf_plus"), ncomp = nc_rf)
rf.scenarios.ker <- expand.grid(fr = c("rf_K"), ncomp = nrow(data))
rf.scenarios <- rbind(rf.scenarios.partial, rf.scenarios.ker)
## First, run random forest

form <- reformulate(covs, response = "Y")
rfmod <- randomForest(form, data = data.c.train, ntree = 100)

rf_obj <- rf_kernel_matrix(model = rfmod, data = data, X = as.matrix(data %>% dplyr::select(-Y)),
                           n_components = max(nc_rf), verbose = TRUE)
elbo_rf <- rf_obj$elbow

data_rf_only_whole <- rf_obj$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
data_rf_plus_whole <- rf_obj$data_rf
# data_rf_mixed_whole <- rf_obj$features_mixed %>% dplyr::mutate(Y = data$Y, Z = data$Z)
data_rf_K <- rf_obj$K 


### KBAL ###
kbal.scenarios <- expand.grid(fr = c("kbal_only", "kbal_plus"), ncomp = nc_rf)
X_unscaled <- model.matrix(reformulate(c(covs, "-1")), data)
# scale to have variance 1 but not mean 0
X <- scale(X_unscaled)
kbal_objs <- lapply(1:length(nc_rf), function(i) {
  nc <- nc_rf[i]
  sim2 <- length(unique(X[,6])) == 2 # X6 in sim2 is bernoulli, we just need a quick check for Sim1 or Sim2
  print(paste("Sim2:", sim2))
  kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,
                         mixed_data = sim2, cat_columns = ifelse(sim2, "X6", NULL))
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
if (verbose) print("KBal Finished!")




bart.scenarios.partial <- expand.grid(fr = c("bart_only", "bart_plus"), ncomp = nc_bart)
bart.scenarios <- bart.scenarios.partial
bart_obj <- bart_kernel_matrix(train = data.c.train, test = data, seed = 1022, verbose = TRUE, covs = covs)
bart.pca <- pca_bart(kernel = bart_obj$kernel, data = data, X = as.matrix(data %>% dplyr::select(-Y)), 
                     n_components = max(nc_bart))
elbo_bart <- bart.pca$elbow
data_bart_only_whole <- bart.pca$features %>% dplyr::mutate(Y = data$Y, Z = data$Z)
data_bart_plus_whole <- bart.pca$data_bart


kbal_pcs <- kbal_objs[[length(kbal_objs)]]$data_kbal_only %>% dplyr::select(-kbal_weights, -treat, -Y)
colnames(kbal_pcs) <- paste0("KBAL_PC", 1:ncol(kbal_pcs))

bart_pcs <- bart.pca$features
colnames(bart_pcs) <- paste0("BART_PC", 1:ncol(bart_pcs))

rf_pcs <- rf_obj$features
colnames(rf_pcs) <- paste0("RF_PC", 1:ncol(rf_pcs))

data <- cbind(bart_pcs, rf_pcs, kbal_pcs, Y = dat$Y, Z = dat$Z)

write_csv(data, "results/kernel_pcs.csv")
# make a ggplot of the first two PCs of RF & BART

andy_theme <- function() {
    theme_minimal() + 
    theme(text = element_text(size = 16))
}


pboth <- data %>% 
  ggplot(aes(x = BART_PC2, y = RF_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "blue") + 
  andy_theme() + 
  labs(x = "PC2 (BART)", y = "PC2 (RF)")
pboth


pbart = data %>% 
  ggplot(aes(x = BART_PC1, y = BART_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "#ff7f0e") + 
  andy_theme() + 
  labs(x = "PC1 (BART)", y = "PC2 (BART)")
pbart

ggsave(
  filename = "paper-figs/bart_pc.pdf",
  plot     = pbart,
  device   = "pdf",      # base grDevices::pdf()
  width    = 6,          # double-column width
  height   = 5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


prf = data %>% 
  ggplot(aes(x = RF_PC1, y = RF_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "#1f77b4") + 
  andy_theme() + 
  labs(x = "PC1 (RF)", y = "PC2 (RF)")
prf

# ggsave(
#   filename = "paper-figs/rf_pc.pdf",
#   plot     = prf,
#   device   = "pdf",      # base grDevices::pdf()
#   width    = 6,          # double-column width
#   height   = 5,        # balanced height
#   units    = "in",
#   useDingbats = FALSE
# )
# 








X <- as.matrix(data %>% dplyr::select(-Y, -Z))
kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = NULL, printprogress = TRUE,
                       mixed_data = FALSE, cat_columns = NULL)


Kkbal <- kbal_obj$svdK$u %*% diag(kbal_obj$svdK$d) %*% t(kbal_obj$svdK$v)
Kbart <- bart_obj$kernel
Kbart <- Kbart / unique(diag(Kbart))
Krf <- data_rf_K
dim(Krf)
dim(Kkbal)
dim(Kbart)
head(Kbart[, 1:8])
head(Krf[, 1:8])

Z = dat["Z"]


# sample an index of Z such that Z == 1
treat.ind <- which(dat$Z == 1)
target <- sample(treat.ind, size = 1)
target

source <- setdiff(1:nrow(dat), target)

# want to look at k(i, j) for all j in source, where i is the target index
bart_vals = Kbart[target, source]
rf_vals = Krf[target, source]
kbal_vals = Kkbal[target, source]
rf_vals = c(rf_vals) %>% unlist()

df <- data.frame(bart = bart_vals, rf = rf_vals, kbal = kbal_vals)
rownames(df) <- NULL
head(df)


write_csv(df, "results/kernel_comparison.csv")


df %>% 
  # dplyr::filter(rf > 0) %>% 
  ggplot(aes(x = bart, y = rf)) +
  geom_point(alpha = 0.6, color = "forestgreen", size = 2) +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "BART", y = "RF") + 
  andy_theme()



# Overlapping histograms of BART vs RF (assumes `df` in memory)

df_hist <- df %>%
  dplyr::select(bart, rf) %>%                     # only the two columns
  tidyr::pivot_longer(everything(),
                      names_to = "method",
                      values_to = "value") %>%
  dplyr::filter(is.finite(value))

# common binwidth for aligned histograms
rng <- range(df_hist$value)
binwidth <- diff(rng) / 40   # ~40 bins; adjust as needed

ggplot(df_hist, aes(x = value, fill = method)) +
  geom_histogram(
    position = "identity",
    alpha = 0.45,
    binwidth = binwidth,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(bart = "#1B9E77", rf = "#D95F02"),
    name = NULL,
    labels = c("BART", "RF")
  ) +
  labs(x = "Value", y = "Count") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 11)
  )


