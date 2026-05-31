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
library(randomForest)
library(kbal)
library(geepack)
library(irlba)

rm(list=ls())

setwd("~/Desktop/BalWeights/forest-kbal/lalonde")
setwd("~/Desktop/BalWeights/forest-kbal/lalonde")

set.seed(2025-10-20)

source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

load("data/lalonde.RData")
data_pilot <- ldw_cps %>% dplyr::rename(Z = treat, Y = re78, hisp = hispanic, educ = education, nodegr = nodegree) %>% 
  dplyr::select(-sample, -data_id)
head(data_pilot)

pilot.dat <- data_pilot %>% dplyr::filter(Z == 0) #%>% sample_n(2000)
dim(pilot.dat)

pilot.dat.log <- pilot.dat
pilot.dat.log$Y <- log1p(pilot.dat.log$Y)



data(lalonde, package = "kbal")
data.as <- lalonde 
data.as <- data.as %>% dplyr::rename(Z = nsw, Y = re78) %>% dplyr::select(-race_ethnicity, -u78)

# data.as$Y <- log1p(data.as$Y)  # log(1 + Y) to handle zeros

table(data.as$Z)

dim(data.as)

covs <- names(data.as)[!names(data.as) %in% c("Z", "Y")]
length(covs)
covs

data.as.log <- data.as

naive.dim <- mean(data.as$Y[data.as$Z == 1]) - mean(data.as$Y[data.as$Z == 0])
naive.dim

log_trans <- TRUE

if (log_trans == TRUE) {
  # log transform
  # data.as.log$re74 <- log1p(data.as.log$re74)
  # data.as.log$re75 <- log1p(data.as.log$re75)
  # data.as.log$age <- log1p(data.as.log$age)
  data.as.log$Y <- log1p(data.as.log$Y)
} else {
  # exp transform
  data.as.log <- data.as.log %>%
    mutate(across(
      .cols = where(is.numeric) & !any_of(c("Z", "Y", "re74", "re75")),
      .fns = ~ if (n_distinct(.x) >= 2 && min(.x, na.rm = TRUE) >= 0) {
        exp(.x)
      } else {
        .x
      }
    ))
}

edat <- eval_data(dat = data.as.log, pilot.dat = pilot.dat.log, 
                  treat.true = 0, verbose = TRUE, simulation = FALSE)

# edat <- eval_data(dat = data.as, pilot.dat = pilot.dat, 
#                   treat.true = 0, verbose = TRUE, simulation = FALSE)

save(edat, file = "results/multi-dataset-results-logged.RData")

# out <- lapply(1:length(edat), function(i) {
#   resi <- edat[[i]]
#   elbo_rf <- resi$elbo_rf
#   elbo_bart <- resi$elbo_bart
#   resi_rest <- resi[!names(resi) %in% c("elbo_rf", "elbo_bart")]
#   dplyr::bind_rows(resi_rest) %>% dplyr::mutate(elbo_rf = elbo_rf, elbo_bart = elbo_bart)
# })
# 
# out_df <- dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
# rownames(out_df) <- NULL
# out_df
