library(tidyverse)

rm(list=ls())

setwd("~/Desktop/BalWeights/forest-kbal/lalonde")


library(Matching)
data(lalonde, package = "Matching")
data.obs <- lalonde 
data.obs <- data.obs %>% dplyr::rename(Z = treat, Y = re78)
data.obs <- data.obs %>% dplyr::select(-u74, -u75)


table(data.obs$Z)

head(data.obs)

df <- data.obs
df$logY <- log1p(df$Y)  # log(1 + Y) to handle zeros
# DiM
att <- mean(df$Y[df$Z==1]) - mean(df$Y[df$Z==0])
att

log.att <- mean(df$logY[df$Z==1]) - mean(df$logY[df$Z==0])
log.att


hist(df$logY)
