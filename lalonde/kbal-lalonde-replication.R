# replicating chad's lalonde results
rm(list = ls())
library(tidyverse)
library(kbal)


data(lalonde, package = "kbal")

df.raw <- lalonde %>% dplyr::select(-race_ethnicity, -u78) %>% rename(treat = nsw, Y = re78)
head(df.raw)

dim(df.raw)

table(df.raw$treat)

df.log <- df.raw
df.log$re74 <- log1p(df.log$re74)
df.log$re75 <- log1p(df.log$re75)

df.sq <- df.raw
df.sq$re74_sq <- df.sq$re74^2
df.sq$re75_sq <- df.sq$re75^2
df.sq$age_sq <- df.sq$age^2


# run OLS to benchmark

ols.stnd <- lm(Y ~ treat + age + educ + black + hisp + married + nodegr + re74 + re75, data = df.raw)
summary(ols.stnd)

ols.log <- lm(Y ~ treat + age + educ + black + hisp + married + nodegr + re74 + re75, data = df.log)
summary(ols.log)

ols.sq <- lm(Y ~ treat + age + educ + black + hisp + married + nodegr + re74 + re75 + re74_sq + re75_sq + age_sq, data = df.sq)
summary(ols.sq)


df <- df.sq

X <- scale(df %>% dplyr::select(-treat, -Y))

kbal_obj <- kbal::kbal(X, treatment = df$treat, printprogress = TRUE, b = 10,
                       mixed_data = TRUE, cat_columns = c("black", "hisp", "married", "u74", "u75", "nodegr", "educ"))

w <- kbal_obj$w

length(w)
sum(w)



df$w <- ifelse(df$treat == 1, 1, w)

head(df)

# compute att
att <- with(df, mean(Y[treat == 1]) - weighted.mean(Y[treat == 0], w = w[treat == 0]))
att




