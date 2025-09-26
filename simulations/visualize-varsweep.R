# visualize varsweep

library(tidyverse)

rm(list=ls())


ww <- 20 # plot width
hh <- 8 # plot height

load("results/varsweep-2000-Sep-24-2025.RData")

unique(scenarios$feat_rep)

#### Sim Summary
# number at the end of the string, after an underscore
num_tail <- "(?i)(?<=_)[+-]?(?:\\d*\\.\\d+|\\d+)(?:e[+-]?\\d+)?$"
rm_tail  <- "_(?i)[+-]?(?:\\d*\\.\\d+|\\d+)(?:e[+-]?\\d+)?$"

df <- scenarios %>%
  mutate(
    varsweep    = str_extract(feat_rep, num_tail),     # character or NA
    varsweep    = as.numeric(varsweep),                   # NA stays NA
    feat_rep = str_remove(feat_rep, rm_tail)        # drop the trailing _<number>
  )

df$varsweep
df$logit_p <- round(df$logit_p, 2)


unique(df$logit_p)

library(dplyr)

# make ordered factor: smallest -> largest
levs  <- sort(unique(df$logit_p))
labs  <- levs #ifelse(levs == -5, "NA", format(levs, scientific = TRUE, digits = 3))

df <- df %>%
  mutate(logit_p = factor(logit_p, levels = levs, labels = labs),
         varsweep = factor(varsweep, levels = sort(unique(varsweep)), labels = sort(unique(varsweep))))


df_plot <- df %>% 
  group_by(est, feat_rep, logit_p, varsweep) %>% 
  dplyr::summarise(
    avg_bias = mean(abs(bias)),
    avg_rel_bias = mean(abs(rel.bias)),
    sd_bias = sd(abs(bias)),
    sd_rel_bias = sd(abs(rel.bias)),
    mse = sqrt(mean(bias^2)),
    sd_mse = sd(sqrt(bias^2)),
    rel_mse = sqrt(mean(rel.bias^2)),
    sd_rel_mse = sd(sqrt(rel.bias^2)),
    avg_cvg = mean(cvg),
    avg_pbr = mean(pbr),
    avg_ess = mean(ess),
    .groups = "drop"
  ) %>% 
  mutate(method = recode(est, 
                         bw_l2 = "Bal. Weights - L2", 
                         bw_inf = "Bal. Weights - L-Inf",
                         lasso_inf = "Augmented Linf BW with Lasso",
                         ipw = "Logistic IPW",
                         ols = "OLS Outcome Regression",
                         ols_l2 = "Aug BW with OLS",
                         ols_ipw = "AIPW with OLS",
                         ridge_l2 = "Aug BW with Ridge",
                         bw.simp = "L2 BW (scaled)",
                         rf = "Random Forest Outcome",
                         bal.wgt = "Balancing Weights"),
         feature_rep = recode(feat_rep,
                              raw = "Raw Covs",
                              rf_only = "RF Feats",
                              rf_plus = "RF Feats + Raw Covs",
                              bart_only = "BART Feats",
                              bart_plus = "BART Feats + Raw Covs",
                              rf_mixed = "RF Mixed Feats",
                              bart_mixed = "BART Mixed Feats",
                              kbal_only = "KBal Feats",
                              kbal_plus = "KBal Feats + Raw Covs"))


df_plot %>% 
  filter(!feat_rep %in% c("raw_NA", "rf_only_NA")) %>% 
  ggplot(aes(x = varsweep, y = avg_rel_bias)) +
  geom_point(aes(color = "RF + Raw"), size = 8, alpha = 0.9) +
  theme_bw() +
  geom_hline(aes(yintercept = df_plot$avg_rel_bias[df_plot$feat_rep == "raw_NA"],
                 color = "Raw Covs", linetype = "Raw Covs")) +
  geom_hline(aes(yintercept = df_plot$avg_rel_bias[df_plot$feat_rep == "rf_only_NA"],
                 color = "RF Only", linetype = "RF Only")) +
  theme(text = element_text(size = 20)) +
  scale_color_manual(
    values = c("RF + Raw" = "seagreen", "Raw Covs" = "red", "RF Only" = "blue"),
    breaks = c("RF + Raw", "Raw Covs", "RF Only"),
    name = "Feature Set"
  ) +
  scale_linetype_manual(
    values = c("Raw Covs" = "solid", "RF Only" = "solid"),
    guide = "none"   # <-- remove the extra legend
  ) +
  labs(x = "Relative importance of kernel features",
       y = "Average Relative Bias")


