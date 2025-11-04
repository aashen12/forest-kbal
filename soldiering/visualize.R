
library(tidyverse)

rm(list = ls())
setwd("~/Desktop/BalWeights/forest-kbal/soldiering")
load("results/soldiering-xfit-educ-2010.RData")
#load("results/soldiering-xfit-distress-ror.RData")
full.df <- do.call(rbind, lapply(out, function(x) x$cfit.df))

# full.df$repeat_num <- rep(1:5, rep(22, 5)) 
dim(full.df)

head(full.df)

K <- max(full.df$id)

results.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::group_by(est, feat_rep, nc, trans) %>% 
  dplyr::summarise(
    est.att   = mean(est.att),
    se.xfit   = sqrt(sum(se^2)) / K,
    elbo_rf   = mean(elbo_rf),
    elbo_bart = mean(elbo_bart),
    .groups   = "drop"
  ) %>% 
  dplyr::mutate(
    lcl = est.att - 1.96 * se.xfit,
    ucl = est.att + 1.96 * se.xfit
  ) %>% 
  dplyr::ungroup()

results.df$nc <- as.numeric(results.df$nc)
head(results.df)

text_size <- 20
y_lim <- c(0.5, 1.5)
outcome <- "math"

raw0_lines <- results.df %>%
  dplyr::filter(feat_rep == "raw") %>%
  dplyr::select(est, trans, raw0_est = est.att)

rf0_lines <- results.df %>% 
  dplyr::filter(est == "rf", feat_rep == "raw") %>% 
  dplyr::select(est, trans, feat_rep, rf0_est = est.att)

benchmark_lines <- data.frame(
  est = "bal.wgt",
  trans = c("none", "log"),
  bench_est = c(1.01, 1.01)
)

kbal.only.df <- results.df %>% 
  filter(est == "kbal.bw", !is.na(est.att)) %>% 
  mutate(est = "bal.wgt", feat_rep = "kebal_only")


df_plot <- results.df %>%
  dplyr::filter(est == "bal.wgt") %>%
  #bind_rows(kbal.only.df) %>%
  dplyr::left_join(raw0_lines, by = c("est", "trans")) %>%
  dplyr::left_join(benchmark_lines, by = c("est", "trans")) %>%
  dplyr::filter(!feat_rep %in% c("raw")) %>%
  dplyr::filter(nc <= 20) %>% 
  dplyr::mutate(
    feat_group = dplyr::case_when(
      feat_rep %in% c("kbal_only", "rf_only", "bart_only", "kbal_K", "rf_K", "bart_K") ~ "Kernel Only",
      feat_rep %in% c("kbal_plus", "rf_plus", "kebal_only",  "bart_plus", "kbal_K_plus", "rf_K_plus", "bart_K_plus") ~ "Kernel + Raw",
      TRUE ~ NA_character_
    ),
    feat_group = factor(feat_group, levels = c("Kernel Only", "Kernel + Raw")),
    family = dplyr::case_when(
      grepl("kbal", feat_rep) ~ "KBal",
      grepl("kebal", feat_rep) ~ "KeBal",
      grepl("rf",   feat_rep) ~ "RF",
      grepl("bart", feat_rep) ~ "BART",
      TRUE ~ NA_character_
    ),
    family = factor(family, levels = c("KBal", "KeBal", "RF", "BART")),
    trans = ifelse(trans == "log", "Transformed", "Untransformed")
  ) %>%
  dplyr::filter(!is.na(feat_group), !is.na(family), feat_group %in% c("Kernel + Raw", "Kernel Only"))

head(df_plot[df_plot$feat_rep == "rf_plus",])


full.plot = df_plot %>% 
  dplyr::filter(feat_group=="Kernel + Raw") %>% 
  ggplot(
    aes(
      x = factor(nc),
      y = est.att,
      color = family,
      shape = family,
      group = interaction(family, feat_group)
    )
  ) +
  geom_point(size = 6.5) +
  # geom_errorbar(
  #   aes(ymin = lcl, ymax = ucl),
  #   width = 0.9,
  #   position = position_dodge(width = 0.05)
  # ) +
  geom_line(linewidth = 0.8) +
  # Baselines (mapped to color so they appear in the legend)
  geom_hline(
    aes(yintercept = raw0_est, color = "Raw Covariates"),
    linewidth = 1.0,
    alpha = 0.65,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  # geom_hline(
  #   aes(yintercept = -0.75, color = "Exp. Benchmark"),
  #   inherit.aes = FALSE, color = "red"
  # ) +
  #facet_wrap(~ feat_group, scales = "free_y") +
  labs(
    y = "Estimated ATT",
    x = "Number of Principal Components",
    title = NULL,
    color = "Features",
    shape = "Features"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text = element_text(size = text_size - 3, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_manual(
    name = "Features",
    values = c(
      KBal             = "#d62728",
      KeBal = "#9467bd",
      RF               = "#1f77b4",
      BART             = "#ff7f0e",
      "Raw Covariates" = "black"
    ),
    labels = c(
      KBal             = "Design Kernel (L2)",
      KeBal             = "Design Kernel (Entropy)",
      RF               = "RF Kernel",
      BART             = "BART Kernel",
      
      "Raw Covariates" = "Raw Covariates"
    ),
    breaks = c("KBal", "KeBal", "RF", "BART", "Raw Covariates")
  ) +
  scale_shape_manual(values = c(KBal = 17, KeBal = 17, RF = 16, BART = 16)) +
  guides(
    shape = "none",
    color = guide_legend(
      override.aes = list(
        shape     = c(17, 16, 16, NA),
        linetype  = c("blank", "blank", "blank", "dashed"),
        linewidth = c(NA, NA, NA, 0.9)
      )
    )
  )


full.plot



ggsave(
  filename = "paper-figs/soldiering_full_plot.pdf",
  plot     = full.plot,
  device   = "pdf",      # base grDevices::pdf()
  width    = 10,          # double-column width
  height   = 6,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


text_size = 30; y_lim = c(0.5, 1.5)

# desired top-to-bottom story order:
order_levels <- c(
  "raw_none",
  "kbal_plus_none",
  "rf_plus_none",
  "bart_plus_none",
  "raw_log",
  "kbal_plus_log",
  "rf_plus_log",
  "bart_plus_log"
)


bstars.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::group_by(est, feat_rep, nc, trans) %>% 
  dplyr::summarise(
    est.att = mean(est.att),
    se.xfit = sqrt(sum(se^2)) / K,
    elbo_rf = mean(elbo_rf),
    elbo_bart = mean(elbo_bart),
    .groups = "drop" # This ensures the output is not a grouped data frame
  ) %>% 
  dplyr::mutate(
    lcl = est.att - 1.96 * se.xfit,
    ucl = est.att + 1.96 * se.xfit
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(est == "bal.wgt", feat_rep == "raw" | nc == 5) %>% 
  dplyr::mutate(
    feat_group = dplyr::case_when(
      feat_rep %in% c("kbal_only", "rf_only", "bart_only") ~ "Kernel Only",
      feat_rep %in% c("kbal_plus","rf_plus", "bart_plus") ~ "Kernel + Raw",
      TRUE ~ NA_character_
    ),
    feat_group = factor(feat_group, levels = c("Kernel Only", "Kernel + Raw")),
    family = dplyr::case_when(
      grepl("kbal", feat_rep) ~ "KBal",
      grepl("rf",   feat_rep) ~ "RF",
      grepl("bart", feat_rep) ~ "BART",
      grepl("raw", feat_rep) ~ "Raw",
      TRUE ~ NA_character_
    ),
    family = factor(family, levels = c("KBal","RF", "BART", "Raw"))
  )



tf.plot = bstars.df %>% 
  mutate(
    feat_rep = factor(
      feat_rep,
      levels = rev(c("raw",
                     "kbal_only","kbal_plus",
                     "bart_only","bart_plus",
                     "rf_only","rf_plus"))
    )
  ) %>% 
  dplyr::filter(!feat_rep %in% c("bart_only", "kbal_only", "rf_only")) %>% 
  ggplot(aes(x = est.att, y = feat_rep, color = family)) +
  geom_point(size = 10) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.4, linewidth = 1.5) +
  theme_minimal() +
  #geom_vline(xintercept = 0.79, linetype = "dashed", linewidth = 1, color = "black") + 
  scale_color_manual(
    name   = "Features",
    values = c(
      KBal             = "#d62728",
      RF               = "#1f77b4",
      BART             = "#ff7f0e",
      Raw              = "gray33",
      "Exp. Benchmark" = "firebrick2"
    ),
    # labels = c(
    #   KBal = "Design Kernel",
    #   RF   = "RF Kernel",
    #   BART = "BART Kernel",
    #   Raw  = "Raw Covariates",
    #   "Exp. Benchmark" = "Exp. Benchmark"
    # ),
    breaks = c("KBal","RF","BART","Raw","Exp. Benchmark")
  ) +
  scale_y_discrete(
    labels = c(
      rf_plus    = "RF + Raw",
      rf_only    = "RF Only",
      kbal_plus  = "KBal + Raw",
      kbal_only  = "KBal Only",
      bart_plus  = "BART + Raw",
      bart_only  = "BART Only",
      raw        = "Raw"
    )
  ) +
  xlim(c(-1, -0.35)) + 
  labs(x = "Estimated ATT", y = "",
       color = "Features") +
  scale_shape_manual(values = c(KBal = 17, RF = 16, BART = 16, Raw = 15)) + 
  theme(
    text = element_text(size = text_size),
    strip.text = element_text(size = text_size - 7, face = "bold"),
    legend.position = "none"
  )

tf.plot 

ggsave(
  filename = "paper-figs/soldiering_tf_plot.pdf",
  plot     = tf.plot,
  device   = "pdf",      # base grDevices::pdf()
  width    = 9,          # double-column width
  height   = 6.5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)

scree.df <- do.call(rbind, lapply(out, function(x) x$expl_var.df))

num_per_id <- scree.df %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(n = n())

length_pc <- unique(num_per_id$n)

scree.df$nc <- rep(1:length_pc, nrow(scree.df) / length_pc)


scree.df.summ <- scree.df %>% 
  dplyr::group_by(id, nc) %>% 
  dplyr::summarise(
    rf_cumvar   = mean(rf_expl_var),
    bart_cumvar = mean(bart_expl_var),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = c(rf_cumvar, bart_cumvar),
    names_to = "kernel",
    values_to = "meanvar"
  )

# create spaghetti scree plot with points and lines running through. x-axis is nc,
# y axis is explained variance. one id per color
scree.plot = scree.df.summ %>% 
  filter(nc <= 15) %>% 
  ggplot(aes(x = nc, y = meanvar, color = factor(id))) +
  geom_point(size = 3, alpha = 0.6) +
  geom_line(linewidth = 1) +
  facet_wrap(~ kernel, ncol = 1,
             labeller = as_labeller(c(rf_cumvar = "RF Kernel",
                                      bart_cumvar = "BART Kernel"))) +
  labs(x = "Principal Component",
       y = "Explained Variance") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size = 18))

scree.plot

ggsave(
  filename = "paper-figs/soldiering_scree_plot.pdf",
  plot     = scree.plot,
  device   = "pdf",      # base grDevices::pdf()
  width    = 8,          # double-column width
  height   = 6,        # balanced height
  units    = "in",
  useDingbats = FALSE
)

unique(scree.df.summ$nc)

### PBR/ESS plots ###

metrics.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::group_by(est, feat_rep, nc, trans) %>% 
  dplyr::summarise(
    pbr = mean(pbr),
    ess = mean(ess),
    elbo_rf = mean(elbo_rf),
    elbo_bart = mean(elbo_bart),
    .groups = "drop" # This ensures the output is not a grouped data frame
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(est == "bal.wgt") %>% 
  dplyr::mutate(encoding = paste(feat_rep, trans, sep = "_")) %>%
  dplyr::mutate(
    feat_group = dplyr::case_when(
      encoding %in% c("kbal_only_log", "rf_only_log", "bart_only_log") ~ "Kernel Only (transformed)",
      encoding %in% c("kbal_plus_log","rf_plus_log", "bart_plus_log") ~ "Kernel + Raw (transformed)",
      encoding %in% c("kbal_only_none", "rf_only_none", "bart_only_none") ~ "Kernel Only",
      encoding %in% c("kbal_plus_none","rf_plus_none", "bart_plus_none") ~ "Kernel + Raw",
      TRUE ~ NA_character_
    ),
    feat_group = factor(feat_group, levels = c("Kernel Only", "Kernel + Raw", "Kernel Only (transformed)", "Kernel + Raw (transformed)")),
    family = dplyr::case_when(
      grepl("kbal", feat_rep) ~ "KBal",
      grepl("rf",   feat_rep) ~ "RF",
      grepl("bart", feat_rep) ~ "BART",
      grepl("raw", feat_rep) ~ "Raw",
      TRUE ~ NA_character_
    ),
    family = factor(family, levels = c("KBal","RF", "BART", "Raw"))
  )

metrics.df$nc <- as.numeric(metrics.df$nc)
metrics.df$nc <- factor(metrics.df$nc)

metrics.df <- metrics.df %>% mutate(trans = ifelse(trans == "none", "No Transformation", "Exponential Transformation"))
metrics.df$trans <- factor(metrics.df$trans, levels = c("No Transformation", "Exponential Transformation"))


table(metrics.df$nc)

head(metrics.df)



brks   <- c("KBal","RF","BART","Raw")
lbls   <- c(KBal="Design Kernel", RF="RF Kernel",
            BART="BART Kernel",   Raw="Raw Covariates")

raw0_metrics_df <- metrics.df %>%
  filter(feat_rep %in% c("raw")) %>% 
  dplyr::select(est, raw0_pbr = pbr, raw0_ess = ess)

metrics.df %>%
  filter(!feat_rep %in% c("raw")) %>%
  left_join(raw0_metrics_df, by = "est") %>% 
  filter(nc == 5) %>% data.frame()


metrics.df %>%
  filter(!feat_rep %in% c("raw")) %>%
  left_join(raw0_metrics_df, by = "est") %>%
  ggplot(aes(x = nc, y = ess,
             color = family, shape = family,
             linetype = family,                       # map linetype too
             group = interaction(family, feat_group))) +
  
  geom_point(size = 5.6) +
  geom_line(linewidth = 0.8) +
  
  geom_hline(aes(yintercept = raw0_ess, color = family, linetype = family),
    linewidth = 0.9, inherit.aes = FALSE
  ) +
  
  theme_bw() +
  facet_wrap(~ feat_group) +
  
  # --- unified legend: same name, breaks, labels across scales ---
  scale_color_manual(name = "Feature Family",
                     values = c(KBal="#d62728", RF="#1f77b4",
                                BART="#ff7f0e", Raw="black"),
                     breaks = brks, labels = lbls) +
  scale_shape_manual(name = "Feature Family",
                     values = c(KBal=17, RF=16, BART=16, Raw=NA),
                     breaks = brks, labels = lbls) +
  scale_linetype_manual(name = "Feature Family",
                        values = c(KBal="solid", RF="solid",
                                   BART="solid", Raw="dashed"),
                        breaks = brks, labels = lbls) +
  guides(
    # single combined legend appearance
    color = guide_legend(override.aes = list(
      shape    = c(17,16,16,NA),
      linetype = c("solid","solid","solid","dashed"),
      linewidth= c(0.8,0.8,0.8,0.9)
    )),
    shape = "none",      # don't make separate shape-only legend
    linetype = "none"    # don't make separate linetype-only legend
  ) + theme(text = element_text(size = 18))


metrics.df %>%
  filter(!feat_rep %in% c("bart_only","kbal_only","rf_only","rf_K"), nc == 5)

metrics.df %>%
  filter(feat_rep ==  "raw")



