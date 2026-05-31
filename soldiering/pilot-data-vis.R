# =============================================================================
# Soldiering Application: Pilot Data Visualization
# =============================================================================
#
# Visualizes results from pilot-data-analysis.R (no cross-fitting).
# Creates: line plot by nc, forest plot at nc=5.
#
# Input: results/multi-dataset-results-educ-dbldip.RData
# =============================================================================

library(tidyverse)
source("../functions/cross-fit.R")

load("results/multi-dataset-results-educ-dbldip.RData")

processed <- process_applied_results(edat, id = 1)
out_df <- processed$results_df
expl_var_df <- processed$expl_var_df

full.df <- out_df %>%
  dplyr::filter(est %in% c("bal.wgt"))

results.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::mutate(
    lcl = est.att - 1.96 * se,
    ucl = est.att + 1.96 * se
  ) %>% 
  dplyr::ungroup()

results.df$nc <- as.numeric(results.df$nc)
text_size <- 25

raw0_lines <- results.df %>%
  dplyr::filter(feat_rep == "raw") %>%
  dplyr::select(est, raw0_est = est.att)

df_plot <- results.df %>%
  dplyr::left_join(raw0_lines, by = c("est")) %>%
  dplyr::filter(!feat_rep %in% c("raw")) %>%
  dplyr::filter(nc <= 20) %>% 
  dplyr::mutate(
    feat_group = dplyr::case_when(
      feat_rep %in% c("kbal_only", "rf_only", "bart_only") ~ "Kernel Only",
      feat_rep %in% c("kbal_plus", "rf_plus", "kebal_only",  "bart_plus") ~ "Kernel + Raw",
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
    family = factor(family, levels = c("KBal", "KeBal", "RF", "BART"))
  ) %>%
  dplyr::filter(!is.na(feat_group), !is.na(family), feat_group == "Kernel + Raw")




df_plot %>% 
  ggplot2::ggplot(
    aes(
      x = nc,
      y = est.att,
      color = family,
      shape = family,
      group = interaction(family, feat_group)
    )
  ) +
  ggplot2::geom_point(size = 7.6) +
  ggplot2::geom_line(linewidth = 0.8) +
  # Baselines (mapped to color so they appear in the legend)
  ggplot2::geom_hline(
    aes(yintercept = raw0_est, color = "Raw Covariates"),
    linewidth = 1.0,
    alpha = 0.65,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggplot2::labs(
    y = "Estimated ATT",
    x = "Number of Principal Components",
    title = NULL,
    color = "Features",
    shape = "Features"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = text_size),
    strip.text = ggplot2::element_text(size = text_size - 3, face = "bold"),
    legend.position = "right"
  ) +
  ggplot2::scale_color_manual(
    name = "Features",
    values = c(
      KBal             = "#d62728",
      KeBal = "#9467bd",
      RF               = "#1f77b4",
      BART             = "#ff7f0e",
      "Raw Covariates" = "black",
      "Exp. Benchmark" = "firebrick2"
    ),
    labels = c(
      KBal             = "Design Kernel (L2)",
      KeBal             = "Design Kernel (Entropy)",
      RF               = "RF Kernel",
      BART             = "BART Kernel",
      
      "Raw Covariates" = "Raw Covariates",
      "Exp. Benchmark" = "Exp. Benchmark"
    ),
    breaks = c("KBal", "KeBal", "RF", "BART", "Raw Covariates", "Exp. Benchmark")
  ) +
  ggplot2::scale_shape_manual(values = c(KBal = 17, KeBal = 17, RF = 16, BART = 16)) +
  ggplot2::guides(
    shape = "none",
    color = ggplot2::guide_legend(
      override.aes = list(
        shape     = c(17, 16, 16, NA),
        linetype  = c("blank", "blank", "blank", "dashed"),
        linewidth = c(NA, NA, NA, 0.9)
      )
    )
  )



bstars.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::group_by(est, feat_rep, nc) %>% 
  dplyr::mutate(
    lcl = est.att - 1.96 * se,
    ucl = est.att + 1.96 * se
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
    labels = c(
      KBal = "Design Kernel",
      RF   = "RF Kernel",
      BART = "BART Kernel",
      Raw  = "Raw Covariates",
      "Exp. Benchmark" = "Exp. Benchmark"
    ),
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
  labs(x = "Estimated ATT", y = "Feature Representation",
       color = "Features") +
  scale_shape_manual(values = c(KBal = 17, RF = 16, BART = 16, Raw = 15)) + 
  theme(
    text = element_text(size = text_size),
    strip.text = element_text(size = text_size - 3, face = "bold"),
    legend.position = "right"
  )

tf.plot 


