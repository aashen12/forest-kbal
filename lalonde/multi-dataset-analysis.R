library(tidyverse)
rm(list = ls())
load("results/multi-dataset-results-logged.RData")


out <- lapply(1:length(edat), function(i) {
  resi <- edat[[i]]
  elbo_rf <- resi$elbo_rf
  elbo_bart <- resi$elbo_bart
  resi_rest <- resi[!names(resi) %in% c("elbo_rf", "elbo_bart")]
  dplyr::bind_rows(resi_rest) %>% dplyr::mutate(elbo_rf = elbo_rf, elbo_bart = elbo_bart)
})

out_df <- dplyr::bind_rows(out)
rownames(out_df) <- NULL
out_df
unique(out_df$est)

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
head(results.df)

text_size <- 25
y_lim <- c(0.5, 1.5)
outcome <- "math"

raw0_lines <- results.df %>%
  dplyr::filter(feat_rep == "raw") %>%
  dplyr::select(est, raw0_est = est.att)

rf0_lines <- results.df %>% 
  dplyr::filter(est == "rf", feat_rep == "raw") %>% 
  dplyr::select(est, feat_rep, rf0_est = est.att)

benchmark_lines <- data.frame(
  est = "bal.wgt",
  trans = c("none", "log"),
  bench_est = c(1.01, 1.01)
)

kbal.only.df <- results.df %>% 
  filter(est == "kbal.bw", !is.na(est.att)) %>% 
  mutate(est = "bal.wgt", feat_rep = "kebal_only")


df_plot <- results.df %>%
  #bind_rows(kbal.only.df) %>%
  dplyr::left_join(raw0_lines, by = c("est")) %>%
  dplyr::filter(!feat_rep %in% c("raw", "rf_K")) %>%
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
  ggplot2::geom_errorbar(
    aes(ymin = lcl, ymax = ucl),
    width = 0.9,
    position = ggplot2::position_dodge(width = 0.05)
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  # Baselines (mapped to color so they appear in the legend)
  ggplot2::geom_hline(
    aes(yintercept = raw0_est, color = "Raw Covariates"),
    linewidth = 1.0,
    alpha = 0.65,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_hline(
    aes(yintercept = 1.01, color = "Exp. Benchmark"),
    inherit.aes = FALSE
  ) +
 # ggplot2::facet_wrap(~ trans, scales = "free_y") +
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
        shape     = c(17, 16, 16, NA, NA),
        linetype  = c("blank", "blank", "blank", "dashed", "solid"),
        linewidth = c(NA, NA, NA, 0.9, 0.9)
      )
    )
  )






  