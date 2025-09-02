# --- Packages ---
library(dplyr)
library(ggplot2)


load("results/wsc-math-xfit.RData")

out 

full.df <- do.call(rbind, out)

# full.df$repeat_num <- rep(1:5, rep(22, 5)) 

K <- max(full.df$id)



results.df <- full.df %>% 
  dplyr::mutate(
    nc = stringr::str_replace_all(feat_rep, "\\D+", ""),
    feat_rep = stringr::str_replace_all(feat_rep, "\\d+", "") %>% stringr::str_remove("_$")
  ) %>% 
  dplyr::relocate(nc, .after = "feat_rep") %>% 
  dplyr::group_by(est, feat_rep, nc, trans) %>% 
  dplyr::summarise(
    est.att = mean(est.att),
    se.xfit = sqrt(sum(se^2)) / K,
    elbo = mean(elbo),
    .groups = "drop" # This ensures the output is not a grouped data frame
  ) %>% 
  dplyr::mutate(
    lcl = est.att - 1.96 * se.xfit,
    ucl = est.att + 1.96 * se.xfit
  ) %>% 
  dplyr::ungroup()

results.df$nc <- as.numeric(results.df$nc)

head(results.df)

summary(full.df$elbo)

create_plot_wsc <- function(trans_level = c("none","log"),
                              text_size = 25,
                              y_lim = c(0.5, 1.5)) {
  trans_level <- match.arg(trans_level)
  outcome <- "math"

  raw0_lines <- results.df %>%
    dplyr::filter(feat_rep == "raw") %>%
    dplyr::select(est, trans, raw0_est = est.att)
  
  rf0_lines <- results.df %>% 
    filter(est == "rf", feat_rep == "raw") %>% 
    dplyr::select(est, trans, feat_rep, rf0_est = est.att)
  
  bench_val <- if (outcome == "math") 0.79 else 2.18
  
  results.df %>%
    dplyr::filter(trans == trans_level, est == "bal.wgt") %>%
    dplyr::left_join(raw0_lines, by = c("est", "trans")) %>%
    dplyr::filter(!feat_rep %in% c("raw","rf_K")) %>%
    dplyr::filter(nc <= 7 & nc >= 3) %>% 
    dplyr::mutate(
      feat_group = dplyr::case_when(
        feat_rep %in% c("kbal_only","rf_only") ~ "Kernel Only",
        feat_rep %in% c("kbal_plus","rf_plus") ~ "Kernel + Raw",
        TRUE ~ NA_character_
      ),
      feat_group = factor(feat_group, levels = c("Kernel Only","Kernel + Raw")),
      family = dplyr::case_when(
        grepl("kbal", feat_rep) ~ "KBal",
        grepl("rf",   feat_rep) ~ "RF",
        TRUE ~ NA_character_
      ),
      family = factor(family, levels = c("KBal","RF"))
    ) %>%
    dplyr::filter(!is.na(feat_group), !is.na(family)) %>%
    ggplot2::ggplot(aes(x = nc, y = est.att,
                        color = family, shape = family,
                        group = interaction(family, feat_group))) +
    geom_point(size = 7.6) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.9,
                  position = position_dodge(width = 0.05)) +
    geom_line(linewidth = 0.8) +
    # Baselines (mapped to color so they appear in the legend)
    geom_hline(aes(yintercept = raw0_est, color = "Raw Covariates"),
               linewidth = 1.0, alpha = 0.65, linetype = "dashed", inherit.aes = FALSE) +
    geom_hline(aes(yintercept = bench_val, color = "Exp. Benchmark"),
               inherit.aes = FALSE) +
    facet_wrap(~ feat_group, scales = "free_y") +
    coord_cartesian(ylim = y_lim) +
    labs(y = "Estimated ATT", x = "Number of Principal Components",
         title = NULL, color = "Features", shape = "Features") +
    theme_bw() +
    theme(text = element_text(size = text_size),
          strip.text = element_text(size = text_size - 3, face = "bold"),
          legend.position = "right") +
    scale_color_manual(
      name   = "Features",
      values = c(
        KBal             = "#d62728",
        RF               = "#1f77b4",
        "Raw Covariates" = "black",
        "Exp. Benchmark" = "firebrick2"
      ),
      labels = c(
        KBal             = "Design Kernel",
        RF               = "Forest Kernel",
        "Raw Covariates" = "Raw Covariates",
        "Exp. Benchmark" = "Exp. Benchmark"
      ),
      breaks = c("KBal","RF","Raw Covariates","Exp. Benchmark")
    ) +
    scale_shape_manual(values = c(KBal = 17, RF = 16)) +
    guides(
      shape = "none",  # hide duplicate shape legend
      color = guide_legend(
        override.aes = list(
          shape    = c(17, 16, NA, NA),              # points for KBal/RF; none for lines
          linetype = c("blank", "blank", "dashed", "solid"),
          linewidth= c(NA, NA, 0.9, 0.9)
        )
      )
    )
}

# create elbo plot



create_plot_wsc("none")
create_plot_wsc("log")





