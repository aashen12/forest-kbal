# --- Packages ---
library(dplyr)
library(ggplot2)

rm(list = ls())

setwd("~/Desktop/BalWeights/forest-kbal/wsc")

transform <- "log"


if (transform == "exp") {
  load("results/wsc-math-xfit-exp.RData")
} else if (transform == "log") {
  load("results/wsc-math-xfit.RData")
}




full.df.trans <- do.call(rbind, lapply(out, function(x) x$out.log$cfit.df))
full.df.notrans <- do.call(rbind, lapply(out, function(x) x$out.notrans$cfit.df))
full.df <- rbind(full.df.trans, full.df.notrans) 

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
  dplyr::ungroup()

results.df$nc <- as.numeric(results.df$nc)

head(results.df)

full.df %>% group_by(trans) %>% 
  summarise(elbo_rf = mean(elbo_rf),
            elbo_bart = mean(elbo_bart))


text_size = 18; y_lim = c(0.5, 1.5)
outcome <- "math"

raw0_lines <- results.df %>%
  dplyr::filter(feat_rep == "raw") %>%
  dplyr::select(est, trans, raw0_est = est.att)

rf0_lines <- results.df %>% 
  filter(est == "rf", feat_rep == "raw") %>% 
  dplyr::select(est, trans, feat_rep, rf0_est = est.att)

bench_val <- if (outcome == "math") 0.79 else 2.18

kbal.only.df <- results.df %>% 
  filter(est == "kbal.bw", !is.na(est.att)) %>% 
  mutate(est = "bal.wgt", feat_rep = "kebal_only")

df_plot <- results.df %>%
  dplyr::filter(est %in% c("bal.wgt")) %>%
  #dplyr::bind_rows(kbal.only.df) %>%
  dplyr::left_join(raw0_lines, by = c("est", "trans")) %>%
  dplyr::filter(!feat_rep %in% c("raw","rf_K")) %>%
  dplyr::filter(nc <= 20) %>% 
  dplyr::mutate(
    feat_group = dplyr::case_when(
      feat_rep %in% c("kbal_only", "rf_only", "bart_only") ~ "Kernel Only",
      feat_rep %in% c("kbal_plus","rf_plus", "bart_plus", "kebal_only") ~ "Kernel + Raw",
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
  dplyr::filter(!is.na(feat_group), !is.na(family), feat_group == "Kernel + Raw") %>% 
  dplyr::mutate(trans = ifelse(trans == "log", ifelse(transform == "log", "Logarithmic", "Exponential"), "Untransformed"))

df_plot$trans <- factor(df_plot$trans, levels = c("Untransformed", ifelse(transform == "log", "Logarithmic", "Exponential")))


num_frs <- length(unique(df_plot$feat_rep))

p.full <- df_plot %>% 
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
             linewidth = 0.7,
             inherit.aes = FALSE) +
  facet_wrap(~ trans, scales = "free_y") +
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
      KeBal = "#e377c2",
      RF               = "#1f77b4",
      BART = "#ff7f0e",
      "Raw Covariates" = "black",
      "Exp. Benchmark" = "#2A873D"
    ),
    labels = c(
      KBal             = "Design Kernel",
      KeBal = "Entropy KBal",
      RF               = "RF Kernel",
      BART = "BART Kernel",
      "Raw Covariates" = "Raw Covariates",
      "Exp. Benchmark" = "Exp. Benchmark"
    ),
    breaks = c("KBal", "KeBal", "RF", "BART", "Raw Covariates","Exp. Benchmark")
  ) +
  scale_shape_manual(values = c(KBal = 17, Kebal = 17, RF = 16, BART = 16)) +
  guides(
    shape = "none",  # hide duplicate shape legend
    color = guide_legend(
      override.aes = list(
        shape    = c(17, 16, 16, NA, NA),              # points for KBal/RF; none for lines
        linetype = c("blank", "blank", "blank", "dashed", "solid"),
        linewidth= c(NA, NA, NA, 0.9, 0.9)
      )
    )
  )


p.full

ggsave(
  filename = paste0("paper-figs/wsc_full_", transform, ".pdf"),
  plot     = p.full,
  device   = "pdf",      # base grDevices::pdf()
  width    = 11,          # double-column width
  height   = 6,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


# full.df$est.att[full.df$feat_rep == "raw_0"] <- 1.008314
# full.df$se[full.df$feat_rep == "raw_0"] <- 0.3080238

full.df %>% filter(trans == "log", feat_rep == "raw_0", est == "bal.wgt")

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
  dplyr::mutate(encoding = paste(feat_rep, trans, sep = "_")) %>%
  dplyr::mutate(
    feat_group = dplyr::case_when(
      encoding %in% c("kbal_only_log", "rf_only_log", "bart_only_log") ~ "Kernel Only (transformed)",
      encoding %in% c("kbal_plus_log","rf_plus_log", "bart_plus_log") ~ "Kernel + Raw (transformed)",
      encoding %in% c("kbal_only_none", "rf_only_none", "bart_only_none") ~ "Kernel Only (untransformed)",
      encoding %in% c("kbal_plus_none","rf_plus_none", "bart_plus_none") ~ "Kernel + Raw (untransformed)",
      TRUE ~ NA_character_
    ),
    feat_group = factor(feat_group, levels = c("Kernel Only (untransformed)", "Kernel + Raw (untransformed)", "Kernel Only (transformed)", "Kernel + Raw (transformed)")),
    family = dplyr::case_when(
      grepl("kbal", feat_rep) ~ "KBal",
      grepl("rf",   feat_rep) ~ "RF",
      grepl("bart", feat_rep) ~ "BART",
      grepl("raw", feat_rep) ~ "Raw",
      TRUE ~ NA_character_
    ),
    family = factor(family, levels = c("KBal","RF", "BART", "Raw"))
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

bench_val <- 0.79
sd_val    <- 0.28
z_val     <- qnorm(0.975)

ci_low  <- bench_val - z_val * sd_val
ci_high <- bench_val + z_val * sd_val

shade_df <- data.frame(
  xmin = c(-Inf, ci_high),
  xmax = c(ci_low, Inf),
  ymin = -Inf,
  ymax = Inf,
  region = "Outside 95% CI"
)


bstars.df %>%
  mutate(
    encoding = factor(encoding, levels = rev(order_levels))  # put first items at top
  ) %>%
  filter(!feat_rep %in% c("bart_only","kbal_only","rf_only")) %>%
  ggplot(aes(x = est.att, y = encoding, color = family)) +
  # geom_rect(data = shade_df,
  #           aes(xmin = xmin, xmax = xmax,
  #               ymin = ymin, ymax = ymax,
  #               fill = region),
  #           inherit.aes = FALSE, alpha = 0.5) +
  geom_vline(xintercept = bench_val, linetype = "dashed",
             linewidth = 1, color = "black") + 
  scale_fill_manual(
    name = "Benchmark CI",
    values = c("Outside 95% CI" = "grey50")
  ) +
  geom_point(size = 6.5) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.55, linewidth = 0.9) +
  #geom_vline(xintercept = 0.79, linetype = "dashed", linewidth = 1, color = "black") +
  theme_minimal() +
  scale_color_manual(
    name   = "Features",
    values = c(KBal="#d62728", RF="#1f77b4", BART="#ff7f0e", Raw="gray33",
               "Exp. Benchmark"="firebrick2"),
    labels = c(KBal="Design Kernel", RF="RF Kernel", BART="BART Kernel",
               Raw="Raw Covariates", "Exp. Benchmark"="Exp. Benchmark"),
    breaks = c("KBal","RF","BART","Raw","Exp. Benchmark")
  ) +
  scale_y_discrete(labels = c(
    raw_none       = "Raw (untransformed)",
    kbal_plus_none = "KBal + Raw (untransformed)",
    rf_plus_none   = "RF + Raw (untransformed)",
    bart_plus_none = "BART + Raw (untransformed)",
    raw_log        = "Raw (transformed)",
    kbal_plus_log  = "KBal + Raw (transformed)",
    rf_plus_log    = "RF + Raw (transformed)",
    bart_plus_log  = "BART + Raw (transformed)"
  )) +
  labs(x = "Estimated ATT", y = "Feature Representation", color = "Features") +
  scale_shape_manual(values = c(KBal=17, RF=16, BART=16, Raw=15)) +
  theme(text = element_text(size = 23), legend.position = "right") #+xlim(0.4, 1.5)

ggsave(
  filename = "paper-figs/wsc_main.pdf",
  device   = "pdf",      # base grDevices::pdf()
  width    = 11,          # double-column width
  height   = 5.5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


bstars.df %>%
  dplyr::filter(trans == "log") %>% 
  mutate(
    encoding = factor(encoding, levels = rev(order_levels))  # put first items at top
  ) %>%
  filter(!feat_rep %in% c("bart_only","kbal_only","rf_only")) %>%
  ggplot(aes(x = est.att, y = encoding, color = family)) +
  # geom_rect(data = shade_df,
  #           aes(xmin = xmin, xmax = xmax,
  #               ymin = ymin, ymax = ymax,
  #               fill = region),
  #           inherit.aes = FALSE, alpha = 0.5) +
  geom_vline(xintercept = bench_val, linetype = "dashed",
             linewidth = 1, color = "black") + 
  scale_fill_manual(
    name = "Benchmark CI",
    values = c("Outside 95% CI" = "grey50")
  ) +
  geom_point(size = 6.5) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.55, linewidth = 0.9) +
  #geom_vline(xintercept = 0.79, linetype = "dashed", linewidth = 1, color = "black") +
  theme_minimal() +
  scale_color_manual(
    name   = "Features",
    values = c(KBal="#d62728", RF="#1f77b4", BART="#ff7f0e", Raw="gray33",
               "Exp. Benchmark"="firebrick2"),
    labels = c(KBal="Design Kernel", RF="RF Kernel", BART="BART Kernel",
               Raw="Raw Covariates", "Exp. Benchmark"="Exp. Benchmark"),
    breaks = c("KBal","RF","BART","Raw","Exp. Benchmark")
  ) +
  scale_y_discrete(labels = c(
    raw_none       = "Raw (untransformed)",
    kbal_plus_none = "KBal + Raw (untransformed)",
    rf_plus_none   = "RF + Raw (untransformed)",
    bart_plus_none = "BART + Raw (untransformed)",
    raw_log        = "Raw",
    kbal_plus_log  = "KBal + Raw",
    rf_plus_log    = "RF + Raw",
    bart_plus_log  = "BART + Raw"
  )) +
  labs(x = "Estimated ATT", y = "Feature Representation", color = "Features") +
  scale_shape_manual(values = c(KBal=17, RF=16, BART=16, Raw=15)) +
  theme(text = element_text(size = 23), legend.position = "right") #+xlim(0.4, 1.5)

ggsave(
  filename = "paper-figs/wsc_main_transformed.pdf",
  device   = "pdf",      # base grDevices::pdf()
  width    = 9,          # double-column width
  height   = 5.5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)




scree.df <- do.call(rbind, lapply(out, function(x) x$out.log$expl_var.df))

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
  labs(x = "Number of Principal Components",
       y = "Explained Variance") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size = 18))

scree.plot

ggsave(
  filename = "paper-figs/wsc_scree_plot.pdf",
  plot     = scree.plot,
  device   = "pdf",      # base grDevices::pdf()
  width    = 8,          # double-column width
  height   = 6,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


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
      encoding %in% c("kbal_only_none", "rf_only_none", "bart_only_none") ~ "Kernel Only (untransformed)",
      encoding %in% c("kbal_plus_none","rf_plus_none", "bart_plus_none") ~ "Kernel + Raw (untransformed)",
      TRUE ~ NA_character_
    ),
    feat_group = factor(feat_group, levels = c("Kernel Only (untransformed)", "Kernel + Raw (untransformed)", "Kernel Only (transformed)", "Kernel + Raw (transformed)")),
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

metrics.df %>%
  filter(!feat_rep %in% c("bart_only","kbal_only","rf_only","rf_K")) %>%
  ggplot(aes(x = nc, y = pbr,
             color = family, shape = family,
             linetype = family,                       # map linetype too
             group = interaction(family, feat_group))) +
  
  geom_point(size = 5.6) +
  geom_line(linewidth = 0.8) +
  
  geom_hline(
    data = metrics.df %>% filter(feat_rep == "raw"),
    aes(yintercept = pbr, color = family, linetype = family),
    linewidth = 0.9, inherit.aes = FALSE
  ) +
  
  theme_bw() +
  facet_wrap(~ trans) +
  
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

  
  
  
  

