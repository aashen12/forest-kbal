# =============================================================================
# Simulation 1: Paper Figures (Figure 3 and Figures B.1-B.4)
# =============================================================================
#
# Reproduces the Simulation 1 figures in Shen et al. (2025) from a single
# results file produced by simulation1.R. All plotting logic lives in
# ../functions/sim-plot-utils.R; this script just loads, summarizes, and saves.
#
# Figure map:
#   Figure 3  (main)  : Relative RMSE,        balancing weights -> sim1_relrmse.pdf
#   Figure B.1        : Absolute relative bias, balancing weights -> sim1_relbias.pdf
#   Figure B.2        : Percent bias reduction + effective sample size
#                       -> sim1_pbr.pdf, sim1_ess.pdf
#   Figure B.3        : Absolute relative bias, trimmed IPW -> sim1_relbias_ipw.pdf
#   Figure B.4        : Relative RMSE,          trimmed IPW -> sim1_relrmse_ipw.pdf
#
# Usage:
#   Rscript visualize-sim1.R
# (Edit `results_file` below to point at the run you want to plot.)
# =============================================================================

library(tidyverse)
source("../functions/sim-plot-utils.R")

# --- Inputs / outputs --------------------------------------------------------
results_file <- "results/simulation1-1000-Nov-12-2025.RData"
out_dir      <- "paper-figs"

load(results_file)              # provides `scenarios`
summ <- summarise_sim_results(scenarios)

# --- Figure 3 (main): relative RMSE, balancing weights -----------------------
fig3 <- plot_sim_metric(summ, "rel_rmse", estimator = "bal.wgt", ylim = c(0, 0.5))
ggsave(file.path(out_dir, "sim1_relrmse.pdf"), fig3,
       width = 10, height = 4.5, useDingbats = FALSE)

# --- Figure B.1: relative bias, balancing weights ----------------------------
b1 <- plot_sim_metric(summ, "rel_bias", estimator = "bal.wgt", ylim = c(0, 0.5))
ggsave(file.path(out_dir, "sim1_relbias.pdf"), b1,
       width = 10, height = 4.5, useDingbats = FALSE)

# --- Figure B.2: percent bias reduction and effective sample size ------------
pbr <- plot_sim_metric(summ, "pbr", estimator = "bal.wgt")
ggsave(file.path(out_dir, "sim1_pbr.pdf"), pbr,
       width = 9, height = 4.5, useDingbats = FALSE)

ess <- plot_sim_metric(summ, "ess", estimator = "bal.wgt")
ggsave(file.path(out_dir, "sim1_ess.pdf"), ess,
       width = 9, height = 4.5, useDingbats = FALSE)

# --- Figures B.3 / B.4: trimmed IPW instead of balancing weights -------------
b3 <- plot_sim_metric(summ, "rel_bias", estimator = "ipw", ylim = c(0, 2))
ggsave(file.path(out_dir, "sim1_relbias_ipw.pdf"), b3,
       width = 10, height = 4.5, useDingbats = FALSE)

b4 <- plot_sim_metric(summ, "rel_rmse", estimator = "ipw", ylim = c(0, 3))
ggsave(file.path(out_dir, "sim1_relrmse_ipw.pdf"), b4,
       width = 10, height = 4.5, useDingbats = FALSE)
