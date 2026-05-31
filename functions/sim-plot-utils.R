# =============================================================================
# Simulation Plotting Utilities
# =============================================================================
#
# Shared helpers for summarizing and plotting Simulation 1 results
# (Figure 3 and Figures B.1-B.4 of Shen et al. 2025), plus the Simulation 1
# revision that sweeps the covariate dimension q.
#
# Sourced by:
#   simulations/visualize-sim1.R           (original figures)
#   simulations/visualize-sim1-revision.R  (dimensionality sweep figures)
#
# Requires tidyverse (dplyr, ggplot2, stringr) to be attached by the caller.
#
# Contents:
#   SIM_METRIC_SPEC        - metric -> {summary column, axis label}
#   SIM_FAMILY_COLORS      - kernel family -> color
#   SIM_FAMILY_SHAPES      - kernel family -> point shape
#   summarise_sim_results  - collapse raw reps into per-cell means
#   plot_sim_metric        - one faceted metric plot (optionally overlaying q)
# =============================================================================


# --- Visual constants --------------------------------------------------------

# Kernel families. "Design Kernel" is the Gaussian design-based kernel (kbal);
# RF and BART are the forest kernels. Colors match the paper figures.
SIM_FAMILY_COLORS <- c(
  "Design Kernel"   = "#d62728",
  "RF"              = "#1f77b4",
  "BART"            = "#ff7f0e",
  "Raw Covariates"  = "black"
)
SIM_FAMILY_SHAPES <- c("Design Kernel" = 17, "RF" = 16, "BART" = 16)
SIM_FAMILY_LEVELS <- c("Design Kernel", "RF", "BART")

# Linetypes used to distinguish covariate dimension q when overlaying.
SIM_Q_LINETYPES <- c("solid", "longdash", "dotted", "dotdash", "twodash")

# Metric -> summary column and axis label. One entry per panel type.
SIM_METRIC_SPEC <- list(
  bias     = list(col = "avg_bias",     label = "Absolute Bias"),
  rel_bias = list(col = "avg_rel_bias", label = "Absolute Relative Bias"),
  rmse     = list(col = "rmse",         label = "Absolute RMSE"),
  rel_rmse = list(col = "rel_rmse",     label = "Relative RMSE"),
  pbr      = list(col = "avg_pbr",      label = "Percent Bias Reduction"),
  ess      = list(col = "avg_ess",      label = "Effective Sample Size")
)


# --- Summarize ---------------------------------------------------------------

#' Collapse raw per-replication results into per-cell means.
#'
#' Splits the trailing principal-component count out of `feat_rep` into a
#' numeric `ncomp`, relabels the full-kernel representations (e.g. "rf_K") as
#' their truncated counterparts so they appear as the largest-PC point, and
#' averages the per-rep metrics within each estimator x feature x PC cell. If a
#' `q` column is present (Simulation 1 revision), it is added to the grouping.
#'
#' @param scenarios Raw results data frame from a simulation run
#' @return One row per (est, feat_rep, ncomp[, q]) with averaged metrics
summarise_sim_results <- function(scenarios) {
  has_q <- "q" %in% names(scenarios)
  has_n <- "n" %in% names(scenarios)
  group_vars <- c("est", "feat_rep", "ncomp", if (has_q) "q", if (has_n) "n")

  scenarios %>%
    mutate(
      ncomp    = as.numeric(stringr::str_extract(feat_rep, "[0-9]+")),
      feat_rep = stringr::str_remove(feat_rep, "_[0-9]+"),
      # Full kernel (all components) becomes the largest-PC point of its family
      feat_rep = recode(feat_rep,
                        rf_K = "rf_only",     rf_K_plus = "rf_plus",
                        bart_K = "bart_only", bart_K_plus = "bart_plus",
                        kbal_K = "kbal_only", kbal_K_plus = "kbal_plus")
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      avg_bias     = mean(abs(bias),     na.rm = TRUE),
      avg_rel_bias = mean(abs(rel.bias), na.rm = TRUE),
      rmse         = sqrt(mean(bias^2,     na.rm = TRUE)),
      rel_rmse     = sqrt(mean(rel.bias^2, na.rm = TRUE)),
      avg_pbr      = mean(pbr, na.rm = TRUE),
      avg_ess      = mean(ess, na.rm = TRUE),
      .groups = "drop"
    )
}


# --- Plot --------------------------------------------------------------------

#' Plot one simulation metric across the number of principal components.
#'
#' Faceted by feature group (Kernel Only vs Kernel + Raw), with one line per
#' kernel family (Design Kernel, RF, BART) and a dashed reference line for
#' balancing on raw covariates alone. This reproduces the structure of
#' Figure 3 and Figures B.1-B.4.
#'
#' To show multiple covariate dimensions q (the revision), set
#' `q_display = "overlay"` (q encoded by linetype, all on one panel) or
#' `q_display = "facet"` (q as an extra facet column).
#'
#' @param summ Output of summarise_sim_results()
#' @param metric One of names(SIM_METRIC_SPEC)
#' @param estimator Weighting estimator to show ("bal.wgt" or "ipw")
#' @param q_display "none" (default), "overlay", or "facet"
#' @param ylim Optional y-axis limits via coord_cartesian (does not drop data)
#' @param text_size Base font size
#' @return A ggplot object
plot_sim_metric <- function(summ, metric, estimator = "bal.wgt",
                            q_display = c("none", "overlay", "facet"),
                            ylim = NULL, text_size = 18) {
  metric    <- match.arg(metric, names(SIM_METRIC_SPEC))
  q_display <- match.arg(q_display)
  spec      <- SIM_METRIC_SPEC[[metric]]
  y_col     <- rlang::sym(spec$col)

  if (q_display != "none" && !("q" %in% names(summ))) {
    stop("q_display = '", q_display, "' requires a 'q' column in `summ`.")
  }

  # Raw-covariate baseline (one value per estimator, or per q if overlaying)
  raw_ref <- summ %>%
    filter(feat_rep == "raw", est == estimator) %>%
    select(any_of("q"), raw0 = !!y_col)

  plot_df <- summ %>%
    filter(est == estimator, feat_rep != "raw") %>%
    mutate(
      feat_group = case_when(
        feat_rep %in% c("kbal_only", "rf_only", "bart_only") ~ "Kernel Only",
        feat_rep %in% c("kbal_plus", "rf_plus", "bart_plus") ~ "Kernel + Raw",
        TRUE ~ NA_character_
      ),
      feat_group = factor(feat_group, levels = c("Kernel Only", "Kernel + Raw")),
      family = case_when(
        grepl("kbal", feat_rep) ~ "Design Kernel",
        grepl("rf",   feat_rep) ~ "RF",
        grepl("bart", feat_rep) ~ "BART",
        TRUE ~ NA_character_
      ),
      family = factor(family, levels = SIM_FAMILY_LEVELS)
    ) %>%
    filter(!is.na(feat_group), !is.na(family))

  overlay <- identical(q_display, "overlay")
  if (overlay) {
    plot_df$q <- factor(plot_df$q)
    raw_ref$q <- factor(raw_ref$q)
  }

  # Base mapping: family always drives color + shape; q (if overlaid) drives linetype
  base_aes <- if (overlay) {
    aes(x = factor(ncomp), y = !!y_col, color = family, shape = family,
        linetype = q, group = interaction(family, feat_group, q))
  } else {
    aes(x = factor(ncomp), y = !!y_col, color = family, shape = family,
        group = interaction(family, feat_group))
  }

  p <- ggplot(plot_df, base_aes) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 3)

  # Raw-covariate reference line(s). geom_hline never inherits the main aes,
  # so it is safe to feed it a data frame without the x/family columns.
  if (overlay) {
    p <- p + geom_hline(data = raw_ref, aes(yintercept = raw0, linetype = q),
                        color = "black", linewidth = 0.5)
  } else {
    p <- p + geom_hline(data = raw_ref, aes(yintercept = raw0),
                        color = "black", linetype = "dashed", linewidth = 0.6)
  }

  # Faceting
  if (identical(q_display, "facet")) {
    p <- p + facet_grid(feat_group ~ q,
                        labeller = labeller(q = function(x) paste0("q = ", x)))
  } else {
    p <- p + facet_grid(~ feat_group)
  }

  # Scales, labels, theme
  p <- p +
    scale_color_manual(values = SIM_FAMILY_COLORS, breaks = SIM_FAMILY_LEVELS) +
    scale_shape_manual(values = SIM_FAMILY_SHAPES, breaks = SIM_FAMILY_LEVELS) +
    labs(x = "Number of principal components", y = spec$label,
         color = "Kernel", shape = "Kernel") +
    theme_bw(base_size = text_size)

  if (overlay) {
    n_q <- nlevels(plot_df$q)
    p <- p +
      scale_linetype_manual(values = SIM_Q_LINETYPES[seq_len(n_q)]) +
      labs(linetype = "Covariates (q)") +
      # color/shape duplicate each other -> single combined "Kernel" legend
      guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
             shape = guide_legend(order = 1),
             linetype = guide_legend(order = 2))
  } else {
    # color and shape duplicate -> merge into one legend
    p <- p + guides(color = guide_legend(), shape = guide_legend())
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  p
}


#' Plot one metric across PCs for a single feature group, faceted by q (and n).
#'
#' Used for the Simulation 1 revision. Produces one figure for a single feature
#' group (e.g. "Kernel + Raw" for the main text, "Kernel Only" for the appendix),
#' with the covariate dimension q as facet columns (left to right) and one
#' colored line per kernel family. When the summary spans more than one sample
#' size n, those become facet rows (one per n). A thin grey dashed line marks the
#' raw-covariate baseline within each facet (it differs across both q and n).
#' Unlike plot_sim_metric(q_display = "overlay"), q is NOT encoded by linetype,
#' so each panel shows only three lines.
#'
#' @param summ Output of summarise_sim_results() (must contain a `q` column)
#' @param metric One of names(SIM_METRIC_SPEC)
#' @param feat_group "Kernel + Raw" or "Kernel Only"
#' @param estimator Weighting estimator to show ("bal.wgt" or "ipw")
#' @param ylim Optional y-axis limits via coord_cartesian (shared across facets)
#' @param text_size Base font size
#' @return A ggplot object
plot_sim_metric_byq <- function(summ, metric,
                                feat_group = c("Kernel + Raw", "Kernel Only"),
                                estimator = "bal.wgt", ylim = NULL, text_size = 18) {
  metric     <- match.arg(metric, names(SIM_METRIC_SPEC))
  feat_group <- match.arg(feat_group)
  spec       <- SIM_METRIC_SPEC[[metric]]
  y_col      <- rlang::sym(spec$col)

  if (!("q" %in% names(summ))) {
    stop("plot_sim_metric_byq() requires a 'q' column in `summ`.")
  }

  # Order q so facets render left to right (e.g. 10, 50, 100)
  q_levels <- sort(unique(summ$q))

  # Facet by sample size n (one row per n) only when more than one n is present.
  # With a single n (e.g. the per-setting previews) the figure stays a single
  # row faceted by q, exactly as before.
  has_n    <- "n" %in% names(summ) && dplyr::n_distinct(summ$n) > 1
  n_levels <- if (has_n) sort(unique(summ$n)) else NULL

  # Raw-covariate baseline: one value per facet (per q, and per n when faceting
  # by n as well, since the baseline differs across both).
  raw_ref <- summ %>%
    filter(feat_rep == "raw", est == estimator) %>%
    select(any_of("n"), q, raw0 = !!y_col) %>%
    mutate(q = factor(q, levels = q_levels))
  if (has_n) raw_ref <- mutate(raw_ref, n = factor(n, levels = n_levels))

  plot_df <- summ %>%
    filter(est == estimator, feat_rep != "raw") %>%
    mutate(
      group = case_when(
        feat_rep %in% c("kbal_only", "rf_only", "bart_only") ~ "Kernel Only",
        feat_rep %in% c("kbal_plus", "rf_plus", "bart_plus") ~ "Kernel + Raw",
        TRUE ~ NA_character_
      ),
      family = case_when(
        grepl("kbal", feat_rep) ~ "Design Kernel",
        grepl("rf",   feat_rep) ~ "RF",
        grepl("bart", feat_rep) ~ "BART",
        TRUE ~ NA_character_
      ),
      family = factor(family, levels = SIM_FAMILY_LEVELS),
      q = factor(q, levels = q_levels)
    ) %>%
    filter(group == feat_group, !is.na(family))
  if (has_n) plot_df <- mutate(plot_df, n = factor(n, levels = n_levels))

  # Shared x axis. The standard PC grid tops out at 100; the only larger value
  # is the full-kernel rep (ncomp == n), which differs across sample sizes. Map
  # every full-kernel point onto a single shared "K" tick so the n = 1000 and
  # n = 2000 rows line up on one discrete axis instead of showing separate
  # "1000" / "2000" ticks crammed against each other.
  grid_levels <- sort(unique(plot_df$ncomp[plot_df$ncomp <= 100]))
  plot_df <- plot_df %>%
    mutate(x = factor(ifelse(ncomp > 100, "K", as.character(ncomp)),
                      levels = c(as.character(grid_levels), "K")))

  # Facet: q across columns always; n across rows when more than one is present.
  facet_layer <- if (has_n) {
    facet_grid(n ~ q, labeller = labeller(q = function(x) paste0("q = ", x),
                                          n = function(x) paste0("n = ", x)))
  } else {
    facet_grid(~ q, labeller = labeller(q = function(x) paste0("q = ", x)))
  }

  p <- ggplot(plot_df, aes(x = x, y = !!y_col,
                           color = family, shape = family,
                           group = interaction(family, q))) +
    # Raw-covariate reference per facet -- dark grey (not full black) and
    # slightly thicker than the data lines so it reads as the anchor
    geom_hline(data = raw_ref, aes(yintercept = raw0),
               color = "grey30", linetype = "dashed", linewidth = 0.8) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 3) +
    facet_layer +
    scale_color_manual(values = SIM_FAMILY_COLORS, breaks = SIM_FAMILY_LEVELS) +
    scale_shape_manual(values = SIM_FAMILY_SHAPES, breaks = SIM_FAMILY_LEVELS) +
    labs(x = "Number of principal components", y = spec$label,
         color = "Kernel", shape = "Kernel") +
    theme_bw(base_size = text_size) +
    guides(color = guide_legend(), shape = guide_legend())

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  p
}
