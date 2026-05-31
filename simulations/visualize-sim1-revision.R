# =============================================================================
# Simulation 1 Revision: Dimensionality-Sweep Figures (Reviewer Comment 3)
# =============================================================================
#
# Reads the result files for each (q, n) cell listed in `results_files` below.
# You give only q and n; the matching combined file in results-spawn is found by
# pattern (the rep-count suffix is ignored), except the q = 10 / n = 1000 cell,
# which points at the original Simulation 1 output (no q/n column of its own).
# Each file is stamped with its q and n before binding.
#
# Output (paper-figs/revision/): each metric below produces two files, suffixed
# _kernelraw.pdf (Kernel + Raw) and _kernelonly.pdf (Kernel Only):
#   sim1revision_relrmse_*      Relative RMSE, balancing weights
#                               (kernelraw = Figure 1 main text; kernelonly = Figure 2 appendix)
#   sim1revision_relbias_*      Absolute relative bias, balancing weights
#   sim1revision_pbr_*          Percent bias reduction
#   sim1revision_ess_*          Effective sample size
#   sim1revision_relbias_ipw_*  Absolute relative bias, trimmed IPW
#   sim1revision_relrmse_ipw_*  Relative RMSE, trimmed IPW
#
# Usage:
#   Rscript visualize-sim1-revision.R
# =============================================================================

library(tidyverse)
source("../functions/sim-plot-utils.R")

# --- Inputs / outputs --------------------------------------------------------
out_dir <- "paper-figs/revision"

# Result files to plot -- one row per (q, n) cell of the sweep. `q` is the
# covariate dimension (facet columns) and `n` is the analysis sample size
# (facet rows).
#
# You only specify q and n; you do NOT name the file. Combined files are named
#   n  = 1000 :  sim1-revision-q-<q>-combined-<reps>.RData          (legacy)
#   n != 1000 :  sim1-revision-q-<q>-n-<n>-combined-<reps>.RData
# and the <reps> suffix changes every time you re-combine. So instead of
# hard-coding a filename, each (q, n) is resolved by pattern below with the rep
# count ignored -- the script keeps working as more reps land. The lone
# exception is q = 10 / n = 1000: the original Simulation 1 output, a
# differently named publication file with no q/n column of its own, so it gets
# an explicit path. Leave `file = NA` for anything that lives in results-spawn.
#
# infrastructure-job-arrays/ lives inside simulations/, so spawn_dir is a plain
# subdirectory path relative to this script (no leading "../").
spawn_dir <- "infrastructure-job-arrays/results-spawn"

results_files <- tibble::tribble(
  ~q,  ~n,   ~file,
  10,  1000, "results-revision/simulation1-1000-Nov-12-2025.RData",  # original Sim 1
  50,  1000, NA,
  100, 1000, NA,
  10,  2000, NA,
  50,  2000, NA,
  100, 2000, NA
)

#' Resolve one (q, n) cell to an actual combined-file path.
#'
#' A non-NA `file` is used as given. Otherwise glob `spawn_dir` for the (q, n)
#' naming pattern, ignoring the rep-count suffix. If several rep-count files
#' exist for the same (q, n) -- which happens because combining writes a new file
#' whenever the rep count changes, leaving the older ones behind -- the one with
#' the most reps is chosen.
resolve_result_file <- function(q, n, file, spawn_dir) {
  if (!is.na(file)) return(file)
  pat <- if (n == 1000) {
    sprintf("^sim1-revision-q-%d-combined-([0-9]+)\\.RData$", q)
  } else {
    sprintf("^sim1-revision-q-%d-n-%d-combined-([0-9]+)\\.RData$", q, n)
  }
  hits <- list.files(spawn_dir, pattern = pat)
  if (length(hits) == 0) {
    stop(sprintf("No combined file for q = %d, n = %d in %s (pattern: %s)",
                 q, n, spawn_dir, pat))
  }
  reps <- as.integer(sub(pat, "\\1", hits))
  file.path(spawn_dir, hits[which.max(reps)])
}

# Resolve every row to a concrete path, then report what was picked so it is
# obvious which rep count is being plotted.
results_files <- results_files %>%
  mutate(file = purrr::pmap_chr(list(q, n, file), function(q, n, file)
    resolve_result_file(q, n, file, spawn_dir)))

message("Result files in use:")
purrr::pwalk(results_files, function(q, n, file)
  message(sprintf("  q = %3d, n = %5d  ->  %s", q, n, basename(file))))

#' Load each resolved file and stamp it with its (q, n), then row-bind. Files may
#' carry their own q/n columns (the array output does); the stamp here is
#' authoritative so the legacy q = 10 file lines up too.
load_revision_results <- function(tbl) {
  missing <- !file.exists(tbl$file)
  if (any(missing)) stop("File(s) not found: ", paste(tbl$file[missing], collapse = ", "))
  purrr::pmap_dfr(tbl, function(q, n, file) {
    e <- new.env()
    load(file, envir = e)           # each file provides a `scenarios` object
    dplyr::mutate(e$scenarios, q = q, n = n)
  })
}

scenarios <- load_revision_results(results_files)
summ <- summarise_sim_results(scenarios)

# Figure size: q always spans the columns; when more than one sample size n is
# present, plot_sim_metric_byq() adds one facet row per n, so grow the height to
# keep each row from getting squished.
fig_w  <- 12
n_rows <- dplyr::n_distinct(scenarios$n)
fig_h  <- if (n_rows > 1) 3.6 * n_rows + 1.2 else 5

# --- Figures: one per (metric, feature group), faceted by n (rows) x q (cols) -
# Each plot_sim_metric_byq() call makes a grid with one colored line per kernel
# family per panel and a thin grey dashed raw-covariate baseline in each panel.
# The y-axis is shared across panels for comparison.

# Figure 1 (main text): Kernel + Raw, relative RMSE
fig_main <- plot_sim_metric_byq(summ, "rel_rmse", "Kernel + Raw", estimator = "bal.wgt")
fig_main
ggsave(file.path(out_dir, "sim1revision_relrmse_kernelraw.pdf"), fig_main,
       width = fig_w, height = fig_h, useDingbats = FALSE)

# Figure 2 (appendix): Kernel Only, relative RMSE
fig_appendix <- plot_sim_metric_byq(summ, "rel_rmse", "Kernel Only", estimator = "bal.wgt")
ggsave(file.path(out_dir, "sim1revision_relrmse_kernelonly.pdf"), fig_appendix,
       width = fig_w, height = fig_h, useDingbats = FALSE)

# Remaining appendix metrics, same structure -- a Kernel + Raw and a Kernel Only
# figure for each.
more_figs <- tibble::tribble(
  ~metric,    ~estimator, ~stem,
  "rel_bias", "bal.wgt",  "sim1revision_relbias",
  "pbr",      "bal.wgt",  "sim1revision_pbr",
  "ess",      "bal.wgt",  "sim1revision_ess",
  "rel_bias", "ipw",      "sim1revision_relbias_ipw",
  "rel_rmse", "ipw",      "sim1revision_relrmse_ipw"
)

for (i in seq_len(nrow(more_figs))) {
  m <- more_figs$metric[i]; est <- more_figs$estimator[i]; stem <- more_figs$stem[i]
  ggsave(file.path(out_dir, paste0(stem, "_kernelraw.pdf")),
         plot_sim_metric_byq(summ, m, "Kernel + Raw", estimator = est),
         width = fig_w, height = fig_h, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(stem, "_kernelonly.pdf")),
         plot_sim_metric_byq(summ, m, "Kernel Only", estimator = est),
         width = fig_w, height = fig_h, useDingbats = FALSE)
}
