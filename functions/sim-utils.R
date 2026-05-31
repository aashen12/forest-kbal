# =============================================================================
# Shared Simulation Utilities
# =============================================================================
#
# Common helper functions used across simulation scripts.
# Eliminates duplication between simulation1.R, simulation2.R, and varsweep.R.
#
# Contents:
#   load_sim_packages()      - Load all required packages
#   find_elbow(y)            - Scree plot elbow detection
#   create_log(prefix)       - Create a timestamped log file
#   log_progress(id, ...)    - Log simulation progress at key intervals
#   process_eval_results()   - Flatten eval_data() output into a data frame
#   setup_parallel()         - Configure parallel workers from SLURM env
# =============================================================================


# --- Package Loading --------------------------------------------------------

#' Load all packages required for forest-kbal simulations.
#' Call this once at the top of each simulation script.
load_sim_packages <- function() {
  suppressPackageStartupMessages({
    library(MASS)
    library(Matrix)
    library(tidyverse)
    library(randomForest)
    library(stochtree)
    library(balancer)
    library(geepack)
    library(glmnet)
    library(kbal)
    library(sandwich)
    library(kernlab)
    library(rsample)
    library(irlba)
    library(parallel)
    library(future)
    library(future.apply)
  })
  options(dplyr.summarise.inform = FALSE)
}


# --- Elbow Detection --------------------------------------------------------

#' Find the elbow point of a scree plot using perpendicular distance.
#'
#' Draws a line from the first to the last eigenvalue, then returns
#' the index of the point with maximum perpendicular distance to that line.
#'
#' Used by both rf_kernel_matrix() and pca_bart() for selecting the
#' number of principal components.
#'
#' @param y Vector of eigenvalues (decreasing order)
#' @return Index of the elbow point
find_elbow <- function(y) {
  n <- length(y)
  x <- seq_len(n)

  # Endpoints of the line

  x1 <- x[1]; y1 <- y[1]
  x2 <- x[n]; y2 <- y[n]

  # Perpendicular distance from each point to the line
  num   <- abs((y2 - y1) * (x - x1) - (x2 - x1) * (y - y1))
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)

  which.max(num / denom)
}


# --- Logging ----------------------------------------------------------------

#' Create a timestamped log file.
#'
#' @param prefix Log file prefix (e.g., "sim1", "sim2-30", "varsweep")
#' @param log_dir Directory for log files (default: "logs")
#' @return Path to the created log file
create_log <- function(prefix, log_dir = "logs") {
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  path <- file.path(log_dir, paste0(prefix, "-", format(Sys.time(), "%b-%d-%H%M-%Y"), ".txt"))
  write("", path, append = FALSE)
  path
}

#' Log simulation progress at key intervals.
#'
#' Writes a progress message for the first 25 reps, the last 25 reps,
#' and every 25th rep in between.
#'
#' @param id Current simulation rep number
#' @param sim_reps Total number of reps
#' @param out_filename Path to log file
log_progress <- function(id, sim_reps, out_filename) {
  if (id <= 25 || id > (sim_reps - 25) || id %% 25 == 0) {
    cat(paste("Starting simulation", id, "at", Sys.time(), "\n"),
        file = out_filename, append = TRUE)
  }
}


# --- Result Processing ------------------------------------------------------

#' Flatten eval_data() output into a tidy data frame for one simulation rep.
#'
#' Takes the nested list returned by eval_data() or eval_data_varsweep()
#' and collapses it into a single data frame, extracting scalar fields
#' (like elbow indices) as columns.
#'
#' @param edat List returned by eval_data() or eval_data_varsweep()
#' @param id Simulation rep identifier (added as a column)
#' @param extra_fields Character vector of scalar fields to extract and keep
#'   as columns (default: elbo_rf, elbo_bart)
#' @param drop_fields Character vector of fields to remove before row-binding
#'   (default: elbo_rf, elbo_bart, rf_expl_var, bart_expl_var, kbal_expl_var)
#' @return A single data frame with all estimator results for this rep
process_eval_results <- function(edat, id,
                                  extra_fields = c("elbo_rf", "elbo_bart"),
                                  drop_fields = c("elbo_rf", "elbo_bart",
                                                   "rf_expl_var", "bart_expl_var",
                                                   "kbal_expl_var")) {
  out <- lapply(seq_along(edat), function(i) {
    resi <- edat[[i]]
    # Extract scalar fields (e.g., elbow indices)
    extras <- lapply(extra_fields, function(f) resi[[f]])
    names(extras) <- extra_fields
    # Remove non-data-frame fields before binding
    resi_rest <- resi[!names(resi) %in% drop_fields]
    row <- dplyr::bind_rows(resi_rest)
    # Add scalar fields as columns
    for (f in names(extras)) row[[f]] <- extras[[f]]
    row
  })

  out_df <- dplyr::bind_rows(out)
  out_df$id <- id
  rownames(out_df) <- NULL
  out_df
}


# --- Parallel Setup ---------------------------------------------------------

#' Configure parallel workers from SLURM environment.
#'
#' Reads SLURM_CPUS_PER_TASK to determine available cores.
#' Falls back to 4 cores if not running under SLURM (e.g., local testing).
#'
#' @param reserve Number of cores to reserve for the main process (default: 1)
#' @return Total number of cores available (before reserving)
setup_parallel <- function(reserve = 1) {
  num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
  plan(multisession, workers = num_cores)
  num_cores
}
