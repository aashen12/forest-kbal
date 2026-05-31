# =============================================================================
# Simulation 1: Forest Kernel Balancing — Tarr & Imai (2025) DGP
# =============================================================================
#
# Main simulation from Section 4 of Shen et al. (2025).
# Compares ATT estimation across feature representations (raw, RF kernel,
# BART kernel, Gaussian kernel) with multiple estimators (balancing weights,
# IPW, OLS, augmented/doubly robust).
#
# DGP: 10 nonlinearly-transformed covariates, logistic treatment assignment,
#      heterogeneous treatment effects. n=1000 per rep.
#      Pilot sample: 1000 additional units, controls only (~500).
#
# Usage:
#   Rscript simulation1.R <num_reps>
#
# Example:
#   sbatch sim1.sh 1000
#
# Arguments:
#   num_reps  Number of Monte Carlo replications
#
# Output:
#   results/simulation1-<num_reps>-<date>.RData
# =============================================================================

rm(list = ls())

# --- Setup -------------------------------------------------------------------
# Source shared utilities first (provides load_sim_packages, find_elbow, etc.)
source("../functions/sim-utils.R")
load_sim_packages()

source("../functions/dgp.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")

# --- Parse command-line arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
sim_reps <- as.numeric(args[1])

# --- Logging -----------------------------------------------------------------
out_filename <- create_log("sim1")
cat(paste("Number of reps:", sim_reps, "\n"), file = out_filename, append = TRUE)

# --- Parallel setup ----------------------------------------------------------
num_cores <- setup_parallel()

# --- Run simulation ----------------------------------------------------------
set.seed(23967)

cat(paste("Starting Simulation 1 at", Sys.time(), "\n"), file = out_filename, append = TRUE)

scenarios <- dplyr::bind_rows(mclapply(1:sim_reps, function(id) {
  log_progress(id, sim_reps, out_filename)

  # Generate analysis sample and pilot sample (controls only)
  bdat.obj <- make_data_sim1(1000)
  bdat <- bdat.obj$out.df
  pilot.dat <- make_data_sim1(1000)$out.df %>% dplyr::filter(Z == 0)

  # Run all estimators across all feature representations
  edat <- eval_data(dat = bdat, pilot.dat = pilot.dat,
                    treat.true = bdat.obj$true.att, verbose = FALSE,
                    dataset = "simulation")

  process_eval_results(edat, id)
}, mc.set.seed = TRUE, mc.cores = num_cores - 1))

# --- Save results ------------------------------------------------------------
save(scenarios, file = "simulation1-temp.RData")
cat("temp file saved\n", file = out_filename, append = TRUE)

filename <- paste0("results/simulation1-", sim_reps, "-", format(Sys.time(), "%b-%d-%Y"), ".RData")
save(scenarios, file = filename)
cat("saved file\n", file = out_filename, append = TRUE)
