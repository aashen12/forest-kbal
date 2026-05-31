# =============================================================================
# Varsweep: Sensitivity to Feature Weight Parameter
# =============================================================================
#
# Sweeps over the mixing parameter p that controls the relative weight
# between raw covariates and kernel features in the balancing objective:
#   X_combined = cbind((1-p) * X_raw, p * X_kernel)
# where p=0 is raw-only and p=1 is kernel-only.
#
# Uses the same DGP as Simulation 1 (Tarr & Imai, 2025).
# See Section 2.3 of Shen et al. (2025) for the normalization discussion.
#
# Usage:
#   Rscript varsweep.R <num_reps>
#
# Example:
#   sbatch varsweep.sh 2000
#
# Arguments:
#   num_reps  Number of Monte Carlo replications
#
# Output:
#   results/varsweep-<num_reps>-<date>.RData
# =============================================================================

rm(list = ls())

# --- Setup -------------------------------------------------------------------
source("../functions/sim-utils.R")
load_sim_packages()

source("../functions/dgp.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")
source("../functions/varsweep-funcs.R")

# --- Parse command-line arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
sim_reps <- as.numeric(args[1])

# --- Logging -----------------------------------------------------------------
out_filename <- create_log("varsweep")
cat(paste("Number of reps:", sim_reps, "\n"), file = out_filename, append = TRUE)

# --- Parallel setup ----------------------------------------------------------
num_cores <- setup_parallel()

# --- Run simulation ----------------------------------------------------------
set.seed(23967)

cat(paste("Starting Varsweep at", Sys.time(), "\n"), file = out_filename, append = TRUE)

scenarios <- dplyr::bind_rows(mclapply(1:sim_reps, function(id) {
  log_progress(id, sim_reps, out_filename)

  # Generate analysis sample and pilot sample (controls only)
  bdat.obj <- make_data_sim1(1000)
  bdat <- bdat.obj$out.df
  pilot.dat <- make_data_sim1(1000)$out.df %>% dplyr::filter(Z == 0)

  edat <- eval_data_varsweep(dat = bdat, pilot.dat = pilot.dat,
                              treat.true = bdat.obj$true.att, verbose = FALSE,
                              simulation = TRUE)

  process_eval_results(edat, id,
                       extra_fields = c("elbo_rf", "elbo_bart"),
                       drop_fields = c("elbo_rf", "elbo_bart"))
}, mc.set.seed = TRUE, mc.cores = num_cores - 1))

# --- Save results ------------------------------------------------------------
save(scenarios, file = "varsweep-temp.RData")
cat("temp file saved\n", file = out_filename, append = TRUE)

filename <- paste0("results/varsweep-", sim_reps, "-", format(Sys.time(), "%b-%d-%Y"), ".RData")
save(scenarios, file = filename)
cat("saved file\n", file = out_filename, append = TRUE)
