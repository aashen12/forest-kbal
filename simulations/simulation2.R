# =============================================================================
# Simulation 2: Forest Kernel Balancing — Kim et al. (2024) DGP
# =============================================================================
#
# Appendix simulation from Shen et al. (2025).
# Uses Kim et al. (2024) DGP with varying overlap conditions.
# True ATT = 0 by construction.
#
# DGP: 6 covariates (3 correlated normals, 1 uniform, 1 chi-squared,
#      1 Bernoulli). Overlap controlled by sig.ep noise parameter.
#      n=1000 per rep. Pilot sample: controls only (~500).
#
# Usage:
#   Rscript simulation2.R <overlap> <num_reps>
#
# Example:
#   sbatch sim2.sh 30 1000    # moderate overlap
#   sbatch sim2.sh 100 1000   # strong overlap
#
# Arguments:
#   overlap   Noise parameter for treatment assignment (30 or 100)
#   num_reps  Number of Monte Carlo replications
#
# Output:
#   results/simulation2-overlap-<overlap>-numsim-<num_reps>-<date>.RData
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

# --- Parse command-line arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
overlap <- as.numeric(args[1])
sim_reps <- as.numeric(args[2])

# --- Logging -----------------------------------------------------------------
out_filename <- create_log(paste0("sim2-", overlap))
cat(paste("Number of reps:", sim_reps, "| Overlap:", overlap, "\n"),
    file = out_filename, append = TRUE)

# --- Parallel setup ----------------------------------------------------------
num_cores <- setup_parallel()

# --- Run simulation ----------------------------------------------------------
set.seed(23967)

cat(paste("Starting Simulation 2 at", Sys.time(), "\n"), file = out_filename, append = TRUE)

scenarios <- dplyr::bind_rows(mclapply(1:sim_reps, function(id) {
  log_progress(id, sim_reps, out_filename)

  # Generate analysis sample and pilot sample (controls only)
  bdat.obj <- make_data_sim2(n = 1000, sig.ep = overlap)
  bdat <- bdat.obj$out.df
  pilot.dat <- make_data_sim2(1000, sig.ep = overlap)$out.df %>% dplyr::filter(Z == 0)

  edat <- tryCatch({
    eval_data(dat = bdat, pilot.dat = pilot.dat,
              treat.true = bdat.obj$true.att, verbose = FALSE,
              dataset = "simulation")
  }, error = function(e) {
    cat(paste("eval_data failed on rep", id, ":", e$message, "\n"),
        file = out_filename, append = TRUE)
    save(list("dat" = bdat, "pilot.dat" = pilot.dat),
         file = paste0("eval-data-fail-", overlap, ".RData"))
    return(NULL)
  })

  if (is.null(edat)) return(NULL)
  process_eval_results(edat, id)
}, mc.set.seed = TRUE, mc.cores = num_cores - 1))

# --- Save results ------------------------------------------------------------
save(scenarios, file = "simulations2-temp.RData")
cat("temp file saved\n", file = out_filename, append = TRUE)

filename <- paste0("results/simulation2-overlap-", overlap, "-numsim-", sim_reps, "-",
                   format(Sys.time(), "%b-%d-%Y"), ".RData")
save(scenarios, file = filename)
cat("saved file\n", file = out_filename, append = TRUE)
