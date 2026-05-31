# =============================================================================
# Simulation 1 Revision: Effect of Covariate Dimensionality (Reviewer Comment 3)
# =============================================================================
#
# Revision simulation addressing the reviewer's question: how does the
# dimensionality q of the covariates X affect the performance of forest kernel
# balancing, and how large q can the method handle?
#
# This duplicates Simulation 1 (Tarr & Imai 2025 DGP) but varies the number of
# covariates q. The DGP cycles the nonlinearities AND the treatment/outcome
# model over consecutive blocks of 10 latent Gaussians, so added covariates are
# genuine confounders (4 active + 6 noise per block). The propensity/outcome
# signal is rescaled by sqrt(#blocks) so propensity overlap stays fixed as q
# grows -- isolating the pure dimensionality effect. See make_data_sim1_highdim().
#
# With q = 10 the DGP reduces exactly to Simulation 1, anchoring the sweep.
#
# DGP: q nonlinearly-transformed covariates, logistic treatment assignment,
#      heterogeneous treatment effects. n=1000 per rep.
#      Pilot sample: 1000 additional units, controls only (~500).
#
# Usage:
#   Rscript sim1-revision.R <q> <num_reps>
#
# Example:
#   sbatch sim1-revision.sh 10 1000
#   sbatch sim1-revision.sh 40 1000
#   sbatch sim1-revision.sh 100 1000
#
# Arguments:
#   q         Number of covariates (10, 40, or 100)
#   num_reps  Number of Monte Carlo replications
#
# Output:
#   results/sim1-revision-q-<q>-numsim-<num_reps>-<date>.RData
# =============================================================================
print("Starting script")
rm(list = ls())

# --- Setup -------------------------------------------------------------------
source("../functions/sim-utils.R")
load_sim_packages()
print("Packages loaded")

source("../functions/dgp.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")
print("Functions loaded")

# --- Parse command-line arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
q <- as.numeric(args[1])
sim_reps <- as.numeric(args[2])
print(q)
print(sim_reps)

# --- Logging -----------------------------------------------------------------
out_filename <- create_log(paste0("sim1-revision-q", q))
cat(paste("Number of reps:", sim_reps, "| q:", q, "\n"),
    file = out_filename, append = TRUE)

# --- Parallel setup ----------------------------------------------------------
num_cores <- setup_parallel()

# --- Run simulation ----------------------------------------------------------
set.seed(23967)

cat(paste("Starting Simulation 1 revision (q =", q, ") at", Sys.time(), "\n"),
    file = out_filename, append = TRUE)

scenarios <- dplyr::bind_rows(mclapply(1:sim_reps, function(id) {
  log_progress(id, sim_reps, out_filename)

  # Generate analysis sample and pilot sample (controls only)
  bdat.obj <- make_data_sim1_highdim(n = 1000, q = q)
  bdat <- bdat.obj$out.df
  pilot.dat <- make_data_sim1_highdim(n = 1000, q = q)$out.df %>% dplyr::filter(Z == 0)

  edat <- tryCatch({
    eval_data(dat = bdat, pilot.dat = pilot.dat,
              treat.true = bdat.obj$true.att, verbose = FALSE,
              dataset = "simulation")
  }, error = function(e) {
    cat(paste("eval_data failed on rep", id, ":", e$message, "\n"),
        file = out_filename, append = TRUE)
    return(NULL)
  })

  if (is.null(edat)) return(NULL)

  # Tag each result row with q so results can be combined across the sweep
  process_eval_results(edat, id) %>% dplyr::mutate(q = q)
}, mc.set.seed = TRUE, mc.cores = max(1, num_cores - 1)))

# --- Save results ------------------------------------------------------------
save(scenarios, file = "sim1-revision-temp.RData")
cat("temp file saved\n", file = out_filename, append = TRUE)

filename <- paste0("results/sim1-revision-q-", q, "-numsim-", sim_reps, "-",
                   format(Sys.time(), "%b-%d-%Y"), ".RData")
save(scenarios, file = filename)
cat("saved file\n", file = out_filename, append = TRUE)
