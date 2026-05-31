# =============================================================================
# Simulation 1 Revision -- SLURM Job-Array Worker
# =============================================================================
#
# Array-friendly version of ../simulations/sim1-revision.R. Instead of running
# all Monte Carlo reps inside one job, this script runs a SMALL batch of reps
# (reps_per_job) as a single SLURM array task. Submitting it with
#   sbatch --array=1-100 sim1-revision-array.sh <q> <reps_per_job> [n]
# spawns 100 independent tasks; combine-spawn-results.R stitches their outputs
# back into one data frame per (q, n).
#
# Each task:
#   - reads its index from SLURM_ARRAY_TASK_ID,
#   - runs reps_per_job reps with globally-unique simulation ids,
#   - sets an EXPLICIT seed per rep (SEED_BASE + global_id), so every rep across
#     the entire array gets a distinct, reproducible seed -- no chance of two
#     tasks generating identical data,
#   - writes results-spawn/sim1-revision-q-<q>-n-<n>-task-<id>.RData.
#
# Arguments:
#   q             Number of covariates (e.g. 10, 50, 100)
#   reps_per_job  Monte Carlo reps run by this single task
#   n             Sample size for each rep (analysis sample; the pilot sample
#                 is drawn the same way and filtered to controls, ~n/2).
#                 Optional, defaults to 1000 -- so old invocations still work.
#
# The estimation pipeline is identical to sim1-revision.R; only the
# orchestration (seeding, batching, sample size, output path) differs. The
# original sim1-revision.R script and its results are not touched.
#
# Usage (cluster, via the array wrapper):
#   sbatch sim1-revision-array.sh <q> <reps_per_job> [n]
# Local single-task test:
#   SLURM_ARRAY_TASK_ID=1 Rscript sim1-revision-array.R 50 2          # n=1000 (default)
#   SLURM_ARRAY_TASK_ID=1 Rscript sim1-revision-array.R 50 2 5000     # n=5000
# =============================================================================

source("../../functions/sim-utils.R")
load_sim_packages()

source("../../functions/dgp.R")
source("../../functions/BART-features.R")
source("../../functions/randomForestFeatures.R")
source("../../functions/sim-estimation-funcs.R")
source("../../functions/sim-eval-funcs.R")

# --- Parse inputs ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
q            <- as.integer(args[1])
reps_per_job <- as.integer(args[2])
# Sample size (optional 3rd arg; default 1000 preserves the original behavior)
n <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 1000L
stopifnot(!is.na(q), !is.na(reps_per_job), !is.na(n),
          q > 0, reps_per_job > 0, n > 0)

# Array task id (default 1 makes local testing painless)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))

# Cores allocated to this single task; used for within-task mclapply.
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# Per-rep seed base. Every rep across the entire array gets a unique seed of
# the form SEED_BASE + global_id, where global_id is unique by construction.
SEED_BASE <- 23967

out_dir <- "results-spawn"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat(sprintf("[task %d] q=%d, n=%d, reps_per_job=%d, num_cores=%d\n",
            task_id, q, n, reps_per_job, num_cores))

# --- Run this task's batch of reps -------------------------------------------
scenarios <- dplyr::bind_rows(parallel::mclapply(seq_len(reps_per_job), function(local_id) {
  # Globally unique id across the whole array (task 1 reps 1..N, task 2 reps N+1..2N, ...)
  global_id <- (task_id - 1L) * reps_per_job + local_id

  # Distinct, reproducible seed for THIS rep, independent of mclapply internals
  set.seed(SEED_BASE + global_id)

  bdat.obj  <- make_data_sim1_highdim(n = n, q = q)
  bdat      <- bdat.obj$out.df
  # Pilot: same generative process, then keep controls (~n/2 by construction)
  pilot.dat <- make_data_sim1_highdim(n = n, q = q)$out.df %>% dplyr::filter(Z == 0)

  edat <- tryCatch(
    eval_data(dat = bdat, pilot.dat = pilot.dat,
              treat.true = bdat.obj$true.att, verbose = FALSE,
              dataset = "simulation"),
    error = function(e) {
      cat(sprintf("[task %d, global_id %d] eval_data failed: %s\n",
                  task_id, global_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(edat)) return(NULL)

  process_eval_results(edat, global_id) %>%
    dplyr::mutate(q = q, n = n, task_id = task_id)
}, mc.cores = max(1L, num_cores), mc.set.seed = FALSE))   # we seed explicitly above

# --- Save this task's output -------------------------------------------------
filename <- file.path(out_dir,
                      sprintf("sim1-revision-q-%d-n-%d-task-%04d.RData", q, n, task_id))
save(scenarios, file = filename)
cat(sprintf("[task %d] wrote %s (%d rows)\n", task_id, filename, nrow(scenarios)))
