# =============================================================================
# Combine SLURM Job-Array Outputs
# =============================================================================
#
# Reads every results-spawn/sim1-revision-q-<q>-n-<n>-task-*.RData file produced
# by the job array and row-binds them, one combined data frame per (q, n). Feed
# the combined file(s) to ../simulations/visualize-sim1-revision.R via the
# `results_files` named vector.
#
# Naming conventions:
#   * Input task files:
#       - new   :  sim1-revision-q-<q>-n-<n>-task-<id>.RData
#       - legacy:  sim1-revision-q-<q>-task-<id>.RData         (no -n-; treated as n = 1000)
#     Both are picked up automatically.
#   * Output combined files:
#       - n = 1000 :  sim1-revision-q-<q>-combined-<num_reps>.RData      (LEGACY -- matches
#                       pre-existing combined files exactly so they are preserved /
#                       overwritten in place rather than producing confusing duplicates)
#       - any other n (e.g. 5000):
#                       sim1-revision-q-<q>-n-<n>-combined-<num_reps>.RData
#   So the n = 1000 and n = 5000 sweeps live alongside each other under
#   different names with no risk of collision.
#
# Usage:
#   Rscript combine-spawn-results.R
# =============================================================================

library(tidyverse)

spawn_dir <- "results-spawn"

# Match both the new pattern (-n-<n>- present) and the legacy pattern (no -n-).
task_files <- list.files(
  spawn_dir,
  pattern = "^sim1-revision-q-[0-9]+(-n-[0-9]+)?-task-[0-9]+\\.RData$",
  full.names = TRUE
)
if (length(task_files) == 0) {
  stop("No task files found in '", spawn_dir,
       "'. Run the array first: sbatch sim1-revision-array.sh <q> <reps_per_job> [n]")
}

# Pull (q, n) out of each filename. n is missing in legacy files -> default 1000.
info <- tibble(
  file = task_files,
  q    = as.integer(str_match(basename(task_files), "q-([0-9]+)")[, 2]),
  n    = as.integer(str_match(basename(task_files), "-n-([0-9]+)-")[, 2])
)
info$n[is.na(info$n)] <- 1000L

message(sprintf("Found %d task files across %d (q, n) combination(s).",
                nrow(info), nrow(dplyr::distinct(info, q, n))))

# One combined file per (q, n) pair
keys <- info %>% dplyr::distinct(q, n) %>% dplyr::arrange(q, n)

for (i in seq_len(nrow(keys))) {
  qq <- keys$q[i]; nn <- keys$n[i]
  files_qn <- info$file[info$q == qq & info$n == nn]

  scenarios <- map_dfr(files_qn, function(f) {
    e <- new.env()
    load(f, envir = e)               # each file provides `scenarios`
    e$scenarios
  })

  n_reps  <- dplyr::n_distinct(scenarios$id)
  n_tasks <- dplyr::n_distinct(scenarios$task_id)

  # n = 1000 uses the legacy output name (no -n-) so it matches any pre-existing
  # combined files. Any other n uses the new -n-<n>- naming.
  out_name <- if (nn == 1000L) {
    sprintf("sim1-revision-q-%d-combined-%d.RData", qq, n_reps)
  } else {
    sprintf("sim1-revision-q-%d-n-%d-combined-%d.RData", qq, nn, n_reps)
  }
  out <- file.path(spawn_dir, out_name)
  save(scenarios, file = out)

  message(sprintf("  q = %3d, n = %5d:  %3d task files (%d tasks)  ->  %d unique reps  ->  %s",
                  qq, nn, length(files_qn), n_tasks, n_reps, basename(out)))
}

message("Done. Feed the combined file(s) to ../simulations/visualize-sim1-revision.R.")
