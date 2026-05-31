# =============================================================================
# WSC Application: Cross-Fit Analysis (Algorithm 1)
# =============================================================================
#
# Within-study comparison using Byrd (2021) data.
# Outcome: math test scores (Section 5.1 of Shen et al. 2025).
# Experimental benchmark ATT: 0.79 (math), 2.18 (vocab).
#
# Each repeat runs TWO cross-fits:
#   1. Untransformed covariates
#   2. Log-transformed covariates (log(x+1) for continuous non-negative)
# Results are averaged across 10 repeats.
#
# Usage: Rscript wsc-xfit.R
# =============================================================================

source("../functions/sim-utils.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")
source("../functions/cross-fit.R")
source("wsc-setup.R")

load_sim_packages()

# --- Load data ----------------------------------------------------------------

outcome <- "math"
log_trans <- TRUE

wsc <- load_wsc_data(outcome)
data.pre <- wsc$data.pre
covs <- wsc$covs

# Log-transformed version of covariates
data.pre.log <- log_transform_covariates(data.pre)

# Split into control and treated
data.c <- data.pre %>% dplyr::filter(Z == 0)
data.t <- data.pre %>% dplyr::filter(Z == 1)

data.c.log <- data.pre.log %>% dplyr::filter(Z == 0)
data.t.log <- data.pre.log %>% dplyr::filter(Z == 1)

# --- Setup --------------------------------------------------------------------

set.seed(12)

num_cores <- setup_parallel(reserve = 1)
out_filename <- create_log("wsc-xfit")

n_repeat <- 10
seeder <- 2025 - 11 - 3
seeds <- seeder * seq_len(n_repeat)

cat(paste("Starting WSC XFIT with log_trans",
          ifelse(log_trans, "LOG", "EXP"), "\n"),
    file = out_filename, append = TRUE)

# --- Cross-fit estimation (10 random splits) ----------------------------------

out <- parallel::mclapply(seq_len(n_repeat), function(i) {
  cat(paste("Starting repeat", i, "at", Sys.time(), "\n"),
      file = out_filename, append = TRUE)

  set.seed(seeds[i])
  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))

  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]

  data.c1.log <- data.c.log[sample_split, ]
  data.c2.log <- data.c.log[-sample_split, ]

  # Untransformed cross-fit
  out.notrans <- run_cross_fit(pilot1 = data.c1, pilot2 = data.c2,
                               treat.dat = data.t, covs = covs,
                               trans = "none", id = i, dataset = "wsc")

  # Log-transformed cross-fit
  out.log <- run_cross_fit(pilot1 = data.c1.log, pilot2 = data.c2.log,
                           treat.dat = data.t.log, covs = covs,
                           trans = "log", id = i, dataset = "wsc")

  cat(paste("Finished repeat", i, "\n"), file = out_filename, append = TRUE)
  list(out.notrans = out.notrans, out.log = out.log)
}, mc.set.seed = TRUE, mc.cores = num_cores - 1)

# --- Save results -------------------------------------------------------------

if (outcome == "math") {
  if (log_trans) {
    save(out, file = "results/wsc-math-xfit.RData")
  } else {
    save(out, file = "results/wsc-math-xfit-exp.RData")
  }
} else if (outcome == "vocab") {
  if (log_trans) {
    save(out, file = "results/wsc-vocab-xfit.RData")
  } else {
    save(out, file = "results/wsc-vocab-xfit-exp.RData")
  }
}
