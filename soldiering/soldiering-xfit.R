# =============================================================================
# Soldiering Application: Cross-Fit Analysis (Algorithm 1)
# =============================================================================
#
# Blattman & Annan (2010) child soldiering in Uganda.
# Outcome: years of education (Section 5.2 of Shen et al. 2025).
#
# Performs 10-repeat 2-fold cross-fitting:
#   - Split controls in half randomly
#   - Each half takes a turn as pilot (training) sample
#   - Average estimates across folds and repeats
#
# Usage: Rscript soldiering-xfit.R
# =============================================================================

source("../functions/sim-utils.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")
source("../functions/cross-fit.R")
source("soldiering-setup.R")

load_sim_packages()

# --- Load data ----------------------------------------------------------------

dat <- load_soldiering_data("blattman_claude.csv")
covs_list <- define_soldiering_covariates(dat$raw.data)
covs <- covs_list$full_covs

df <- dat$df %>%
  dplyr::filter(!is.na(Z) & !is.na(Y)) %>%
  dplyr::select(dplyr::all_of(c("Y", "Z", covs)))

data.c <- df %>% dplyr::filter(Z == 0)
data.t <- df %>% dplyr::filter(Z == 1)

# --- Setup --------------------------------------------------------------------

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
RNGkind("L'Ecuyer-CMRG")
set.seed(12)

num_cores <- setup_parallel(reserve = 1)
out_filename <- create_log("soldiering-xfit")

n_repeat <- 10
seeder <- 2025 - 11 - 3
seeds <- seeder * seq_len(n_repeat)

cat("Starting soldiering XFIT\n", file = out_filename, append = TRUE)

# --- Cross-fit estimation (10 random splits) ----------------------------------

out <- parallel::mclapply(seq_len(n_repeat), function(i) {
  set.seed(seeds[i])

  cat(paste("Starting repeat", i, "at", Sys.time(), "\n"),
      file = out_filename, append = TRUE)

  sample_split <- sample(nrow(data.c), round(0.5 * nrow(data.c)))
  data.c1 <- data.c[sample_split, ]
  data.c2 <- data.c[-sample_split, ]

  result <- run_cross_fit(pilot1 = data.c1, pilot2 = data.c2,
                          treat.dat = data.t, covs = covs,
                          trans = "none", id = i, seed = seeds[i])

  cat(paste("Finished repeat", i, "\n"), file = out_filename, append = TRUE)
  result
}, mc.set.seed = FALSE, mc.preschedule = TRUE, mc.cores = num_cores - 1)

# --- Save results -------------------------------------------------------------

filename <- paste0("results/soldiering-xfit-educ-", seeder - 1, ".RData")
save(out, file = filename)
