# =============================================================================
# Soldiering Application: Pilot Data Analysis (No Cross-Fitting)
# =============================================================================
#
# Uses ALL controls as the pilot (training) sample — no sample splitting.
# This is the "double-dip" baseline where the same controls appear in both
# the pilot and analysis sets.
#
# Removes covariates with missingness in the control group before estimation.
#
# Usage: Rscript pilot-data-analysis.R
# =============================================================================

source("../functions/sim-utils.R")
source("../functions/BART-features.R")
source("../functions/randomForestFeatures.R")
source("../functions/sim-estimation-funcs.R")
source("../functions/sim-eval-funcs.R")
source("soldiering-setup.R")

load_sim_packages()

# --- Load data ----------------------------------------------------------------

dat <- load_soldiering_data("blattman_claude.csv")
covs_list <- define_soldiering_covariates(dat$raw.data)
covs <- covs_list$full_covs

df <- dat$df

# Use controls to identify covariates with missingness, then exclude them
df.pilot <- df %>% dplyr::filter(Z == 0)
na_cols <- apply(df.pilot, 2, function(x) sum(is.na(x)))
na_cols <- na_cols[na_cols > 0]
covs <- covs[!covs %in% names(na_cols)]

# Remove rows with missing Z or Y, select analysis columns
df <- df %>%
  dplyr::filter(!is.na(Z) & !is.na(Y)) %>%
  dplyr::select(dplyr::all_of(c("Y", "Z", covs)))

df.pilot <- df.pilot %>%
  dplyr::filter(!is.na(Z) & !is.na(Y)) %>%
  dplyr::select(dplyr::all_of(c("Y", "Z", covs)))

# --- Estimation (no cross-fitting) --------------------------------------------

set.seed(12)

edat <- eval_data(dat = df, pilot.dat = df.pilot,
                  treat.true = 0, verbose = TRUE,
                  covs = covs, dataset = "soldiering-dbldip")

# --- Save results -------------------------------------------------------------

save(edat, file = "results/multi-dataset-results-educ-dbldip.RData")
