# =============================================================================
# Cross-Fitting Utilities for Applied Data Analysis
# =============================================================================
#
# Functions for 2-fold cross-fitting (Algorithm 1 in Shen et al. 2025)
# applied to real datasets (soldiering, WSC). Handles the extra output
# fields (explained variance) that applied analyses need beyond simulations.
#
# Contents:
#   process_applied_results()  - Flatten eval_data() output, keeping expl_var
#   run_cross_fit()            - 2-fold cross-fit with averaged estimates
# =============================================================================


# --- Result Processing --------------------------------------------------------

#' Flatten eval_data() output for applied data analysis.
#'
#' Like process_eval_results() in sim-utils.R, but also extracts the
#' explained variance vectors (rf_expl_var, bart_expl_var) into a separate
#' data frame for scree plots.
#'
#' @param edat List returned by eval_data()
#' @param id Repeat identifier (added as a column)
#' @return List with results_df and expl_var_df
process_applied_results <- function(edat, id = 999) {
  out <- lapply(seq_along(edat), function(i) {
    resi <- edat[[i]]
    expl_var_df <- data.frame(rf_expl_var = resi$rf_expl_var,
                              bart_expl_var = resi$bart_expl_var)
    resi_rest <- resi[!names(resi) %in% c("elbo_rf", "elbo_bart",
                                           "rf_expl_var", "bart_expl_var",
                                           "kbal_expl_var")]
    results_df <- dplyr::bind_rows(resi_rest) %>%
      dplyr::mutate(elbo_rf = resi$elbo_rf, elbo_bart = resi$elbo_bart)
    list(results_df = results_df, expl_var_df = expl_var_df)
  })
  results_df <- dplyr::bind_rows(lapply(out, function(x) x$results_df)) %>%
    dplyr::mutate(id = id)
  expl_var_df <- dplyr::bind_rows(lapply(out, function(x) x$expl_var_df)) %>%
    dplyr::mutate(id = id)
  list(results_df = results_df, expl_var_df = expl_var_df)
}


# --- Cross-Fitting ------------------------------------------------------------

#' Two-fold cross-fitting for applied data (Algorithm 1).
#'
#' Splits controls into two halves. Each half takes a turn as the pilot
#' (training) sample while the other enters the analysis set alongside
#' the treated units. Results are averaged across folds; standard errors
#' use the cross-fit formula: se = sqrt(se1^2/2 + se2^2/2).
#'
#' @param pilot1 First control half (data frame)
#' @param pilot2 Second control half (data frame)
#' @param treat.dat Treated units (data frame)
#' @param covs Character vector of covariate names
#' @param trans Transformation label (default: "none")
#' @param id Repeat identifier
#' @param seed Random seed passed to eval_data()
#' @param dataset Dataset identifier for eval_data() (default: "soldiering")
#' @return List with cfit.df (averaged results) and expl_var.df (averaged explained variance)
run_cross_fit <- function(pilot1, pilot2, treat.dat, covs,
                          trans = "none", id = 999, seed = 1,
                          dataset = "soldiering") {
  # Fold 1: pilot2 trains features, pilot1 joins analysis set
  analysis1 <- dplyr::bind_rows(treat.dat, pilot1)
  fit.2.1.raw <- eval_data(dat = analysis1, pilot.dat = pilot2,
                            treat.true = 0, verbose = TRUE,
                            covs = covs, dataset = dataset, seed = seed)
  print("Finished first fit")

  # Fold 2: pilot1 trains features, pilot2 joins analysis set
  analysis2 <- dplyr::bind_rows(treat.dat, pilot2)
  fit.1.2.raw <- eval_data(dat = analysis2, pilot.dat = pilot1,
                            treat.true = 0, verbose = TRUE,
                            covs = covs, dataset = dataset, seed = seed)
  print("Finished second fit")

  # Process both folds (once each, not redundantly)
  fit.2.1 <- process_applied_results(fit.2.1.raw, id = id)
  fit.1.2 <- process_applied_results(fit.1.2.raw, id = id)

  # Average explained variance across folds
  expl_var.df <- fit.1.2$expl_var_df
  expl_var.df$rf_expl_var <- (fit.1.2$expl_var_df$rf_expl_var +
                               fit.2.1$expl_var_df$rf_expl_var) / 2
  expl_var.df$bart_expl_var <- (fit.1.2$expl_var_df$bart_expl_var +
                                 fit.2.1$expl_var_df$bart_expl_var) / 2

  # Average numeric results across folds
  cfit.df <- fit.1.2$results_df
  num_cols <- sapply(cfit.df, is.numeric)
  cfit.df[, num_cols] <- (fit.1.2$results_df[, num_cols] +
                           fit.2.1$results_df[, num_cols]) / 2

  # Cross-fit SE formula
  cfit.df$se <- sqrt(fit.1.2$results_df$se^2 / 2 +
                      fit.2.1$results_df$se^2 / 2)
  cfit.df$lcl <- cfit.df$est.att - 1.96 * cfit.df$se
  cfit.df$hcl <- cfit.df$est.att + 1.96 * cfit.df$se
  cfit.df$trans <- trans

  list(cfit.df = cfit.df, expl_var.df = expl_var.df)
}
