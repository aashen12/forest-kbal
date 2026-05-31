# =============================================================================
# WSC Application: Data Loading, Covariate Definitions, Transformations
# =============================================================================
#
# Within-study comparison using Byrd (2021) data.
# Two outcomes: math test scores (mathPost) and vocabulary (vocabPost).
# The experimental benchmark ATT is 0.79 for math and 2.18 for vocab.
#
# Shared by wsc-xfit.R.
# =============================================================================


#' Load and prepare the WSC observational study data.
#'
#' @param outcome "math" or "vocab" — determines which treatment and outcome
#' @return List with data.pre (prepared data frame) and covs (covariate names)
load_wsc_data <- function(outcome = "math") {
  load("data/wsc-obs-study-data.RData", envir = environment())

  if (outcome == "math") {
    data.obs$Y <- data.obs$mathPost
    data.pre <- data.obs
    data.pre$std_math <- scale(data.pre$mathPost)
    data.pre <- data.pre %>% dplyr::rename(Z = m.treat)
  } else if (outcome == "vocab") {
    data.obs$Y <- data.obs$vocabPost
    data.pre <- data.obs
    data.pre$std_math <- scale(data.pre$vocabPost)
    data.pre <- data.pre %>% dplyr::rename(Z = v.treat)
  }

  covs <- c(
    "female", "white", "black", "asian", "hisp", "married", "logAge", "income",
    "collegeS", "collegeM", "collegeD", "calc", "logBooks", "mathLike", "big5O", "big5C",
    "big5E", "big5A", "big5N", "AMAS", "logBDI", "MCS", "GSES", "vocabPre",
    "mathPre"
  )

  list(data.pre = data.pre, covs = covs)
}


#' Apply log transformation to continuous covariates.
#'
#' Transforms numeric covariates with 3+ distinct values and non-negative
#' range using log(x + 1). Binary/categorical and Z/Y are left unchanged.
#'
#' @param data.pre Data frame to transform
#' @return Transformed data frame
log_transform_covariates <- function(data.pre) {
  data.pre %>%
    dplyr::mutate(dplyr::across(
      .cols = where(is.numeric) & !dplyr::any_of(c("Z", "Y")),
      .fns = ~ if (dplyr::n_distinct(.x) >= 3 && min(.x, na.rm = TRUE) >= 0) {
        log(.x + 1)
      } else {
        .x
      }
    ))
}
