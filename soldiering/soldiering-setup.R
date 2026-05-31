# =============================================================================
# Soldiering Application: Data Loading and Covariate Definitions
# =============================================================================
#
# Blattman & Annan (2010) child soldiering in Uganda.
# Shared by soldiering-xfit.R and pilot-data-analysis.R.
#
# Covariate groups:
#   AC_vars  — age and geographic region indicators
#   ACG_vars — AC_vars plus G1-G8 war experience indicators
#   ctrls    — household and parental controls (polynomial expansions)
#   ctrl_hh  — household-level controls (raw)
#   ctrl_i   — individual-level parental/orphan controls
#   full_covs — union of all groups
# =============================================================================


#' Load and prepare the Blattman & Annan soldiering dataset.
#'
#' Reads blattman_claude.csv, renames variables to standard names
#' (Z = treatment, Y = outcome).
#'
#' @param path Path to blattman_claude.csv
#' @return List with df (prepared data frame) and raw.data (original data)
load_soldiering_data <- function(path = "blattman_claude.csv") {
  raw.data <- readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::rename(log.wage = lwage_mo)

  df <- raw.data %>%
    dplyr::mutate(wage_mo = exp(log.wage) - 1) %>%
    dplyr::rename(Z = abd, Y = educ)

  list(df = df, raw.data = raw.data)
}


#' Define covariate groups from the raw Blattman data.
#'
#' @param raw.data The raw data frame (needed for column name matching)
#' @return Named list of covariate vectors: AC_vars, ACG_vars, ctrls,
#'   ctrl_hh, ctrl_i, full_covs
define_soldiering_covariates <- function(raw.data) {
  AC_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$",
                  names(raw.data), value = TRUE)
  ACG_vars <- grep("^A1[4-9]$|^A2[0-9]$|^C_(ach|lan|kit|pad|amo|oru|paj)$|^G[1-8]_",
                   names(raw.data), value = TRUE)

  ctrls <- c("fthr_ed0", "fthr_ed4", "fthr_ed7", "mthr_ed0", "mthr_ed4", "mthr_ed7",
             "no_fthr96", "no_mthr96", "orphan96", "hh_fthr_frm", "hh_size96",
             "hh_size96_2", "hh_size96_3", "hh_size96_4", "hh_land", "hh_land_2",
             "hh_land_3", "hh_land_4", "hh_cattle", "hh_cattle_2", "hh_cattle_3",
             "hh_cattle_4", "hh_stock", "hh_stock_2", "hh_stock_3", "hh_stock_4", "hh_plow")

  ctrl_hh <- c("age", "hh_fthr_frm", "hh_size96", "hh_land", "landrich",
               "hh_cattle", "hh_stock", "hh_plow")

  ctrl_i <- c("fthr_ed0", "fthr_ed", "mthr_ed0", "mthr_ed", "no_fthr96",
              "no_mthr96", "orphan96")

  full_covs <- c(AC_vars, ACG_vars, ctrls, ctrl_hh, ctrl_i)

  list(AC_vars = AC_vars, ACG_vars = ACG_vars, ctrls = ctrls,
       ctrl_hh = ctrl_hh, ctrl_i = ctrl_i, full_covs = full_covs)
}
