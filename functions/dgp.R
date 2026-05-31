# =============================================================================
# Data-Generating Processes for Simulation Studies
# =============================================================================
#
# Defines the DGPs used in the forest-kbal paper simulations.
# See Section 4 of Shen et al. (2025) for details.
#
# DGP 1 (make_data_sim1): Tarr & Imai (2025) / Wong & Chan (2018)
#   - 10 covariates with nonlinear transformations of latent Gaussians
#   - Logistic treatment assignment, heterogeneous treatment effects
#   - Used in Simulation 1 (main paper) and varsweep analysis
#
# DGP 1 high-dim (make_data_sim1_highdim): Simulation 1 revision
#   - Generalizes DGP 1 to arbitrary covariate dimension q
#   - Cycles the nonlinearities and the treatment/outcome model over
#     consecutive blocks of 10 latent Gaussians
#   - Signal rescaled by sqrt(#blocks) so overlap stays fixed as q grows
#   - Used in Simulation 1 revision (reviewer comment 3: effect of dim(X))
#
# DGP 2 (make_data_sim2): Kim et al. (2024)
#   - 6 covariates (3 correlated normals, 1 uniform, 1 chi-squared, 1 Bernoulli)
#   - Overlap controlled by sig.ep noise parameter
#   - True ATT = 0 by construction
#   - Used in Simulation 2 (appendix)
# =============================================================================


#' DGP 1: Tarr & Imai (2025) simulation
#'
#' Generates data with 10 nonlinearly transformed covariates,
#' logistic treatment assignment, and heterogeneous treatment effects.
#' True ATT varies across simulation draws.
#'
#' @param n Sample size
#' @return List with out.df (data frame with X1-X10, Z, Y) and true.att
make_data_sim1 <- function(n) {
  # Latent variables: W1, ..., W10 iid standard Gaussian
  ZZ <- MASS::mvrnorm(n, mu = rep(0, 10), Sigma = diag(10))
  colnames(ZZ) <- paste0("ZZ", 1:10)

  # Nonlinear covariate transformations (Section 4.1)
  X <- data.frame(
    X1  = exp(ZZ[, 1] / 2),
    X2  = ZZ[, 2] / (1 + exp(ZZ[, 1])),
    X3  = (ZZ[, 1] * ZZ[, 3] / 25 + 0.6)^3,
    X4  = (ZZ[, 2] + ZZ[, 4] + 20)^2,
    X5  = ZZ[, 5],
    X6  = ZZ[, 6],
    X7  = ZZ[, 7],
    X8  = ZZ[, 8],
    X9  = ZZ[, 9],
    X10 = ZZ[, 10]
  )

  # Treatment assignment: P(Z=1 | W) = logit^{-1}(-W1 - 0.1*W4)
  logits <- -ZZ[, 1] - 0.1 * ZZ[, 4]
  prob_Z <- exp(logits) / (1 + exp(logits))
  Z <- rbinom(n, 1, prob_Z)

  # Potential outcomes: Y(Z) = 200 + 10Z + (1.5Z - 0.5)(27.4W1 + 13.7W2 + 13.7W3 + 13.7W4) + eps
  interaction_term <- 27.4 * ZZ[,1] + 13.7 * ZZ[,2] + 13.7 * ZZ[,3] + 13.7 * ZZ[,4]
  Y0 <- 200 - 0.5 * interaction_term + rnorm(n)
  Y1 <- Y0 + 10 + 1.5 * interaction_term
  Y  <- Z * Y1 + (1 - Z) * Y0

  out.df <- data.frame(X, Z = Z, Y = Y)
  true.att <- mean(Y1[Z == 1]) - mean(Y0[Z == 1])

  list(out.df = out.df, true.att = true.att)
}


#' DGP 1 (high-dimensional variant): Tarr & Imai (2025) with q covariates
#'
#' Generalizes make_data_sim1() to arbitrary covariate dimension q. The data
#' are built from q latent standard Gaussians grouped into consecutive blocks
#' of 10. Within each block, the same nonlinear transformations as DGP 1 are
#' applied (positions 1-4 nonlinear and confounding; positions 5-10 identity
#' noise), and each block contributes the same terms to the propensity and
#' outcome models. The summed propensity and outcome linear predictors are
#' rescaled by sqrt(number of blocks) so that propensity overlap and outcome
#' signal-to-noise stay approximately constant as q grows -- this isolates the
#' effect of dimensionality rather than confounding strength.
#'
#' A final partial block (when q is not a multiple of 10) is allowed: only the
#' positions that exist are constructed. With q = 10 this reduces exactly to
#' make_data_sim1().
#'
#' Used in the Simulation 1 revision (reviewer comment 3: how does covariate
#' dimensionality q affect performance?).
#'
#' @param n Sample size
#' @param q Number of covariates (built in blocks of 10; partial final block allowed)
#' @return List with out.df (data frame with X1..Xq, Z, Y) and true.att
make_data_sim1_highdim <- function(n, q = 10) {
  n_blocks <- ceiling(q / 10)

  # Latent variables: W1, ..., Wq iid standard Gaussian
  W <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = diag(q))
  W <- matrix(W, nrow = n, ncol = q)  # guard against vector result when q == 1

  X <- matrix(NA_real_, nrow = n, ncol = q)
  ps_logit <- numeric(n)
  interaction_term <- numeric(n)

  for (g in seq_len(n_blocks)) {
    offset <- (g - 1) * 10
    bsize  <- min(10, q - offset)      # block size (10, or smaller for partial block)
    idx    <- offset + seq_len(bsize)  # global covariate indices for this block
    w      <- W[, idx, drop = FALSE]   # latent Gaussians for this block

    # Nonlinear covariate transformations (cycled per block).
    # Positions 5..bsize stay identity, so initialize with the raw latents.
    Xb <- w
    if (bsize >= 1) Xb[, 1] <- exp(w[, 1] / 2)
    if (bsize >= 2) Xb[, 2] <- w[, 2] / (1 + exp(w[, 1]))
    if (bsize >= 3) Xb[, 3] <- (w[, 1] * w[, 3] / 25 + 0.6)^3
    if (bsize >= 4) Xb[, 4] <- (w[, 2] + w[, 4] + 20)^2
    X[, idx] <- Xb

    # Treatment uses positions 1 and 4; outcome uses positions 1-4 (cycled per block).
    if (bsize >= 1) ps_logit <- ps_logit - w[, 1]
    if (bsize >= 4) ps_logit <- ps_logit - 0.1 * w[, 4]

    if (bsize >= 1) interaction_term <- interaction_term + 27.4 * w[, 1]
    if (bsize >= 2) interaction_term <- interaction_term + 13.7 * w[, 2]
    if (bsize >= 3) interaction_term <- interaction_term + 13.7 * w[, 3]
    if (bsize >= 4) interaction_term <- interaction_term + 13.7 * w[, 4]
  }

  # Rescale accumulated signal so overlap / signal-to-noise stay fixed across q
  ps_logit         <- ps_logit / sqrt(n_blocks)
  interaction_term <- interaction_term / sqrt(n_blocks)

  X <- as.data.frame(X)
  names(X) <- paste0("X", seq_len(q))

  # Treatment assignment: P(Z=1 | W) = logit^{-1}(rescaled block sum)
  prob_Z <- exp(ps_logit) / (1 + exp(ps_logit))
  Z <- rbinom(n, 1, prob_Z)

  # Potential outcomes: same form as DGP 1, with the rescaled interaction term
  Y0 <- 200 - 0.5 * interaction_term + rnorm(n)
  Y1 <- Y0 + 10 + 1.5 * interaction_term
  Y  <- Z * Y1 + (1 - Z) * Y0

  out.df <- data.frame(X, Z = Z, Y = Y)
  true.att <- mean(Y1[Z == 1]) - mean(Y0[Z == 1])

  list(out.df = out.df, true.att = true.att)
}


#' DGP 2: Kim et al. (2024) simulation with varying overlap
#'
#' Generates data with 6 covariates. Overlap between treated and control
#' groups is controlled by sig.ep (larger = more overlap).
#' True ATT = 0 by construction.
#'
#' @param n Sample size
#' @param sig.ep Noise in treatment assignment (30 = moderate overlap, 100 = strong overlap)
#' @return List with out.df (data frame with X1-X6, Z, Y) and true.att (always 0)
make_data_sim2 <- function(n, sig.ep) {
  # Correlated normal covariates (X1, X2, X3)
  sig.123 <- diag(c(2, 1, 1))
  sig.123[1, 2] <- 1; sig.123[1, 3] <- -1; sig.123[2, 3] <- -0.5
  sig.123 <- Matrix::forceSymmetric(sig.123)

  X.123 <- as.matrix(MASS::mvrnorm(n, mu = rep(0, 3), Sigma = sig.123))
  colnames(X.123) <- paste0("X", 1:3)
  X4 <- runif(n, -3, 3)
  X5 <- rchisq(n, 1)
  X6 <- rbinom(n, 1, 0.5)
  X <- cbind(X.123, X4, X5, X6)

  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]

  # Treatment: nonlinear function of covariates + noise controlled by sig.ep
  expression_value <- X1^2 + 2*X2^2 - 2*X3^2 - (X4 + 1)^3 - 0.5*log(X5 + 10) + X6 - 1.5
  Z <- ifelse(expression_value + rnorm(n, 0, sig.ep) > 0, 1, 0)

  # Outcome (no treatment effect by construction)
  Y <- (X.123[, 1] + X.123[, 2] + X5)^2 + rnorm(n, 0, 1)

  out.df <- data.frame(X, Z = Z, Y = Y)

  list(out.df = out.df, true.att = 0)
}
