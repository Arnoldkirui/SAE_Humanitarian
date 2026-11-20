# ==============================================================================
# 00_functions.R
#
# Small area estimation models, implemented directly rather than called from a
# package.
#
# WHY IMPLEMENTED HERE
#
# The `sae` package does all of this and is the right choice for production
# work. These implementations exist so the mechanics are visible: the shrinkage
# weight, the variance components, the MSE decomposition. An analyst who has
# seen gamma computed can explain to a programme manager why one LGA's estimate
# moved a long way and another's barely moved at all, which is the conversation
# that actually decides whether SAE gets used.
#
# Everything here is checked against the standard references in
# docs/references.md.
# ==============================================================================

# ==============================================================================
# FAY-HERRIOT AREA-LEVEL MODEL
# ==============================================================================
#
# The model has two levels:
#
#   sampling:   y_i   = theta_i + e_i      e_i ~ N(0, psi_i),  psi_i KNOWN
#   linking:    theta_i = x_i'beta + v_i   v_i ~ N(0, sigma_v2)
#
# y_i is the direct survey estimate for area i, psi_i its sampling variance,
# and theta_i the quantity you actually want. The model says the truth is
# predictable from auxiliary data up to an area-specific departure v_i.
#
# What makes it work is that psi_i is treated as known. That is both the
# model's strength and its weakest assumption: psi_i comes from the survey, and
# in a small area it is estimated from very few observations, so it is itself
# uncertain. The model does not propagate that uncertainty.

#' REML log-likelihood for the Fay-Herriot variance component
#'
#' Restricted rather than full maximum likelihood because ML is biased downward
#' for variance components, and with 40-60 areas that bias is not negligible.
#' A downward-biased sigma_v2 pulls every estimate too hard toward the
#' regression line, which is exactly the failure an SAE sceptic will look for.
fh_reml_loglik <- function(sigma_v2, y, X, psi) {
  V_diag <- sigma_v2 + psi
  Vinv <- 1 / V_diag

  XtVinvX <- crossprod(X, X * Vinv)
  XtVinvy <- crossprod(X, y * Vinv)

  beta <- tryCatch(solve(XtVinvX, XtVinvy), error = function(e) NULL)
  if (is.null(beta)) return(-Inf)

  resid <- as.vector(y - X %*% beta)

  -0.5 * (sum(log(V_diag)) +
          determinant(XtVinvX, logarithm = TRUE)$modulus[1] +
          sum(resid^2 * Vinv))
}

#' Fit a Fay-Herriot model
#'
#' @param y   direct estimates, one per area
#' @param X   design matrix of area-level auxiliary variables (include intercept)
#' @param psi sampling variances of the direct estimates, assumed known
#' @return list with EBLUPs, shrinkage weights, MSE and the fitted components
fit_fay_herriot <- function(y, X, psi, max_sigma = NULL) {

  stopifnot(length(y) == length(psi), nrow(X) == length(y))
  if (any(psi <= 0, na.rm = TRUE)) {
    stop("Sampling variances must be strictly positive. An area with psi = 0 ",
         "is claiming a direct estimate with no error, which no finite sample ",
         "supports.")
  }

  m <- length(y)
  p <- ncol(X)

  if (is.null(max_sigma)) max_sigma <- 10 * var(y)

  opt <- optimize(fh_reml_loglik, interval = c(0, max_sigma), maximum = TRUE,
                  y = y, X = X, psi = psi)
  sigma_v2 <- opt$maximum

  # A variance component estimated at (or very near) zero is a real result, not
  # a numerical failure: it says the auxiliary variables explain the between-
  # area variation entirely, so the model has no area-specific departure to
  # add. Every EBLUP then collapses onto the regression line. Flag it, because
  # it usually means the auxiliaries are too good to be true - often one of
  # them is derived from the outcome.
  sigma_at_zero <- sigma_v2 < 1e-8

  V_diag <- sigma_v2 + psi
  Vinv   <- 1 / V_diag

  XtVinvX <- crossprod(X, X * Vinv)
  XtVinvX_inv <- solve(XtVinvX)
  beta <- as.vector(XtVinvX_inv %*% crossprod(X, y * Vinv))

  x_beta <- as.vector(X %*% beta)

  # The shrinkage weight. gamma near 1 means the direct estimate is trusted;
  # near 0 means the model takes over. This single number is the most useful
  # diagnostic the method produces.
  gamma <- sigma_v2 / (sigma_v2 + psi)

  eblup <- gamma * y + (1 - gamma) * x_beta

  # --- MSE, Prasad-Rao / Datta-Lahiri decomposition --------------------------
  #
  # g1  the irreducible part: how much better than the direct estimate you can
  #     possibly do given the model
  # g2  the cost of having estimated beta
  # g3  the cost of having estimated sigma_v2
  #
  # Reporting the MSE matters more here than in most modelling. The whole claim
  # of SAE is "this estimate is more reliable than the direct one" - and that
  # claim is only checkable if the MSE is on the table next to psi.

  g1 <- gamma * psi
  g2 <- sapply(seq_len(m), function(i) {
    xi <- X[i, , drop = FALSE]
    as.numeric((1 - gamma[i])^2 * (xi %*% XtVinvX_inv %*% t(xi)))
  })

  # Asymptotic variance of the REML estimator of sigma_v2
  var_sigma_v2 <- 2 / sum(1 / V_diag^2)
  g3 <- psi^2 / (sigma_v2 + psi)^3 * var_sigma_v2

  mse <- g1 + g2 + 2 * g3

  list(
    eblup        = eblup,
    gamma        = gamma,
    beta         = beta,
    sigma_v2     = sigma_v2,
    sigma_at_zero = sigma_at_zero,
    mse          = mse,
    g1 = g1, g2 = g2, g3 = g3,
    x_beta       = x_beta,
    loglik       = opt$objective,
    n_areas      = m,
    n_params     = p
  )
}


# ==============================================================================
# BATTESE-HARTER-FULLER UNIT-LEVEL MODEL
# ==============================================================================
#
#   y_ij = x_ij'beta + v_i + e_ij
#
# with v_i ~ N(0, sigma_v2) the area effect and e_ij ~ N(0, sigma_e2) the
# household residual.
#
# The area mean predictor is
#
#   theta_i = Xbar_i'beta + gamma_i * (ybar_i - xbar_i'beta)
#
# where Xbar_i is the POPULATION mean of the covariates in area i and xbar_i
# the sample mean.
#
# THE REQUIREMENT PEOPLE UNDERESTIMATE
#
# Xbar_i has to come from outside the survey - a census, an administrative
# register, a satellite product aggregated to the area. And the covariates have
# to mean the same thing in both sources. A household roster question asked one
# way in your survey and another way in the census is not the same variable,
# and using it as though it were introduces a bias no diagnostic in this file
# will catch.
#
# That requirement, not the mathematics, is why area-level models are more
# common in humanitarian work than unit-level ones.

#' Fit a Battese-Harter-Fuller nested error model
#'
#' @param formula   fixed-effects formula, e.g. y ~ x1 + x2
#' @param data      unit-level sample data
#' @param area_var  name of the area identifier column
#' @param pop_means data frame of population covariate means, one row per area
#' @param pop_sizes data frame with area and population size N_i
fit_bhf <- function(formula, data, area_var, pop_means, pop_sizes) {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for the unit-level model.\n",
         "install.packages('lme4')")
  }

  re_formula <- stats::update(formula,
                              paste(". ~ . + (1 |", area_var, ")"))

  fit <- lme4::lmer(re_formula, data = data, REML = TRUE)

  vc <- as.data.frame(lme4::VarCorr(fit))
  sigma_v2 <- vc$vcov[vc$grp == area_var]
  sigma_e2 <- vc$vcov[vc$grp == "Residual"]

  beta <- lme4::fixef(fit)
  xvars <- names(beta)[names(beta) != "(Intercept)"]

  missing_pop <- setdiff(xvars, names(pop_means))
  if (length(missing_pop) > 0) {
    stop("Population means are missing for: ",
         paste(missing_pop, collapse = ", "),
         "\nThe unit-level model cannot run without a population mean for ",
         "every covariate. This is the binding constraint on the method, not ",
         "a data-tidying problem.")
  }

  areas <- pop_means[[area_var]]

  # Sample means and counts per area
  samp <- data |>
    dplyr::group_by(.data[[area_var]]) |>
    dplyr::summarise(
      n_i = dplyr::n(),
      ybar = mean(.data[[all.vars(formula)[1]]]),
      dplyr::across(dplyr::all_of(xvars), mean, .names = "xbar_{.col}"),
      .groups = "drop"
    )

  out <- pop_means |>
    dplyr::left_join(samp, by = area_var) |>
    dplyr::left_join(pop_sizes, by = area_var)

  Xbar_beta <- beta[["(Intercept)"]] +
    rowSums(sapply(xvars, function(v) out[[v]] * beta[[v]]))

  xbar_beta <- beta[["(Intercept)"]] +
    rowSums(sapply(xvars, function(v) out[[paste0("xbar_", v)]] * beta[[v]]))

  gamma <- sigma_v2 / (sigma_v2 + sigma_e2 / out$n_i)

  eblup <- Xbar_beta + gamma * (out$ybar - xbar_beta)

  # MSE, leading terms only (Rao & Molina 2015, ch. 7). The full expression
  # carries further terms for the estimation of beta and the variance
  # components; a bootstrap is the honest route for a real submission.
  mse_g1 <- gamma * sigma_e2 / out$n_i
  mse <- mse_g1

  data.frame(
    area      = areas,
    n_i       = out$n_i,
    N_i       = out$N_i,
    ybar      = out$ybar,
    synthetic = Xbar_beta,
    gamma     = gamma,
    eblup     = eblup,
    mse       = mse,
    stringsAsFactors = FALSE
  ) |>
    dplyr::rename_with(~ area_var, .cols = "area")
}


# ==============================================================================
# SHARED HELPERS
# ==============================================================================

#' Coefficient of variation of an estimate, as a percentage
#'
#' The standard reportability screen. Many statistical offices suppress or flag
#' estimates above 25-35% CV; the exact threshold is a policy choice and should
#' be stated rather than assumed.
cv_pct <- function(estimate, variance) {
  100 * sqrt(variance) / estimate
}

#' Design effect from a two-stage design
#'
#' deff = 1 + (b - 1) * rho, with b the average cluster size and rho the
#' intra-cluster correlation. Included because the direct variance in
#' R/02 has to account for clustering, and ignoring it is the single most
#' common way a humanitarian survey overstates its own precision.
design_effect <- function(avg_cluster_size, icc) {
  1 + (avg_cluster_size - 1) * icc
}
