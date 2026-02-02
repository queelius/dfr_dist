# Shared test fixtures for dfr.dist tests
#
# These helper functions are automatically loaded by testthat before running
# any tests. They provide common DFR distribution fixtures used across
# multiple test files.

# =============================================================================
# Distribution Fixtures
# =============================================================================

#' Create exponential DFR distribution for testing
#'
#' Exponential has constant hazard h(t) = lambda, H(t) = lambda * t
#'
#' @param lambda Optional rate parameter. If NULL, no default parameter.
#' @return A dfr_dist object
make_exponential_dfr <- function(lambda = NULL) {
  dfr_dist(
    rate = function(t, par, ...) {
      rep(par[1], length(t))  # constant hazard = lambda
    },
    par = lambda
  )
}

#' Create Weibull DFR distribution for testing
#'
#' Weibull has hazard h(t) = (k/sigma) * (t/sigma)^(k-1)
#' Cumulative hazard H(t) = (t/sigma)^k
#'
#' @param shape Optional shape parameter k
#' @param scale Optional scale parameter sigma
#' @return A dfr_dist object
make_weibull_dfr <- function(shape = NULL, scale = NULL) {
  par <- if (!is.null(shape) && !is.null(scale)) c(shape, scale) else NULL
  dfr_dist(
    rate = function(t, par, ...) {
      k <- par[1]      # shape
      sigma <- par[2]  # scale
      (k / sigma) * (t / sigma)^(k - 1)
    },
    par = par
  )
}

# =============================================================================
# Data Frame Fixtures
# =============================================================================

#' Create test data frame for exact observations
#'
#' @param times Vector of observation times
#' @return Data frame with t and delta columns
make_exact_data <- function(times) {
  data.frame(
    t = times,
    delta = rep(1, length(times))  # 1 = exact observation
  )
}

#' Create test data frame for right-censored observations
#'
#' @param times Vector of censoring times
#' @return Data frame with t and delta columns
make_censored_data <- function(times) {
  data.frame(
    t = times,
    delta = rep(0, length(times))  # 0 = right-censored
  )
}

#' Create mixed data (exact + censored)
#'
#' @param exact_times Vector of exact observation times
#' @param censored_times Vector of censoring times
#' @return Data frame with t and delta columns
make_mixed_data <- function(exact_times, censored_times) {
  rbind(
    make_exact_data(exact_times),
    make_censored_data(censored_times)
  )
}

# =============================================================================
# Analytical Functions for Verification
# =============================================================================

#' Analytical exponential log-likelihood for exact observations
#' loglik = sum(log(lambda) - lambda * t_i) = n * log(lambda) - lambda * sum(t)
exp_loglik_exact <- function(times, lambda) {
  n <- length(times)
  n * log(lambda) - lambda * sum(times)
}

#' Analytical exponential log-likelihood for right-censored observations
#' loglik = sum(-lambda * t_i) = -lambda * sum(t)
exp_loglik_censored <- function(times, lambda) {
  -lambda * sum(times)
}

#' Analytical exponential score for exact observations
#' d/dlambda [n * log(lambda) - lambda * sum(t)] = n/lambda - sum(t)
exp_score_exact <- function(times, lambda) {
  n <- length(times)
  n / lambda - sum(times)
}

#' Analytical exponential Hessian for exact observations
#' d^2/dlambda^2 = -n/lambda^2
exp_hessian_exact <- function(times, lambda) {
  n <- length(times)
  matrix(-n / lambda^2, nrow = 1, ncol = 1)
}
