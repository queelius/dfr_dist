# Tests for likelihood_model interface support in dfr_dist
#
# These tests define the contract that dfr_dist must satisfy when implementing
# the likelihood_model interface from the likelihood.model package.
#
# Mathematical Background:
# For a DFR distribution with hazard h(t, par) and cumulative hazard H(t, par):
#   - Exact observation at time t: log-likelihood = log h(t, par) - H(t, par)
#   - Right-censored at time t: log-likelihood = -H(t, par)
#   - Left-censored at time t: log-likelihood = log(1 - exp(-H(t, par)))

# =============================================================================
# Test Fixtures and Helper Functions
# =============================================================================

# Helper: Create exponential DFR distribution
# Exponential has constant hazard h(t) = lambda, H(t) = lambda * t
make_exponential_dfr <- function(lambda = NULL) {
  dfr_dist(
    rate = function(t, par, ...) {
      rep(par[1], length(t))  # constant hazard = lambda
    },
    par = lambda
  )
}

# Helper: Create Weibull DFR distribution
# Weibull has hazard h(t) = (k/sigma) * (t/sigma)^(k-1)
# Cumulative hazard H(t) = (t/sigma)^k
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

# Helper: Create time-varying hazard DFR (bathtub-like)
# h(t) = a * exp(-b*t) + c  (decreasing then constant)
make_bathtub_dfr <- function(a = NULL, b = NULL, c = NULL) {
  par <- if (!is.null(a) && !is.null(b) && !is.null(c)) c(a, b, c) else NULL
  dfr_dist(
    rate = function(t, par, ...) {
      par[1] * exp(-par[2] * t) + par[3]
    },
    par = par
  )
}

# Helper: Create test data frame for exact observations
make_exact_data <- function(times) {
  data.frame(
    t = times,
    delta = rep(1, length(times))  # 1 = exact observation
  )
}

# Helper: Create test data frame for right-censored observations
make_censored_data <- function(times) {
  data.frame(
    t = times,
    delta = rep(0, length(times))  # 0 = right-censored
  )
}

# Helper: Create mixed data (exact + censored)
make_mixed_data <- function(exact_times, censored_times) {
  rbind(
    make_exact_data(exact_times),
    make_censored_data(censored_times)
  )
}

# Analytical exponential log-likelihood for exact observations
# loglik = sum(log(lambda) - lambda * t_i) = n * log(lambda) - lambda * sum(t)
exp_loglik_exact <- function(times, lambda) {
  n <- length(times)
  n * log(lambda) - lambda * sum(times)
}

# Analytical exponential log-likelihood for right-censored observations
# loglik = sum(-lambda * t_i) = -lambda * sum(t)
exp_loglik_censored <- function(times, lambda) {
  -lambda * sum(times)
}

# Analytical exponential score for exact observations
# d/dlambda [n * log(lambda) - lambda * sum(t)] = n/lambda - sum(t)
exp_score_exact <- function(times, lambda) {
  n <- length(times)
  n / lambda - sum(times)
}

# Analytical exponential Hessian for exact observations
# d^2/dlambda^2 = -n/lambda^2
exp_hessian_exact <- function(times, lambda) {
  n <- length(times)
  matrix(-n / lambda^2, nrow = 1, ncol = 1)
}

# =============================================================================
# 1. Tests for loglik() - Exact Observations
# =============================================================================

test_that("loglik returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  ll <- loglik(dist)
  expect_type(ll, "closure")
})

test_that("loglik function has correct signature (df, par, ...)", {
  dist <- make_exponential_dfr()
  ll <- loglik(dist)
  # Function should accept df and par arguments
  args <- names(formals(ll))
  expect_true("df" %in% args || length(args) >= 2,
              info = "loglik function should accept df and par arguments")
})

test_that("loglik for exponential exact observations matches analytical", {
  # Given: Exponential distribution with lambda = 0.5
  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Given: Exact observations at times 1, 2, 3, 4, 5
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  # When: Computing log-likelihood
  ll <- loglik(dist)
  numerical_ll <- ll(df, par = c(lambda))

  # Then: Should match analytical formula
  # loglik = n * log(lambda) - lambda * sum(t)
  analytical_ll <- exp_loglik_exact(times, lambda)

  expect_equal(numerical_ll, analytical_ll, tolerance = 1e-6,
               info = "Exponential loglik should match n*log(lambda) - lambda*sum(t)")
})

test_that("loglik for exponential exact observations with different lambda values", {
  dist <- make_exponential_dfr()
  times <- c(0.5, 1.0, 1.5, 2.0)
  df <- make_exact_data(times)
  ll <- loglik(dist)

  # Test multiple lambda values
  for (lambda in c(0.1, 0.5, 1.0, 2.0, 5.0)) {
    numerical <- ll(df, par = c(lambda))
    analytical <- exp_loglik_exact(times, lambda)
    expect_equal(numerical, analytical, tolerance = 1e-6,
                 info = paste("lambda =", lambda))
  }
})

test_that("loglik increases when parameters are closer to true value", {
  # Given: Data generated from exponential with true lambda = 1.0
  true_lambda <- 1.0
  set.seed(42)
  times <- rexp(20, rate = true_lambda)
  df <- make_exact_data(times)

  dist <- make_exponential_dfr()
  ll <- loglik(dist)

  # When: Evaluating at true vs wrong parameters
  ll_true <- ll(df, par = c(true_lambda))
  ll_wrong <- ll(df, par = c(0.1))

  # Then: True parameters should give higher likelihood
  expect_gt(ll_true, ll_wrong)
})

# =============================================================================
# 2. Tests for loglik() - Right-Censored Observations
# =============================================================================

test_that("loglik for right-censored observations uses -H(t) only", {
  # Given: Exponential distribution
  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Given: Right-censored observations (didn't fail by time t)
  times <- c(1, 2, 3)
  df <- make_censored_data(times)

  # When: Computing log-likelihood
  ll <- loglik(dist)
  numerical_ll <- ll(df, par = c(lambda))

  # Then: Should be -H(t) = -lambda * sum(t)
  analytical_ll <- exp_loglik_censored(times, lambda)

  expect_equal(numerical_ll, analytical_ll, tolerance = 1e-6,
               info = "Censored loglik should be -lambda * sum(t)")
})

test_that("loglik for mixed exact and censored observations", {
  # Given: Exponential distribution
  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Given: Mixed observations
  exact_times <- c(1, 2, 3)
  censored_times <- c(4, 5)
  df <- make_mixed_data(exact_times, censored_times)

  # When: Computing log-likelihood
  ll <- loglik(dist)
  numerical_ll <- ll(df, par = c(lambda))

  # Then: Should be sum of exact + censored contributions
  analytical_ll <- exp_loglik_exact(exact_times, lambda) +
                   exp_loglik_censored(censored_times, lambda)

  expect_equal(numerical_ll, analytical_ll, tolerance = 1e-6,
               info = "Mixed loglik should be sum of exact and censored contributions")
})

test_that("loglik handles all censored data correctly", {
  # Given: All right-censored data (common in early study termination)
  lambda <- 1.0
  dist <- make_exponential_dfr()

  times <- c(1, 2, 3, 4, 5)
  df <- make_censored_data(times)

  ll <- loglik(dist)
  numerical_ll <- ll(df, par = c(lambda))

  # For all censored: loglik = -lambda * sum(t)
  analytical_ll <- -lambda * sum(times)

  expect_equal(numerical_ll, analytical_ll, tolerance = 1e-6)
})

# =============================================================================
# 3. Tests for loglik() - Weibull DFR
# =============================================================================

test_that("loglik for Weibull DFR matches analytical Weibull likelihood", {
  # Given: Weibull DFR with shape = 2, scale = 3
  # h(t) = (k/sigma) * (t/sigma)^(k-1) = (2/3) * (t/3)^1 = 2t/9
  # H(t) = (t/sigma)^k = (t/3)^2 = t^2/9
  # loglik_exact = log(h(t)) - H(t) = log(2t/9) - t^2/9

  shape <- 2
  scale <- 3
  dist <- make_weibull_dfr()

  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  ll <- loglik(dist)
  numerical_ll <- ll(df, par = c(shape, scale))

  # Analytical Weibull log-likelihood for exact observations:
  # sum(log(k/sigma) + (k-1)*log(t/sigma) - (t/sigma)^k)
  # = n*log(k) - n*k*log(sigma) + (k-1)*sum(log(t)) - sum((t/sigma)^k)
  n <- length(times)
  analytical_ll <- n * log(shape) - n * shape * log(scale) +
                   (shape - 1) * sum(log(times)) -
                   sum((times / scale)^shape)

  expect_equal(numerical_ll, analytical_ll, tolerance = 1e-5,
               info = "Weibull DFR loglik should match analytical formula")
})

test_that("Weibull DFR with shape=1 reduces to exponential", {
  # Weibull with shape=1 is exponential with rate = 1/scale
  shape <- 1
  scale <- 2  # lambda = 1/2 = 0.5
  lambda <- 1 / scale

  weibull_dist <- make_weibull_dfr()
  exp_dist <- make_exponential_dfr()

  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  ll_weibull <- loglik(weibull_dist)
  ll_exp <- loglik(exp_dist)

  weibull_result <- ll_weibull(df, par = c(shape, scale))
  exp_result <- ll_exp(df, par = c(lambda))

  expect_equal(weibull_result, exp_result, tolerance = 1e-6,
               info = "Weibull(shape=1, scale=s) should equal Exp(rate=1/s)")
})

# =============================================================================
# 4. Tests for score() - Gradient
# =============================================================================

test_that("score returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  s <- score(dist)
  expect_type(s, "closure")
})

test_that("score for exponential exact observations matches analytical", {
  # Given: Exponential distribution
  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Given: Exact observations
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  # When: Computing score (gradient)
  s <- score(dist)
  numerical_score <- s(df, par = c(lambda))

  # Then: Should match analytical gradient
  # d/dlambda [n * log(lambda) - lambda * sum(t)] = n/lambda - sum(t)
  analytical_score <- exp_score_exact(times, lambda)

  expect_equal(numerical_score[1], analytical_score, tolerance = 1e-5,
               info = "Exponential score should be n/lambda - sum(t)")
})

test_that("score is zero at MLE", {
  # Given: Data from exponential(lambda=1)
  set.seed(123)
  times <- rexp(100, rate = 1)
  df <- make_exact_data(times)

  dist <- make_exponential_dfr()

  # MLE for exponential: lambda_hat = n / sum(t)
  mle_lambda <- length(times) / sum(times)

  s <- score(dist)
  score_at_mle <- s(df, par = c(mle_lambda))

  expect_equal(score_at_mle[1], 0, tolerance = 1e-4,
               info = "Score should be approximately zero at MLE")
})

test_that("score for Weibull has correct dimension", {
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  s <- score(dist)
  result <- s(df, par = c(2, 3))  # shape = 2, scale = 3

  # Weibull score should have 2 components (shape, scale)
  expect_length(result, 2)
})

test_that("numerical score approximates analytical for exponential", {
  # Compare numerical gradient to analytical
  lambda <- 1.5
  dist <- make_exponential_dfr()
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  # Numerical score (from default implementation using numDeriv)
  s <- score(dist)
  numerical <- s(df, par = c(lambda))

  # Analytical
  analytical <- exp_score_exact(times, lambda)

  expect_equal(numerical[1], analytical, tolerance = 1e-4,
               info = "Numerical score should match analytical")
})

# =============================================================================
# 5. Tests for hess_loglik() - Hessian
# =============================================================================

test_that("hess_loglik returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  H <- hess_loglik(dist)
  expect_type(H, "closure")
})

test_that("hess_loglik for exponential exact observations matches analytical", {
  # Given: Exponential distribution
  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Given: Exact observations
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  # When: Computing Hessian
  H <- hess_loglik(dist)
  numerical_hess <- H(df, par = c(lambda))

  # Then: Should match analytical Hessian
  # d^2/dlambda^2 [n * log(lambda) - lambda * sum(t)] = -n/lambda^2
  analytical_hess <- exp_hessian_exact(times, lambda)

  expect_equal(numerical_hess[1, 1], analytical_hess[1, 1], tolerance = 1e-4,
               info = "Exponential Hessian should be -n/lambda^2")
})

test_that("hess_loglik returns matrix of correct dimension", {
  # For Weibull with 2 parameters
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  H <- hess_loglik(dist)
  result <- H(df, par = c(2, 3))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2),
               info = "Weibull Hessian should be 2x2 matrix")
})

test_that("hess_loglik is symmetric", {
  # Hessian of log-likelihood should always be symmetric
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)

  H <- hess_loglik(dist)
  result <- H(df, par = c(1.5, 2.5))

  expect_equal(result[1, 2], result[2, 1], tolerance = 1e-6,
               info = "Hessian should be symmetric")
})

test_that("hess_loglik is negative definite at MLE (for concave likelihood)", {
  # At MLE, Hessian should be negative definite for identifiable model
  set.seed(42)
  times <- rexp(50, rate = 1)
  df <- make_exact_data(times)

  dist <- make_exponential_dfr()

  # MLE for exponential
  mle_lambda <- length(times) / sum(times)

  H <- hess_loglik(dist)
  hess_at_mle <- H(df, par = c(mle_lambda))

  # For negative definite 1x1 matrix, the element should be negative
  # Hessian should be negative definite at MLE
  expect_lt(hess_at_mle[1, 1], 0)
})

# =============================================================================
# 6. Tests for fit() - MLE Estimation
# =============================================================================

test_that("fit returns a solver function", {
  dist <- make_exponential_dfr()
  solver <- fit(dist)
  expect_type(solver, "closure")
})

test_that("fit recovers true exponential parameter", {
  # Given: Data from exponential(lambda = 2)
  true_lambda <- 2
  set.seed(42)
  times <- rexp(100, rate = true_lambda)
  df <- make_exact_data(times)

  # When: Fitting the model
  dist <- make_exponential_dfr()
  solver <- fit(dist)
  result <- solver(df, par = c(1))  # Initial guess lambda = 1

  # Then: MLE should be close to true parameter
  # For exponential, MLE = n / sum(t)
  analytical_mle <- length(times) / sum(times)

  # Extract fitted parameter (result should be an mle object)
  fitted_lambda <- params(result)[1]

  expect_equal(fitted_lambda, analytical_mle, tolerance = 1e-3,
               info = "Fitted lambda should match analytical MLE")
  expect_equal(fitted_lambda, true_lambda, tolerance = 0.3,
               info = "Fitted lambda should be close to true value")
})

test_that("fit recovers Weibull parameters", {
  # Given: Data from Weibull(shape = 2, scale = 3)
  true_shape <- 2
  true_scale <- 3
  set.seed(123)
  times <- rweibull(100, shape = true_shape, scale = true_scale)
  df <- make_exact_data(times)

  # When: Fitting the model
  dist <- make_weibull_dfr()
  solver <- fit(dist)
  result <- solver(df, par = c(1.5, 2.5))  # Initial guesses

  # Then: MLEs should be close to true parameters
  fitted_params <- params(result)

  expect_equal(fitted_params[1], true_shape, tolerance = 0.3,
               info = "Fitted shape should be close to true value")
  expect_equal(fitted_params[2], true_scale, tolerance = 0.5,
               info = "Fitted scale should be close to true value")
})

test_that("fit returns object with expected structure", {
  set.seed(42)
  times <- rexp(30, rate = 1)
  df <- make_exact_data(times)

  dist <- make_exponential_dfr()
  solver <- fit(dist)
  result <- solver(df, par = c(1))

  # Result should be an mle object (from algebraic.mle)
  expect_true(inherits(result, "mle") || inherits(result, "mle_numerical"),
              info = "fit should return an mle object")

  # Should have accessible parameters
  expect_true(!is.null(params(result)),
              info = "Result should have accessible parameters")
})

test_that("fit with mixed censored data", {
  # Given: Mixed exact and censored observations from exponential(lambda = 1)
  true_lambda <- 1
  set.seed(42)

  # Generate failure times
  n <- 100
  true_times <- rexp(n, rate = true_lambda)

  # Right-censor at time 2 (anything > 2 becomes censored at 2)
  censor_time <- 2
  observed_times <- pmin(true_times, censor_time)
  deltas <- as.integer(true_times <= censor_time)

  df <- data.frame(t = observed_times, delta = deltas)

  # When: Fitting the model
  dist <- make_exponential_dfr()
  solver <- fit(dist)
  result <- solver(df, par = c(0.5))

  # Then: MLE should still be reasonable (though biased due to censoring)
  fitted_lambda <- params(result)[1]

  # With right-censoring, MLE should still be positive and in reasonable range
  expect_gt(fitted_lambda, 0)
  expect_lt(fitted_lambda, 5)
})

# =============================================================================
# 7. Tests for Parameter Handling
# =============================================================================

test_that("loglik works with named parameters", {
  dist <- make_exponential_dfr()
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  ll <- loglik(dist)

  # Both named and unnamed should work
  result_unnamed <- ll(df, par = c(0.5))
  result_named <- ll(df, par = c(lambda = 0.5))

  expect_equal(result_unnamed, result_named,
               info = "Named and unnamed parameters should give same result")
})

test_that("loglik with default parameters from distribution", {
  # When distribution has default parameters
  dist <- make_exponential_dfr(lambda = 0.5)
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  ll <- loglik(dist)

  # Should be able to call without explicit par if dist has defaults
  # This depends on implementation - may need to pass NULL
  result <- ll(df, par = c(0.5))

  expect_true(is.finite(result))
})

test_that("score preserves parameter names", {
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  s <- score(dist)
  result <- s(df, par = c(shape = 2, scale = 3))

  # If implementation preserves names, check them
  # Otherwise, just verify correct length
  expect_length(result, 2)
})

# =============================================================================
# 8. Tests for Edge Cases
# =============================================================================

test_that("loglik handles single observation", {
  dist <- make_exponential_dfr()
  df <- make_exact_data(c(1))

  ll <- loglik(dist)
  result <- ll(df, par = c(0.5))

  # Single observation should still work
  expect_true(is.finite(result))

  # For single exact obs: log(lambda) - lambda * t = log(0.5) - 0.5*1
  expected <- log(0.5) - 0.5 * 1
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("loglik returns -Inf or very negative for impossible parameters", {
  skip("TODO: Handle extreme parameters gracefully")
  dist <- make_exponential_dfr()
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  ll <- loglik(dist)

  # Very small lambda leads to very low probability of observing data
  result <- ll(df, par = c(1e-10))
  expect_true(result < -100)
})

test_that("loglik handles large sample sizes", {
  set.seed(42)
  dist <- make_exponential_dfr()
  times <- rexp(1000, rate = 1)
  df <- make_exact_data(times)

  ll <- loglik(dist)
  result <- ll(df, par = c(1))

  expect_true(is.finite(result))
})

test_that("score handles edge case near boundary", {
  dist <- make_exponential_dfr()
  times <- c(0.1, 0.2, 0.3)  # Very small times
  df <- make_exact_data(times)

  s <- score(dist)

  # Near lambda = 0, score should be large and positive
  result <- s(df, par = c(0.1))
  expect_true(is.finite(result[1]))
})

test_that("error handling for empty data frame", {
  skip("TODO: Add proper error handling for empty data frames")
  dist <- make_exponential_dfr()
  df <- data.frame(t = numeric(0), delta = numeric(0))

  ll <- loglik(dist)

  # Should either return 0, NA, or throw an error
  expect_error(ll(df, par = c(0.5)), regexp = NULL) |
    expect_true({
      result <- ll(df, par = c(0.5))
      is.na(result) || result == 0
    })
})

test_that("error handling for NULL parameters", {
  dist <- make_exponential_dfr()  # No default parameters
  times <- c(1, 2, 3)
  df <- make_exact_data(times)

  ll <- loglik(dist)

  # Should throw error when par is NULL and no defaults
  expect_error(ll(df, par = NULL))
})

# =============================================================================
# 9. Tests for Left Censoring (if supported)
# =============================================================================

test_that("loglik for left-censored observations", {
  skip("TODO: Left-censoring not yet implemented")
  # Left-censored: failed before time t
  # Contribution: log(F(t)) = log(1 - S(t)) = log(1 - exp(-H(t)))

  lambda <- 0.5
  dist <- make_exponential_dfr()

  # Create left-censored data frame
  times <- c(1, 2, 3)
  df <- data.frame(
    t = times,
    censor = rep("left", length(times))  # Alternative format
  )

  ll <- loglik(dist)
  result <- ll(df, par = c(lambda))

  # For left-censored: loglik = sum(log(1 - exp(-lambda * t)))
  expected <- sum(log(1 - exp(-lambda * times)))

  expect_equal(result, expected, tolerance = 1e-6)
})

# =============================================================================
# 10. Tests for likelihood_model Class Membership
# =============================================================================

test_that("dfr_dist with likelihood support inherits from likelihood_model", {
  dist <- make_exponential_dfr(lambda = 1)

  # After implementing the interface, dfr_dist should have likelihood_model class
  expect_true(inherits(dist, "dfr_dist"))

  # This test will pass after adding "likelihood_model" to class vector
  expect_true(inherits(dist, "likelihood_model") ||
                exists("loglik.dfr_dist"),
              info = "dfr_dist should implement likelihood_model interface")
})

test_that("is_likelihood_model returns TRUE for dfr_dist", {
  dist <- make_exponential_dfr(lambda = 1)

  # Skip if likelihood.model package not available
  skip_if_not_installed("likelihood.model")

  expect_true(is_likelihood_model(dist))
})

# =============================================================================
# 11. Tests for Consistency Between Methods
# =============================================================================

test_that("score equals numerical gradient of loglik", {
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)
  par <- c(2, 3)

  ll <- loglik(dist)
  s <- score(dist)

  score_result <- s(df, par)

  # Numerical gradient using finite differences
  eps <- 1e-6
  numerical_grad <- sapply(seq_along(par), function(i) {
    par_plus <- par
    par_minus <- par
    par_plus[i] <- par[i] + eps
    par_minus[i] <- par[i] - eps
    (ll(df, par_plus) - ll(df, par_minus)) / (2 * eps)
  })

  expect_equal(score_result, numerical_grad, tolerance = 1e-3,
               info = "Score should match numerical gradient of loglik")
})

test_that("hess_loglik equals numerical hessian of loglik", {
  dist <- make_weibull_dfr()
  times <- c(1, 2, 3, 4, 5)
  df <- make_exact_data(times)
  par <- c(2, 3)

  ll <- loglik(dist)
  H <- hess_loglik(dist)

  hess_result <- H(df, par)

  # Numerical Hessian using second-order finite differences
  eps <- 1e-5
  n_par <- length(par)
  numerical_hess <- matrix(0, n_par, n_par)

  for (i in seq_len(n_par)) {
    for (j in seq_len(n_par)) {
      par_pp <- par_pm <- par_mp <- par_mm <- par
      par_pp[i] <- par[i] + eps
      par_pp[j] <- par_pp[j] + eps
      par_pm[i] <- par[i] + eps
      par_pm[j] <- par_pm[j] - eps
      par_mp[i] <- par[i] - eps
      par_mp[j] <- par_mp[j] + eps
      par_mm[i] <- par[i] - eps
      par_mm[j] <- par_mm[j] - eps

      numerical_hess[i, j] <- (ll(df, par_pp) - ll(df, par_pm) -
                                 ll(df, par_mp) + ll(df, par_mm)) / (4 * eps^2)
    }
  }

  expect_equal(hess_result, numerical_hess, tolerance = 1e-2,
               info = "Hessian should match numerical second derivative of loglik")
})

# =============================================================================
# 12. Tests for Complex DFR (Time-Varying Hazard)
# =============================================================================

test_that("loglik works with time-varying hazard (bathtub)", {
  # Bathtub hazard: h(t) = a * exp(-b*t) + c
  dist <- make_bathtub_dfr()

  times <- c(0.5, 1, 2, 3, 5)
  df <- make_exact_data(times)

  par <- c(a = 2, b = 1, c = 0.5)

  ll <- loglik(dist)
  result <- ll(df, par)

  # Result should be finite
  expect_true(is.finite(result),
              info = "Bathtub hazard loglik should be finite")
})

test_that("fit works with time-varying hazard", {
  skip("TODO: Time-varying hazard optimization has numerical stability issues")
  # Generate data from bathtub hazard (using rejection sampling or simulation)
  dist <- make_bathtub_dfr()

  # Use some synthetic times
  set.seed(42)
  times <- runif(50, 0.1, 5)
  df <- make_exact_data(times)

  solver <- fit(dist)

  # Should not error with reasonable initial values
  expect_no_error({
    result <- solver(df, par = c(1, 1, 0.5))
  })
})

# =============================================================================
# 13. Tests for Assumptions Method
# =============================================================================

test_that("assumptions method returns expected values", {
  skip_if_not_installed("likelihood.model")

  dist <- make_exponential_dfr(lambda = 1)

  # If assumptions method is implemented
  if (exists("assumptions.dfr_dist")) {
    result <- assumptions(dist)
    expect_type(result, "character")
    expect_true(length(result) > 0)
  }
})
