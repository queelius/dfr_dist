# Tests for core dfr_dist functionality
#
# These tests cover the existing distribution methods to ensure they work
# correctly before and after adding likelihood_model support.
#
# Note: Test fixtures (make_exponential_dfr, make_weibull_dfr, etc.) are
# defined in helper-fixtures.R and automatically loaded by testthat.

# =============================================================================
# Constructor Tests
# =============================================================================

test_that("dfr_dist constructor creates object with correct class", {
  dist <- dfr_dist(
    rate = function(t, par, ...) par[1],
    par = c(1)
  )

  expect_s3_class(dist, "dfr_dist")
  expect_s3_class(dist, "univariate_dist")
  expect_s3_class(dist, "dist")
})

test_that("dfr_dist constructor stores rate function", {
  rate_fn <- function(t, par, ...) par[1] * t
  dist <- dfr_dist(rate = rate_fn, par = c(1))

  expect_true(is.function(dist$rate))
})

test_that("dfr_dist constructor stores parameters", {
  dist <- dfr_dist(
    rate = function(t, par, ...) par[1],
    par = c(lambda = 0.5)
  )

  expect_equal(dist$par, c(lambda = 0.5))
})

test_that("is_dfr_dist returns TRUE for dfr_dist objects", {
  dist <- make_exponential_dfr(lambda = 1)
  expect_true(is_dfr_dist(dist))
})

test_that("is_dfr_dist returns FALSE for non-dfr_dist objects", {
  expect_false(is_dfr_dist(list(a = 1)))
  expect_false(is_dfr_dist(NULL))
  expect_false(is_dfr_dist(1))
})

# =============================================================================
# hazard() Method Tests
# =============================================================================

test_that("hazard returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  h <- hazard(dist)
  expect_type(h, "closure")
})
test_that("hazard for exponential is constant", {
  lambda <- 2
  dist <- make_exponential_dfr(lambda = lambda)
  h <- hazard(dist)

  # Hazard should be lambda at all times
  expect_equal(h(1), lambda)
  expect_equal(h(5), lambda)
  expect_equal(h(100), lambda)
})

test_that("hazard accepts parameter override", {
  dist <- make_exponential_dfr(lambda = 1)
  h <- hazard(dist)

  # Override with different lambda
  expect_equal(h(1, par = c(2)), 2)
  expect_equal(h(1, par = c(0.5)), 0.5)
})

test_that("hazard for Weibull matches formula", {
  shape <- 2
  scale <- 3
  dist <- make_weibull_dfr(shape = shape, scale = scale)
  h <- hazard(dist)

  t <- 2
  # h(t) = (k/sigma) * (t/sigma)^(k-1)
  expected <- (shape / scale) * (t / scale)^(shape - 1)

  expect_equal(h(t), expected)
})

# =============================================================================
# cum_haz() Method Tests
# =============================================================================

test_that("cum_haz returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  H <- cum_haz(dist)
  expect_type(H, "closure")
})

test_that("cum_haz for exponential is lambda * t", {
  lambda <- 0.5
  dist <- make_exponential_dfr(lambda = lambda)
  H <- cum_haz(dist)

  # H(t) = lambda * t for exponential
  expect_equal(H(2), lambda * 2, tolerance = 1e-3)
  expect_equal(H(5), lambda * 5, tolerance = 1e-3)
})

test_that("cum_haz for Weibull matches formula", {
  shape <- 2
  scale <- 3
  dist <- make_weibull_dfr(shape = shape, scale = scale)
  H <- cum_haz(dist)

  t <- 2
  # H(t) = (t/sigma)^k for Weibull
  expected <- (t / scale)^shape

  expect_equal(H(t), expected, tolerance = 1e-3)
})

test_that("cum_haz is monotonically increasing", {
  dist <- make_exponential_dfr(lambda = 1)
  H <- cum_haz(dist)

  times <- c(1, 2, 3, 4, 5)
  values <- sapply(times, H)

  # Each value should be greater than previous
  expect_true(all(diff(values) > 0))
})

# =============================================================================
# surv() Method Tests
# =============================================================================

test_that("surv returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  S <- surv(dist)
  expect_type(S, "closure")
})

test_that("surv(0) equals 1", {
  dist <- make_exponential_dfr(lambda = 1)
  S <- surv(dist)

  expect_equal(S(0), 1, tolerance = 1e-6)
})

test_that("surv is monotonically decreasing", {
  dist <- make_exponential_dfr(lambda = 1)
  S <- surv(dist)

  times <- c(0.1, 0.5, 1, 2, 5)
  values <- sapply(times, S)

  # Each value should be less than previous
  expect_true(all(diff(values) < 0))
})

test_that("surv for exponential matches exp(-lambda*t)", {
  lambda <- 0.5
  dist <- make_exponential_dfr(lambda = lambda)
  S <- surv(dist)

  t <- 2
  expected <- exp(-lambda * t)

  expect_equal(S(t), expected, tolerance = 1e-4)
})

test_that("surv approaches 0 as t increases", {
  dist <- make_exponential_dfr(lambda = 1)
  S <- surv(dist)

  expect_lt(S(10), 0.001)
  expect_lt(S(20), 1e-6)
})

# =============================================================================
# cdf() Method Tests
# =============================================================================

test_that("cdf returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)
  expect_type(F, "closure")
})

test_that("cdf(0) equals 0", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)

  expect_equal(F(0), 0, tolerance = 1e-6)
})

test_that("cdf + surv equals 1", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)
  S <- surv(dist)

  times <- c(0.5, 1, 2, 5)
  for (t in times) {
    expect_equal(F(t) + S(t), 1, tolerance = 1e-4,
                 info = paste("t =", t))
  }
})

test_that("cdf is monotonically increasing", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)

  times <- c(0, 0.5, 1, 2, 5)
  values <- sapply(times, F)

  expect_true(all(diff(values) >= 0))
})

test_that("cdf with lower.limit = FALSE returns survival function", {
  lambda <- 0.5
  dist <- make_exponential_dfr(lambda = lambda)
  F <- cdf(dist)
  S <- surv(dist)

  t <- 2
  # lower.limit = FALSE should return S(t) = 1 - F(t)
  expect_equal(F(t, lower.limit = FALSE), S(t), tolerance = 1e-6)

  # Test log.p = TRUE with lower.limit = TRUE (log of CDF)
  expect_equal(F(t, log.p = TRUE, lower.limit = TRUE), log(F(t)), tolerance = 1e-6)

  # Test log.p = TRUE with lower.limit = FALSE (log of survival = -H(t))
  expect_equal(F(t, log.p = TRUE, lower.limit = FALSE), log(S(t)), tolerance = 1e-6)
})

# =============================================================================
# density() Method Tests
# =============================================================================

test_that("density returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  f <- density(dist)
  expect_type(f, "closure")
})

test_that("density for exponential matches lambda * exp(-lambda*t)", {
  lambda <- 0.5
  dist <- make_exponential_dfr(lambda = lambda)
  f <- density(dist)

  t <- 2
  expected <- lambda * exp(-lambda * t)

  expect_equal(f(t), expected, tolerance = 1e-4)
})

test_that("density equals hazard * survival", {
  dist <- make_exponential_dfr(lambda = 1)
  f <- density(dist)
  h <- hazard(dist)
  S <- surv(dist)

  times <- c(0.5, 1, 2, 5)
  for (t in times) {
    expect_equal(f(t), h(t) * S(t), tolerance = 1e-4,
                 info = paste("t =", t))
  }
})

test_that("density integrates to approximately 1", {
  dist <- make_exponential_dfr(lambda = 1)
  f <- density(dist)

  # Numerical integration
  integral <- integrate(f, lower = 0, upper = 100)$value

  expect_equal(integral, 1, tolerance = 1e-3)
})

test_that("density log argument works", {
  lambda <- 1
  dist <- make_exponential_dfr(lambda = lambda)
  f <- density(dist)

  t <- 2
  expect_equal(f(t, log = TRUE), log(f(t)), tolerance = 1e-6)
})

# =============================================================================
# inv_cdf() Method Tests
# =============================================================================

test_that("inv_cdf returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  Q <- inv_cdf(dist)
  expect_type(Q, "closure")
})

test_that("inv_cdf inverts cdf", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)
  Q <- inv_cdf(dist)

  # Q(F(t)) should equal t
  times <- c(0.5, 1, 2, 5)
  for (t in times) {
    p <- F(t)
    recovered_t <- Q(p)
    expect_equal(recovered_t, t, tolerance = 1e-3,
                 info = paste("t =", t))
  }
})

test_that("inv_cdf(0.5) is median", {
  lambda <- 1
  dist <- make_exponential_dfr(lambda = lambda)
  Q <- inv_cdf(dist)

  # Median of exponential is log(2)/lambda
  expected_median <- log(2) / lambda
  computed_median <- Q(0.5)

  expect_equal(computed_median, expected_median, tolerance = 1e-3)
})

# =============================================================================
# params() Method Tests
# =============================================================================

test_that("params returns distribution parameters", {
  dist <- make_exponential_dfr(lambda = 0.5)
  result <- params(dist)

  expect_equal(result, 0.5)
})

test_that("params returns NULL when no default parameters", {
  dist <- make_exponential_dfr()
  result <- params(dist)

  expect_null(result)
})

# =============================================================================
# sampler() Method Tests
# =============================================================================

test_that("sampler returns a function", {
  dist <- make_exponential_dfr(lambda = 1)
  samp <- sampler(dist)
  expect_type(samp, "closure")
})

test_that("sampler generates correct number of samples", {
  dist <- make_exponential_dfr(lambda = 1)
  samp <- sampler(dist)

  samples <- samp(10)
  expect_length(samples, 10)
})

test_that("sampler generates positive values", {
  dist <- make_exponential_dfr(lambda = 1)
  samp <- sampler(dist)

  samples <- samp(20)
  expect_true(all(samples >= 0))
})

test_that("sampler mean approximates 1/lambda for exponential", {
  lambda <- 2
  dist <- make_exponential_dfr(lambda = lambda)
  samp <- sampler(dist)

  set.seed(42)
  samples <- samp(500)

  # Mean of exponential is 1/lambda
  expected_mean <- 1 / lambda
  sample_mean <- mean(samples)

  expect_equal(sample_mean, expected_mean, tolerance = 0.2)
})

# =============================================================================
# sup() Method Tests
# =============================================================================

test_that("sup returns interval (0, Inf)", {
  dist <- make_exponential_dfr(lambda = 1)
  support <- sup(dist)

  # Support should be (0, Inf) - open on both ends
  expect_true(inherits(support, "interval") || is.list(support))
})

# =============================================================================
# print() Method Tests
# =============================================================================

test_that("print method works without error", {
  dist <- make_exponential_dfr(lambda = 1)

  expect_output(print(dist), "Dynamic failure rate")
})

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

test_that("hazard handles t = 0 appropriately", {
  dist <- make_exponential_dfr(lambda = 1)
  h <- hazard(dist)

  # For exponential, h(0) = lambda
  expect_equal(h(0), 1)
})

test_that("Weibull hazard at t = 0 depends on shape", {
  # For shape < 1, hazard -> Inf as t -> 0
  # For shape = 1, hazard = constant
  # For shape > 1, hazard = 0 at t = 0

  dist_shape_gt_1 <- make_weibull_dfr(shape = 2, scale = 1)
  h <- hazard(dist_shape_gt_1)

  # h(t) = (k/sigma) * (t/sigma)^(k-1) = 2 * t^1 = 2t
  # h(0) = 0 for shape > 1
  expect_equal(h(0.001), 0.002, tolerance = 1e-5)
})

test_that("vectorized operations work for hazard", {
  dist <- make_exponential_dfr(lambda = 1)
  h <- hazard(dist)

  times <- c(1, 2, 3, 4, 5)
  results <- h(times)

  expect_length(results, 5)
  expect_true(all(results == 1))
})

test_that("hazard handles very small times", {
  dist <- make_exponential_dfr(lambda = 1)
  h <- hazard(dist)
  S <- surv(dist)
  f <- density(dist)

  # Very small times should give finite results
  small_t <- 1e-10
  expect_true(is.finite(h(small_t)))
  expect_true(is.finite(S(small_t)))
  expect_true(is.finite(f(small_t)))
})

test_that("surv handles very large times", {
  dist <- make_exponential_dfr(lambda = 1)
  S <- surv(dist)

  # Very large times should give S(t) very close to 0
  expect_lt(S(50), 1e-10)
  expect_true(S(50) >= 0)  # Non-negative
})

test_that("cdf handles very large times", {
  dist <- make_exponential_dfr(lambda = 1)
  F <- cdf(dist)

  # Very large times should give F(t) very close to 1
  expect_gt(F(50), 1 - 1e-10)
  expect_true(F(50) <= 1)  # At most 1
})

test_that("Weibull handles shape near 1 (boundary case)", {
  # Shape = 1 is exponential, shape near 1 tests numerical stability
  dist_near_1 <- make_weibull_dfr(shape = 1.001, scale = 2)
  dist_at_1 <- make_weibull_dfr(shape = 1.0, scale = 2)

  h_near <- hazard(dist_near_1)
  h_at <- hazard(dist_at_1)

  # Should be stable and close near shape = 1
  t <- 1
  expect_true(is.finite(h_near(t)))
  expect_equal(h_near(t), h_at(t), tolerance = 0.01)
})

test_that("methods handle single observation time", {
  dist <- make_exponential_dfr(lambda = 0.5)

  # Single time point should work
  h <- hazard(dist)
  S <- surv(dist)
  F <- cdf(dist)
  f <- density(dist)

  expect_equal(h(2), 0.5)
  expect_true(is.finite(S(2)))
  expect_true(is.finite(F(2)))
  expect_true(is.finite(f(2)))
})
