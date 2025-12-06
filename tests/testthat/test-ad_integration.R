# Tests for AD integration with femtograd

test_that("exponential with analytical score gives correct gradient", {
    # Exponential: h(t) = lambda, H(t) = lambda * t
    # loglik = sum(delta) * log(lambda) - lambda * sum(t)
    # score = sum(delta) / lambda - sum(t)

    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        score_fn = function(df, par, ...) {
            t <- df$t
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, length(t))
            c(sum(delta) / par[1] - sum(t))
        }
    )

    # Test data
    set.seed(42)
    times <- rexp(20, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 20))

    # Get score function
    s <- score(exp_dist)

    # Test at MLE (score should be ~0)
    mle_lambda <- length(times) / sum(times)
    score_at_mle <- s(df, par = c(mle_lambda))
    expect_equal(score_at_mle, 0, tolerance = 1e-10)

    # Test at different parameter
    score_at_1 <- s(df, par = c(1))
    expected <- length(times) / 1 - sum(times)
    expect_equal(score_at_1, expected, tolerance = 1e-10)
})

test_that("exponential with analytical score gives correct Hessian via AD", {
    skip_if_not_installed("femtograd")

    # Exponential Hessian: d^2 loglik / d lambda^2 = -n / lambda^2
    # For AD Jacobian of score to work, score_fn must use femtograd ops

    n <- 20
    sum_t <- 15.5  # Fixed for reproducibility

    # Score function written to be AD-compatible
    # score = n/lambda - sum_t
    # d(score)/d(lambda) = -n/lambda^2
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        score_fn = function(df, par, ...) {
            # For AD compatibility, use scalar ops
            n_obs <- nrow(df)
            sum_t <- sum(df$t)
            c(n_obs / par[1] - sum_t)
        }
    )

    df <- data.frame(t = rep(sum_t / n, n), delta = rep(1, n))

    # Get Hessian function
    H <- hess_loglik(exp_dist)
    hess <- H(df, par = c(1))

    # Expected: -n / lambda^2 = -20 / 1 = -20
    expected_hess <- -n / 1^2
    expect_equal(hess[1, 1], expected_hess, tolerance = 1)
})

test_that("numerical fallback works when femtograd not available", {
    # Create distribution without analytical functions
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t))
    )

    set.seed(42)
    times <- rexp(20, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 20))

    # Score should still work via numerical differentiation
    s <- score(exp_dist)
    score_val <- s(df, par = c(1))

    # Should be approximately n/lambda - sum(t)
    expected <- 20 / 1 - sum(times)
    expect_equal(score_val, expected, tolerance = 0.01)

    # Hessian should also work
    H <- hess_loglik(exp_dist)
    hess <- H(df, par = c(1))

    # Should be approximately -n/lambda^2 = -20
    expect_equal(hess[1, 1], -20, tolerance = 0.5)
})

test_that("Weibull score matches numerical gradient", {
    # This test verifies numerical fallback works correctly for Weibull
    # No cum_haz_rate or score_fn -> uses numerical differentiation

    weibull_dist <- dfr_dist(
        rate = function(t, par, ...) {
            k <- par[1]
            sigma <- par[2]
            (k / sigma) * (t / sigma)^(k - 1)
        }
        # No cum_haz_rate -> forces numerical fallback
    )

    # Simulate data
    set.seed(123)
    true_k <- 2
    true_sigma <- 3
    u <- runif(30)
    times <- true_sigma * (-log(u))^(1/true_k)
    df <- data.frame(t = times, delta = rep(1, 30))

    # Compare score to numerical gradient
    ll <- loglik(weibull_dist)
    s <- score(weibull_dist)

    par_test <- c(1.5, 2.5)

    # Numerical gradient (ground truth)
    num_grad <- numDeriv::grad(function(p) ll(df, par = p), par_test)

    # Score function (numerical fallback)
    score_val <- s(df, par = par_test)

    expect_equal(score_val, num_grad, tolerance = 0.01)
})

test_that("ad_hessian helper function works correctly", {
    skip_if_not_installed("femtograd")

    # Simple quadratic: f(x) = x^2, gradient = 2x, Hessian = 2
    grad_fn <- function(par) {
        2 * par[1]
    }

    hess <- ad_hessian(grad_fn, c(3))
    expect_equal(hess[1, 1], 2, tolerance = 1e-10)
})

test_that("ad_jacobian helper function works correctly", {
    skip_if_not_installed("femtograd")

    # f(x, y) = (x + y, x * y)
    # Jacobian: [[1, 1], [y, x]]
    f <- function(par) {
        x <- par[[1]]
        y <- par[[2]]
        c(x + y, x * y)
    }

    jac <- ad_jacobian(f, c(2, 3))

    # At (2, 3): [[1, 1], [3, 2]]
    expect_equal(jac[1, 1], 1, tolerance = 1e-10)
    expect_equal(jac[1, 2], 1, tolerance = 1e-10)
    expect_equal(jac[2, 1], 3, tolerance = 1e-10)
    expect_equal(jac[2, 2], 2, tolerance = 1e-10)
})
