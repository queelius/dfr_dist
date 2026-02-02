# Tests for derivative computation (score and Hessian)
#
# These tests verify the 2-tier fallback chains:
#   Score:   score_fn → numDeriv::grad
#   Hessian: hess_fn  → numDeriv::hessian

# =============================================================================
# Score Function Tests
# =============================================================================

test_that("score uses analytical score_fn when provided", {
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        score_fn = function(df, par, ...) {
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
            c(sum(delta) / par[1] - sum(df$t))
        }
    )

    set.seed(42)
    times <- rexp(20, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 20))

    s <- score(exp_dist)

    # At MLE, score should be ~0
    mle_lambda <- length(times) / sum(times)
    expect_equal(s(df, par = c(mle_lambda)), 0, tolerance = 1e-10)

    # At different parameter
    score_at_1 <- s(df, par = c(1))
    expected <- 20 / 1 - sum(times)
    expect_equal(score_at_1, expected, tolerance = 1e-10)
})

test_that("score falls back to numDeriv::grad when no score_fn provided", {
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t))
    )

    set.seed(42)
    times <- rexp(20, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 20))

    s <- score(exp_dist)
    score_val <- s(df, par = c(1))

    # Should approximate analytical score: n/lambda - sum(t)
    expected <- 20 / 1 - sum(times)
    expect_equal(score_val, expected, tolerance = 0.01)
})

test_that("Weibull score via numDeriv matches analytical", {
    # Distribution with no score_fn → forces numerical fallback
    weibull_dist <- dfr_dist(
        rate = function(t, par, ...) {
            k <- par[1]
            sigma <- par[2]
            (k / sigma) * (t / sigma)^(k - 1)
        }
    )

    set.seed(123)
    true_k <- 2
    true_sigma <- 3
    u <- runif(30)
    times <- true_sigma * (-log(u))^(1/true_k)
    df <- data.frame(t = times, delta = rep(1, 30))

    ll <- loglik(weibull_dist)
    s <- score(weibull_dist)
    par_test <- c(1.5, 2.5)

    num_grad <- numDeriv::grad(function(p) ll(df, par = p), par_test)
    score_val <- s(df, par = par_test)

    expect_equal(score_val, num_grad, tolerance = 0.01)
})

# =============================================================================
# Hessian Function Tests
# =============================================================================

test_that("hess_loglik uses analytical hess_fn when provided", {
    n <- 20
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        hess_fn = function(df, par, ...) {
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
            n_events <- sum(delta == 1)
            matrix(-n_events / par[1]^2, nrow = 1, ncol = 1)
        }
    )

    df <- data.frame(t = rep(0.775, n), delta = rep(1, n))

    H <- hess_loglik(exp_dist)
    hess <- H(df, par = c(1))

    # Expected: -n / lambda^2 = -20 / 1 = -20
    expect_equal(hess[1, 1], -20, tolerance = 1e-10)
})

test_that("hess_loglik falls back to numDeriv::hessian when no hess_fn provided", {
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t))
    )

    set.seed(42)
    times <- rexp(20, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 20))

    H <- hess_loglik(exp_dist)
    hess <- H(df, par = c(1))

    # Should approximate -n/lambda^2 = -20
    expect_equal(hess[1, 1], -20, tolerance = 0.5)
})

# =============================================================================
# Combined score_fn + hess_fn Tests
# =============================================================================

test_that("score_fn and hess_fn both used when both provided", {
    # Track which functions are called
    score_called <- FALSE
    hess_called <- FALSE

    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        score_fn = function(df, par, ...) {
            score_called <<- TRUE
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
            c(sum(delta) / par[1] - sum(df$t))
        },
        hess_fn = function(df, par, ...) {
            hess_called <<- TRUE
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
            n_events <- sum(delta == 1)
            matrix(-n_events / par[1]^2, nrow = 1, ncol = 1)
        }
    )

    df <- data.frame(t = c(1, 2, 3), delta = c(1, 1, 1))

    s <- score(exp_dist)
    s(df, par = c(0.5))
    expect_true(score_called)

    H <- hess_loglik(exp_dist)
    H(df, par = c(0.5))
    expect_true(hess_called)
})

test_that("score_fn provided without hess_fn: score analytical, hess numerical", {
    exp_dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        score_fn = function(df, par, ...) {
            delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
            c(sum(delta) / par[1] - sum(df$t))
        }
        # No hess_fn → numerical fallback
    )

    set.seed(42)
    times <- rexp(50, rate = 1)
    df <- data.frame(t = times, delta = rep(1, 50))

    # Score should be exact
    s <- score(exp_dist)
    mle_lambda <- length(times) / sum(times)
    expect_equal(s(df, par = c(mle_lambda)), 0, tolerance = 1e-10)

    # Hessian should be numerical but correct
    H <- hess_loglik(exp_dist)
    hess <- H(df, par = c(1))
    expect_equal(hess[1, 1], -50, tolerance = 0.5)
})

# =============================================================================
# Weibull Multi-Parameter Tests
# =============================================================================

test_that("Weibull analytical Hessian matches numerical", {
    weibull_dist <- dfr_weibull()

    set.seed(42)
    n <- 50
    times <- 3 * (-log(runif(n)))^(1/2)
    df <- data.frame(t = times, delta = rep(1, n))

    test_par <- c(1.8, 2.8)

    # Analytical Hessian from dfr_weibull's hess_fn
    H <- hess_loglik(weibull_dist)
    hess_analytical <- H(df, par = test_par)

    # Numerical Hessian for comparison
    ll <- loglik(weibull_dist)
    hess_num <- numDeriv::hessian(function(p) ll(df, par = p), test_par)

    expect_equal(hess_analytical, hess_num, tolerance = 0.01)
})

# =============================================================================
# Log-Normal (Pure Numerical Fallback)
# =============================================================================

test_that("log-normal works with pure numerical fallback", {
    lognormal_dist <- dfr_dist(
        rate = function(t, par, ...) {
            mu <- par[1]
            sigma <- par[2]
            z <- (log(t) - mu) / sigma
            dnorm(z) / (t * sigma * pnorm(-z))
        }
    )

    set.seed(321)
    times <- rlnorm(50, meanlog = 1, sdlog = 0.5)
    df <- data.frame(t = times, delta = rep(1, 50))

    ll <- loglik(lognormal_dist)
    s <- score(lognormal_dist)
    H <- hess_loglik(lognormal_dist)

    ll_val <- ll(df, par = c(1, 0.5))
    expect_type(ll_val, "double")
    expect_true(is.finite(ll_val))

    score_val <- s(df, par = c(1, 0.5))
    expect_length(score_val, 2)

    hess_val <- H(df, par = c(1, 0.5))
    expect_equal(dim(hess_val), c(2, 2))
})

# =============================================================================
# Helper Constructor Derivative Tests
# =============================================================================

test_that("dfr_exponential has working score_fn and hess_fn", {
    dist <- dfr_exponential(lambda = 1)

    set.seed(123)
    times <- rexp(50, rate = 0.5)
    df <- data.frame(t = times, delta = rep(1, 50))

    # Score at MLE should be ~0
    s <- score(dist)
    mle_lambda <- length(times) / sum(times)
    expect_equal(s(df, par = c(mle_lambda)), 0, tolerance = 1e-10)

    # Hessian should be -n/lambda^2
    H <- hess_loglik(dist)
    hess <- H(df, par = c(1))
    expect_equal(hess[1, 1], -50, tolerance = 1e-10)
})

test_that("dfr_weibull has working score_fn and hess_fn", {
    dist <- dfr_weibull(shape = 2, scale = 3)

    set.seed(123)
    times <- rweibull(100, shape = 2, scale = 3)
    df <- data.frame(t = times, delta = rep(1, 100))

    s <- score(dist)
    score_val <- s(df, par = c(2, 3))
    expect_length(score_val, 2)

    H <- hess_loglik(dist)
    hess_val <- H(df, par = c(2, 3))
    expect_equal(dim(hess_val), c(2, 2))

    # Verify analytical Hessian matches numerical
    ll <- loglik(dist)
    hess_num <- numDeriv::hessian(function(p) ll(df, par = p), c(2, 3))
    expect_equal(hess_val, hess_num, tolerance = 0.01)
})

test_that("dfr_gompertz score_fn works, hess_fn uses numerical fallback", {
    dist <- dfr_gompertz(a = 0.01, b = 0.5)

    set.seed(456)
    a_true <- 0.01
    b_true <- 0.5
    u <- runif(50)
    times <- (1/b_true) * log(1 - (b_true/a_true) * log(1-u))
    times <- times[is.finite(times) & times > 0]
    df <- data.frame(t = times, delta = rep(1, length(times)))

    s <- score(dist)
    score_val <- s(df, par = c(a_true, b_true))
    expect_length(score_val, 2)

    H <- hess_loglik(dist)
    hess_val <- H(df, par = c(a_true, b_true))
    expect_equal(dim(hess_val), c(2, 2))
})

test_that("dfr_loglogistic score_fn works, hess_fn uses numerical fallback", {
    dist <- dfr_loglogistic(alpha = 2, beta = 3)

    set.seed(789)
    u <- runif(60)
    times <- 2 * ((1 - u) / u)^(1/3)
    df <- data.frame(t = times, delta = rep(1, length(times)))

    s <- score(dist)
    score_val <- s(df, par = c(2, 3))
    expect_length(score_val, 2)

    H <- hess_loglik(dist)
    hess_val <- H(df, par = c(2, 3))
    expect_equal(dim(hess_val), c(2, 2))
})
