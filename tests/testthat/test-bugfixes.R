# Tests for bug fixes from left-censoring and vectorization review (v0.5.0)

# =============================================================================
# Cumulative hazard numerical path: vector t input
# =============================================================================

test_that("cum_haz numerical path accepts vector t", {
    # No cum_haz_rate forces numerical integration path
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        par = c(0.5)
    )
    H <- cum_haz(dist)

    t_vec <- c(1, 2, 3, 4, 5)
    result <- H(t_vec)

    expect_length(result, 5)

    for (i in seq_along(t_vec)) {
        expect_equal(result[i], H(t_vec[i]), tolerance = 1e-3,
                     info = paste("t =", t_vec[i]))
    }

    # Constant rate 0.5: H(t) = 0.5 * t
    expect_equal(result, 0.5 * t_vec, tolerance = 1e-3)
})

test_that("cum_haz numerical path matches analytical for vectors", {
    lambda <- 0.5

    # Analytical path
    dist_analytical <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        par = c(lambda)
    )

    # Numerical path (no cum_haz_rate)
    dist_numerical <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        par = c(lambda)
    )

    H_analytical <- cum_haz(dist_analytical)
    H_numerical <- cum_haz(dist_numerical)

    t_vec <- c(0.5, 1, 2, 5, 10)
    expect_equal(H_numerical(t_vec), H_analytical(t_vec), tolerance = 1e-3)
})

# =============================================================================
# Score/Hessian with custom ob_col and delta_col
# =============================================================================

test_that("dfr_exponential score works with custom ob_col", {
    # Use "obs" (not "time") to avoid R's $ partial matching with "t"
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        ob_col = "obs",
        score_fn = dfr_exponential()$score_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    s <- score(dist)
    result <- s(df, par = c(0.5))

    expect_true(is.finite(result))
    # n_events/lambda - sum(obs) = 3/0.5 - 15 = -9
    expect_equal(result, c(3 / 0.5 - sum(df$obs)), tolerance = 1e-6)
})

test_that("dfr_exponential hessian works with custom ob_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        ob_col = "obs",
        hess_fn = dfr_exponential()$hess_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    H <- hess_loglik(dist)
    result <- H(df, par = c(0.5))

    # -n_events / lambda^2 = -3 / 0.25 = -12
    expect_equal(result[1, 1], -12, tolerance = 1e-6)
})

test_that("dfr_exponential score works with custom delta_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        delta_col = "status",
        score_fn = dfr_exponential()$score_fn
    )

    df <- data.frame(t = c(1, 2, 3, 4, 5), status = c(1, 1, 1, 0, 0))
    s <- score(dist)

    result <- s(df, par = c(0.5))
    expect_true(is.finite(result))

    # n_events=3, lambda=0.5, sum_t=15 => 3/0.5 - 15 = -9
    expect_equal(result, c(-9), tolerance = 1e-6)
})

test_that("dfr_weibull score works with custom ob_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) {
            k <- par[[1]]; sigma <- par[[2]]
            (k / sigma) * (t / sigma)^(k - 1)
        },
        cum_haz_rate = function(t, par, ...) (t / par[[2]])^par[[1]],
        ob_col = "obs",
        score_fn = dfr_weibull()$score_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    s <- score(dist)
    result <- s(df, par = c(2, 3))
    expect_length(result, 2)
    expect_true(all(is.finite(result)))
})

test_that("dfr_weibull hessian works with custom ob_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) {
            k <- par[[1]]; sigma <- par[[2]]
            (k / sigma) * (t / sigma)^(k - 1)
        },
        cum_haz_rate = function(t, par, ...) (t / par[[2]])^par[[1]],
        ob_col = "obs",
        hess_fn = dfr_weibull()$hess_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    H <- hess_loglik(dist)
    result <- H(df, par = c(2, 3))

    expect_equal(dim(result), c(2, 2))
    expect_true(all(is.finite(result)))
})

test_that("dfr_gompertz score works with custom ob_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) par[[1]] * exp(par[[2]] * t),
        cum_haz_rate = function(t, par, ...) (par[[1]] / par[[2]]) * (exp(par[[2]] * t) - 1),
        ob_col = "obs",
        score_fn = dfr_gompertz()$score_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    s <- score(dist)
    result <- s(df, par = c(0.1, 0.2))

    expect_length(result, 2)
    expect_true(all(is.finite(result)))
})

test_that("dfr_loglogistic score works with custom ob_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) {
            alpha <- par[[1]]; beta <- par[[2]]
            tr <- t / alpha
            (beta / alpha) * tr^(beta - 1) / (1 + tr^beta)
        },
        cum_haz_rate = function(t, par, ...) log(1 + (t / par[[1]])^par[[2]]),
        ob_col = "obs",
        score_fn = dfr_loglogistic()$score_fn
    )

    df <- data.frame(obs = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))
    s <- score(dist)
    result <- s(df, par = c(3, 2))

    expect_length(result, 2)
    expect_true(all(is.finite(result)))
})

test_that("score at MLE is zero with custom ob_col and delta_col", {
    set.seed(42)
    n <- 100
    times <- rexp(n, rate = 0.5)
    df <- data.frame(obs = times, status = rep(1, n))

    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        ob_col = "obs",
        delta_col = "status",
        score_fn = dfr_exponential()$score_fn,
        hess_fn = dfr_exponential()$hess_fn
    )

    s <- score(dist)
    mle <- n / sum(times)

    expect_equal(s(df, par = c(mle)), 0, tolerance = 1e-6)
})

# =============================================================================
# inv_cdf: vector p input
# =============================================================================

test_that("inv_cdf accepts vector of probabilities", {
    dist <- dfr_exponential(lambda = 1)
    Q <- inv_cdf(dist)

    p_vec <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    result <- Q(p_vec)

    # Should return a vector of same length
    expect_length(result, 5)

    # Each should match scalar call
    for (i in seq_along(p_vec)) {
        expect_equal(result[i], Q(p_vec[i]), tolerance = 1e-3,
                     info = paste("p =", p_vec[i]))
    }
})

test_that("inv_cdf vector result matches analytical quantiles", {
    lambda <- 2
    dist <- dfr_exponential(lambda = lambda)
    Q <- inv_cdf(dist)

    p_vec <- c(0.1, 0.5, 0.9)
    result <- Q(p_vec)

    # Exponential quantile: Q(p) = -log(1-p) / lambda
    expected <- -log(1 - p_vec) / lambda
    expect_equal(result, expected, tolerance = 1e-3)
})

test_that("inv_cdf vector works with Weibull", {
    dist <- dfr_weibull(shape = 2, scale = 3)
    Q <- inv_cdf(dist)

    p_vec <- c(0.25, 0.5, 0.75)
    result <- Q(p_vec)

    expect_length(result, 3)
    expect_true(all(diff(result) > 0))  # monotonically increasing
})

# =============================================================================
# loglik NaN guard for all censoring types
# =============================================================================

test_that("loglik returns -Inf (not NaN) for invalid params with right-censored data", {
    # cum_haz_rate includes log(par[1]) which produces NaN when par[1] < 0
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(exp(par[1]), length(t)),
        cum_haz_rate = function(t, par, ...) exp(par[1]) * t + log(par[1])
    )

    df <- data.frame(t = c(1, 2, 3), delta = c(0, 0, 0))
    ll <- loglik(dist)

    result <- suppressWarnings(ll(df, par = c(-1)))
    expect_equal(result, -Inf)
})

test_that("loglik returns -Inf (not NaN) for invalid params with left-censored data", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t
    )

    df <- data.frame(t = c(1, 2, 3), delta = c(-1, -1, -1))
    ll <- loglik(dist)

    # Negative lambda -> H(t) < 0 -> log1p(-exp(positive)) -> NaN -> -Inf
    result <- ll(df, par = c(-1))
    expect_equal(result, -Inf)
})

test_that("loglik returns -Inf for NaN in mixed censoring", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t
    )

    df <- data.frame(
        t = c(1, 2, 3, 4),
        delta = c(1, 0, -1, 1)
    )
    ll <- loglik(dist)
    result <- ll(df, par = c(-1))
    expect_equal(result, -Inf)
})

# =============================================================================
# DESCRIPTION mentions left-censoring support
# =============================================================================

test_that("DESCRIPTION mentions left-censoring support", {
    desc <- packageDescription("flexhaz")
    expect_true(
        grepl("left.censored", desc$Description, ignore.case = TRUE),
        info = "DESCRIPTION should mention left-censored data support"
    )
})
