# Tests for diagnostic methods (residuals, plot, qqplot)

test_that("residuals.dfr_dist computes Cox-Snell residuals correctly", {
    # For exponential with lambda, H(t) = lambda * t
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))

    cs_resid <- residuals(exp_dist, df, type = "cox-snell")

    # H(t) = 0.5 * t
    expected <- 0.5 * c(1, 2, 3, 4, 5)
    expect_equal(cs_resid, expected, tolerance = 1e-6)
})

test_that("residuals.dfr_dist computes Martingale residuals correctly", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))

    mart_resid <- residuals(exp_dist, df, type = "martingale")

    # Martingale = delta - H(t)
    expected <- c(1, 1, 1, 0, 0) - 0.5 * c(1, 2, 3, 4, 5)
    expect_equal(mart_resid, expected, tolerance = 1e-6)
})

test_that("Martingale residuals sum to approximately zero at MLE", {
    # At the MLE, martingale residuals should sum to approximately zero
    set.seed(42)
    true_lambda <- 0.5
    df <- data.frame(t = rexp(100, rate = true_lambda), delta = 1)

    # Fit exponential
    solver <- fit(dfr_exponential())
    result <- solver(df, par = c(1.0))

    # Create distribution with fitted parameters
    fitted_dist <- dfr_exponential(lambda = coef(result))

    mart_resid <- residuals(fitted_dist, df, type = "martingale")

    # Should sum to approximately zero
    expect_equal(sum(mart_resid), 0, tolerance = 0.5)
})

test_that("Cox-Snell residuals follow Exp(1) when model is correct", {
    # Generate data from known distribution
    set.seed(123)
    n <- 200
    true_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = sampler(true_dist)(n), delta = 1)

    # Fit the model
    solver <- fit(dfr_exponential())
    result <- solver(df, par = c(1.0))
    fitted_dist <- dfr_exponential(lambda = coef(result))

    # Get Cox-Snell residuals
    cs_resid <- residuals(fitted_dist, df, type = "cox-snell")

    # Should follow Exp(1) - check mean is approximately 1
    expect_equal(mean(cs_resid), 1, tolerance = 0.2)

    # Kolmogorov-Smirnov test against Exp(1)
    ks_test <- ks.test(cs_resid, "pexp", rate = 1)
    expect_true(ks_test$p.value > 0.05)  # Should not reject null
})

test_that("residuals requires parameters", {
    exp_dist <- dfr_exponential()  # No default parameters
    df <- data.frame(t = c(1, 2, 3), delta = c(1, 1, 0))

    expect_error(residuals(exp_dist, df), "Parameters required")

    # But works when par provided
    cs_resid <- residuals(exp_dist, df, par = c(0.5), type = "cox-snell")
    expect_equal(length(cs_resid), 3)
})

test_that("residuals handles missing delta column", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = c(1, 2, 3))  # No delta column

    # Cox-Snell works (doesn't need delta)
    cs_resid <- residuals(exp_dist, df, type = "cox-snell")
    expect_equal(length(cs_resid), 3)

    # Martingale assumes all events
    mart_resid <- residuals(exp_dist, df, type = "martingale")
    expected <- rep(1, 3) - 0.5 * c(1, 2, 3)
    expect_equal(mart_resid, expected, tolerance = 1e-6)
})

test_that("plot.dfr_dist creates survival plot without error", {
    exp_dist <- dfr_exponential(lambda = 0.5)

    # Should not error
    expect_silent({
        result <- plot(exp_dist, what = "survival", xlim = c(0, 5))
    })

    # Returns t and y values
    result <- plot(exp_dist, what = "survival", xlim = c(0, 5))
    expect_true("t" %in% names(result))
    expect_true("y" %in% names(result))
    expect_equal(length(result$t), 200)  # default n
})

test_that("plot.dfr_dist creates hazard plot", {
    weib_dist <- dfr_weibull(shape = 2, scale = 5)

    result <- plot(weib_dist, what = "hazard", xlim = c(0, 10))

    # Weibull k>1 has increasing hazard
    expect_true(result$y[length(result$y)] > result$y[1])
})

test_that("plot.dfr_dist creates cumhaz plot", {
    exp_dist <- dfr_exponential(lambda = 0.5)

    result <- plot(exp_dist, what = "cumhaz", xlim = c(0, 5))

    # Cumulative hazard is linear for exponential
    # H(t) = 0.5 * t, so H(5) should be approximately 2.5
    expect_equal(result$y[length(result$y)], 2.5, tolerance = 0.1)
})

test_that("plot.dfr_dist works with data overlay", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = rexp(50, rate = 0.5), delta = rep(1, 50))

    # Should not error with empirical overlay
    expect_silent({
        result <- plot(exp_dist, data = df, what = "survival", empirical = TRUE)
    })
})

test_that("plot.dfr_dist add parameter works", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    exp_dist2 <- dfr_exponential(lambda = 1.0)

    # First plot
    plot(exp_dist, what = "survival", xlim = c(0, 5), col = "blue")

    # Add second curve - should not error
    expect_silent({
        plot(exp_dist2, what = "survival", add = TRUE, col = "red")
    })
})

test_that("kaplan_meier computes correct estimates", {
    # Simple test case
    time <- c(1, 2, 3, 4, 5)
    delta <- c(1, 0, 1, 1, 0)  # Events at 1, 3, 4; censored at 2, 5

    km <- flexhaz:::kaplan_meier(time, delta)

    # Event times should be 1, 3, 4
    expect_equal(km$time, c(1, 3, 4))

    # At time 1: n=5, d=1, S = 4/5 = 0.8
    expect_equal(km$surv[1], 0.8)

    # At time 3: n=3 (1,2 gone), d=1, S = 0.8 * 2/3 ≈ 0.533
    expect_equal(km$surv[2], 0.8 * 2/3, tolerance = 1e-10)

    # At time 4: n=2, d=1, S = 0.533 * 1/2 ≈ 0.267
    expect_equal(km$surv[3], 0.8 * 2/3 * 1/2, tolerance = 1e-10)
})

test_that("kaplan_meier handles all censored data", {
    time <- c(1, 2, 3)
    delta <- c(0, 0, 0)  # All censored

    km <- flexhaz:::kaplan_meier(time, delta)

    # No event times
    expect_equal(length(km$time), 0)
})

test_that("qqplot_residuals runs without error", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(t = rexp(50, rate = 0.5), delta = rep(1, 50))

    # Should not error
    expect_silent({
        result <- qqplot_residuals(exp_dist, df)
    })

    # Returns Cox-Snell residuals
    expect_equal(length(result), 50)
})

test_that("qqplot_residuals handles censored data", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    df <- data.frame(
        t = c(rexp(30, rate = 0.5), rep(5, 20)),  # 20 censored at t=5
        delta = c(rep(1, 30), rep(0, 20))
    )

    # Should not error
    expect_silent({
        result <- qqplot_residuals(exp_dist, df)
    })

    # Returns all residuals (including censored)
    expect_equal(length(result), 50)
})

test_that("residuals works with Weibull distribution", {
    weib_dist <- dfr_weibull(shape = 2, scale = 5)
    df <- data.frame(t = c(1, 3, 5, 7, 10), delta = c(1, 1, 0, 1, 0))

    cs_resid <- residuals(weib_dist, df, type = "cox-snell")

    # H(t) = (t/sigma)^k = (t/5)^2
    expected <- (c(1, 3, 5, 7, 10) / 5)^2
    expect_equal(cs_resid, expected, tolerance = 1e-4)
})

test_that("residuals works with Gompertz distribution", {
    gomp_dist <- dfr_gompertz(a = 0.01, b = 0.1)
    df <- data.frame(t = c(1, 2, 3), delta = c(1, 1, 1))

    cs_resid <- residuals(gomp_dist, df, type = "cox-snell")

    # H(t) = (a/b) * (exp(b*t) - 1)
    a <- 0.01
    b <- 0.1
    expected <- (a / b) * (exp(b * c(1, 2, 3)) - 1)
    expect_equal(cs_resid, expected, tolerance = 1e-4)
})

test_that("residuals works with log-logistic distribution", {
    ll_dist <- dfr_loglogistic(alpha = 10, beta = 2)
    df <- data.frame(t = c(5, 10, 15), delta = c(1, 1, 1))

    cs_resid <- residuals(ll_dist, df, type = "cox-snell")

    # H(t) = log(1 + (t/alpha)^beta)
    alpha <- 10
    beta <- 2
    expected <- log(1 + (c(5, 10, 15) / alpha)^beta)
    expect_equal(cs_resid, expected, tolerance = 1e-4)
})

# =============================================================================
# Diagnostics with Custom Column Names
# =============================================================================

test_that("residuals errors when ob_col not found in data", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        ob_col = "survival_time",
        par = c(0.5)
    )
    df <- data.frame(t = c(1, 2, 3), delta = c(1, 1, 1))

    expect_error(residuals(dist, df), "Time column 'survival_time' not found")
})

test_that("residuals works with custom ob_col and delta_col", {
    dist <- dfr_dist(
        rate = function(t, par, ...) rep(par[1], length(t)),
        cum_haz_rate = function(t, par, ...) par[1] * t,
        ob_col = "time",
        delta_col = "status",
        par = c(0.5)
    )
    df <- data.frame(time = c(1, 2, 3), status = c(1, 0, 1))

    # Cox-Snell residuals: H(t) = 0.5 * t
    cs_resid <- residuals(dist, df, type = "cox-snell")
    expect_equal(cs_resid, c(0.5, 1.0, 1.5), tolerance = 1e-6)

    # Martingale residuals: delta - H(t)
    mart_resid <- residuals(dist, df, type = "martingale")
    expected <- c(1, 0, 1) - c(0.5, 1.0, 1.5)
    expect_equal(mart_resid, expected, tolerance = 1e-6)
})

test_that("plot with cumhaz and empirical overlay produces cumhaz curve", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    set.seed(42)
    df <- data.frame(
        t = rexp(50, rate = 0.5),
        delta = sample(c(0, 1), 50, replace = TRUE, prob = c(0.3, 0.7))
    )

    # Should not error, and should exercise the cumhaz empirical overlay path
    expect_silent({
        result <- plot(exp_dist, data = df, what = "cumhaz",
                       xlim = c(0, 5), empirical = TRUE)
    })

    expect_true("t" %in% names(result))
    expect_true("y" %in% names(result))
})

test_that("qqplot_residuals works when delta column is absent", {
    exp_dist <- dfr_exponential(lambda = 0.5)
    # Data without delta column - all observations treated as exact
    df <- data.frame(t = rexp(30, rate = 0.5))

    expect_silent({
        result <- qqplot_residuals(exp_dist, df)
    })

    # Should return all residuals (no filtering needed)
    expect_equal(length(result), 30)
})
