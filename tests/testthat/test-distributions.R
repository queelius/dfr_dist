# Tests for helper distribution constructors
# Tests verify analytical formulas against numerical computation and known values

test_that("dfr_exponential creates valid distribution", {
    # With parameter
    exp_dist <- dfr_exponential(lambda = 0.5)
    expect_true(is_dfr_dist(exp_dist))

    # Hazard is constant
    h <- hazard(exp_dist)
    expect_equal(h(1), 0.5)
    expect_equal(h(10), 0.5)
    expect_equal(h(c(1, 5, 10)), c(0.5, 0.5, 0.5))

    # Cumulative hazard is linear
    H <- cum_haz(exp_dist)
    expect_equal(H(2), 1.0)
    expect_equal(H(4), 2.0)

    # Survival matches analytical
    S <- surv(exp_dist)
    expect_equal(S(2), exp(-1.0), tolerance = 1e-6)

    # Without parameter (for fitting)
    exp_null <- dfr_exponential()
    expect_true(is_dfr_dist(exp_null))
    h_null <- hazard(exp_null)
    expect_equal(h_null(1, par = c(2)), 2)
})

test_that("dfr_exponential score function is correct", {
    exp_dist <- dfr_exponential()
    s <- score(exp_dist)

    # Create test data
    df <- data.frame(t = c(1, 2, 3, 4, 5), delta = c(1, 1, 1, 0, 0))

    # At the MLE, score should be approximately zero
    # MLE for exponential: lambda = n_events / sum(t)
    n_events <- sum(df$delta == 1)
    total_time <- sum(df$t)
    lambda_mle <- n_events / total_time

    score_at_mle <- s(df, par = c(lambda_mle))
    expect_equal(score_at_mle, 0, tolerance = 1e-10)
})

test_that("dfr_exponential MLE matches analytical solution", {
    set.seed(42)
    true_lambda <- 0.5
    n <- 100
    times <- rexp(n, rate = true_lambda)
    df <- data.frame(t = times, delta = rep(1, n))

    # Analytical MLE
    analytical_mle <- n / sum(times)

    # Fitted MLE
    solver <- fit(dfr_exponential())
    result <- solver(df, par = c(1.0))

    expect_equal(coef(result), analytical_mle, tolerance = 1e-6)
})

test_that("dfr_weibull creates valid distribution", {
    # With parameters
    weib_dist <- dfr_weibull(shape = 2, scale = 3)
    expect_true(is_dfr_dist(weib_dist))

    # Hazard formula: h(t) = (k/sigma) * (t/sigma)^(k-1)
    h <- hazard(weib_dist)
    t <- 2
    k <- 2
    sigma <- 3
    expected_h <- (k / sigma) * (t / sigma)^(k - 1)
    expect_equal(h(t), expected_h, tolerance = 1e-10)

    # Cumulative hazard: H(t) = (t/sigma)^k
    H <- cum_haz(weib_dist)
    expected_H <- (t / sigma)^k
    expect_equal(H(t), expected_H, tolerance = 1e-6)

    # Survival
    S <- surv(weib_dist)
    expect_equal(S(t), exp(-expected_H), tolerance = 1e-6)

    # Shape = 1 reduces to exponential
    weib_exp <- dfr_weibull(shape = 1, scale = 2)  # lambda = 1/2
    h_weib <- hazard(weib_exp)
    expect_equal(h_weib(5), 0.5, tolerance = 1e-10)
})

test_that("dfr_weibull score function is correct", {
    weib_dist <- dfr_weibull()
    s <- score(weib_dist)

    # Generate Weibull data
    set.seed(123)
    true_k <- 2
    true_sigma <- 3
    n <- 50
    u <- runif(n)
    times <- true_sigma * (-log(u))^(1 / true_k)
    df <- data.frame(t = times, delta = rep(1, n))

    # Fit the model
    solver <- fit(weib_dist)
    result <- solver(df, par = c(1.5, 2.5))

    # Score at MLE should be near zero
    score_at_mle <- s(df, par = coef(result))
    expect_equal(score_at_mle, c(0, 0), tolerance = 1e-3)
})

test_that("dfr_weibull MLE recovers true parameters", {
    set.seed(456)
    true_k <- 2.5
    true_sigma <- 10
    n <- 200

    # Inverse CDF sampling for Weibull
    u <- runif(n)
    times <- true_sigma * (-log(u))^(1 / true_k)
    df <- data.frame(t = times, delta = rep(1, n))

    solver <- fit(dfr_weibull())
    result <- solver(df, par = c(2, 8))

    # Should recover parameters reasonably well
    expect_equal(coef(result)[1], true_k, tolerance = 0.3)
    expect_equal(coef(result)[2], true_sigma, tolerance = 1.0)
})

test_that("dfr_gompertz creates valid distribution", {
    # With parameters
    gomp_dist <- dfr_gompertz(a = 0.01, b = 0.1)
    expect_true(is_dfr_dist(gomp_dist))

    # Hazard formula: h(t) = a * exp(b * t)
    h <- hazard(gomp_dist)
    t <- 5
    a <- 0.01
    b <- 0.1
    expected_h <- a * exp(b * t)
    expect_equal(h(t), expected_h, tolerance = 1e-10)

    # Cumulative hazard: H(t) = (a/b) * (exp(b*t) - 1)
    H <- cum_haz(gomp_dist)
    expected_H <- (a / b) * (exp(b * t) - 1)
    expect_equal(H(t), expected_H, tolerance = 1e-6)

    # Survival
    S <- surv(gomp_dist)
    expect_equal(S(t), exp(-expected_H), tolerance = 1e-6)

    # Initial hazard equals a
    expect_equal(h(0), a, tolerance = 1e-10)
})

test_that("dfr_gompertz hazard doubles correctly", {
    # Set b = log(2)/T so hazard doubles every T units
    T_double <- 10
    b <- log(2) / T_double
    gomp_dist <- dfr_gompertz(a = 0.01, b = b)

    h <- hazard(gomp_dist)

    # Hazard at t=0
    h0 <- h(0)

    # Hazard should double at t = T_double
    expect_equal(h(T_double), 2 * h0, tolerance = 1e-10)

    # And quadruple at 2*T_double
    expect_equal(h(2 * T_double), 4 * h0, tolerance = 1e-10)
})

test_that("dfr_loglogistic creates valid distribution", {
    # With parameters
    ll_dist <- dfr_loglogistic(alpha = 10, beta = 2)
    expect_true(is_dfr_dist(ll_dist))

    # Hazard formula: h(t) = (beta/alpha) * (t/alpha)^(beta-1) / (1 + (t/alpha)^beta)
    h <- hazard(ll_dist)
    t <- 5
    alpha <- 10
    beta <- 2
    t_ratio <- t / alpha
    expected_h <- (beta / alpha) * t_ratio^(beta - 1) / (1 + t_ratio^beta)
    expect_equal(h(t), expected_h, tolerance = 1e-10)

    # Cumulative hazard: H(t) = log(1 + (t/alpha)^beta)
    H <- cum_haz(ll_dist)
    expected_H <- log(1 + t_ratio^beta)
    expect_equal(H(t), expected_H, tolerance = 1e-6)

    # Survival: S(t) = 1 / (1 + (t/alpha)^beta)
    S <- surv(ll_dist)
    expected_S <- 1 / (1 + t_ratio^beta)
    expect_equal(S(t), expected_S, tolerance = 1e-6)
})

test_that("dfr_loglogistic median equals alpha when beta > 1", {
    alpha <- 100
    beta <- 3

    ll_dist <- dfr_loglogistic(alpha = alpha, beta = beta)
    Q <- inv_cdf(ll_dist)

    # Median should equal alpha
    median <- Q(0.5)
    expect_equal(median, alpha, tolerance = 1e-3)
})

test_that("dfr_loglogistic hazard is non-monotonic for beta > 1", {
    ll_dist <- dfr_loglogistic(alpha = 10, beta = 3)
    h <- hazard(ll_dist)

    # Hazard should increase initially, then decrease
    t_early <- 2
    t_mid <- 10
    t_late <- 50

    h_early <- h(t_early)
    h_mid <- h(t_mid)
    h_late <- h(t_late)

    # Early should be less than mid (increasing phase)
    expect_true(h_early < h_mid)

    # Late should be less than mid (decreasing phase)
    expect_true(h_late < h_mid)
})

test_that("all distributions handle censored data", {
    # Test that fitting works with mixed exact/censored data
    df <- data.frame(
        t = c(1, 2, 3, 4, 5, 6, 7, 8),
        delta = c(1, 1, 1, 0, 0, 1, 1, 0)
    )

    # Exponential
    solver_exp <- fit(dfr_exponential())
    result_exp <- solver_exp(df, par = c(0.5))
    expect_true(result_exp$converged)

    # Weibull
    solver_weib <- fit(dfr_weibull())
    result_weib <- solver_weib(df, par = c(1.5, 5))
    expect_true(result_weib$converged)

    # Gompertz
    solver_gomp <- fit(dfr_gompertz())
    result_gomp <- solver_gomp(df, par = c(0.1, 0.1))
    expect_true(result_gomp$converged)

    # Log-logistic
    solver_ll <- fit(dfr_loglogistic())
    result_ll <- solver_ll(df, par = c(5, 2))
    expect_true(result_ll$converged)
})

test_that("distributions work with parameter override", {
    # All distribution methods should accept par argument

    exp_dist <- dfr_exponential(lambda = 1)
    h <- hazard(exp_dist)
    expect_equal(h(1), 1)  # default
    expect_equal(h(1, par = c(2)), 2)  # override

    weib_dist <- dfr_weibull(shape = 2, scale = 1)
    S <- surv(weib_dist)
    S_default <- S(1)
    S_override <- S(1, par = c(3, 2))
    expect_true(S_default != S_override)

    gomp_dist <- dfr_gompertz(a = 0.1, b = 0.1)
    H <- cum_haz(gomp_dist)
    H_default <- H(5)
    H_override <- H(5, par = c(0.2, 0.2))
    expect_true(H_default != H_override)

    ll_dist <- dfr_loglogistic(alpha = 10, beta = 2)
    cdf_fn <- cdf(ll_dist)
    F_default <- cdf_fn(10)
    F_override <- cdf_fn(10, par = c(20, 3))
    expect_true(F_default != F_override)
})

test_that("analytical cumulative hazard matches numerical integration",
          {
              # Test that our analytical cum_haz_rate matches numerical integration

              # Exponential
              exp_dist <- dfr_exponential(lambda = 0.5)
              H_analytical <- cum_haz(exp_dist)
              h <- hazard(exp_dist)
              H_numerical <- function(t) {
                  integrate(h, 0, t)$value
              }
              expect_equal(H_analytical(5), H_numerical(5), tolerance = 1e-3)

              # Weibull
              weib_dist <- dfr_weibull(shape = 2, scale = 3)
              H_analytical <- cum_haz(weib_dist)
              h <- hazard(weib_dist)
              H_numerical <- function(t) {
                  integrate(h, 0, t)$value
              }
              expect_equal(H_analytical(5), H_numerical(5), tolerance = 1e-3)

              # Gompertz
              gomp_dist <- dfr_gompertz(a = 0.01, b = 0.1)
              H_analytical <- cum_haz(gomp_dist)
              h <- hazard(gomp_dist)
              H_numerical <- function(t) {
                  integrate(h, 0, t)$value
              }
              expect_equal(H_analytical(5), H_numerical(5), tolerance = 1e-3)

              # Log-logistic
              ll_dist <- dfr_loglogistic(alpha = 10, beta = 2)
              H_analytical <- cum_haz(ll_dist)
              h <- hazard(ll_dist)
              H_numerical <- function(t) {
                  integrate(h, 0, t)$value
              }
              expect_equal(H_analytical(5), H_numerical(5), tolerance = 1e-3)
          })

test_that("analytical score matches numerical gradient", {
    skip_if_not_installed("numDeriv")

    # Test that analytical score_fn matches numerical gradient

    df <- data.frame(t = c(1, 2, 3, 4, 5), delta = c(1, 1, 0, 1, 0))

    # Exponential
    exp_dist <- dfr_exponential()
    ll_exp <- loglik(exp_dist)
    s_exp <- score(exp_dist)
    par_exp <- c(0.3)

    grad_numerical <- numDeriv::grad(function(p) ll_exp(df, p), par_exp)
    grad_analytical <- s_exp(df, par_exp)
    expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)

    # Weibull
    weib_dist <- dfr_weibull()
    ll_weib <- loglik(weib_dist)
    s_weib <- score(weib_dist)
    par_weib <- c(1.5, 3)

    grad_numerical <- numDeriv::grad(function(p) ll_weib(df, p), par_weib)
    grad_analytical <- s_weib(df, par_weib)
    expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)

    # Gompertz
    gomp_dist <- dfr_gompertz()
    ll_gomp <- loglik(gomp_dist)
    s_gomp <- score(gomp_dist)
    par_gomp <- c(0.1, 0.2)

    grad_numerical <- numDeriv::grad(function(p) ll_gomp(df, p), par_gomp)
    grad_analytical <- s_gomp(df, par_gomp)
    expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)

    # Log-logistic
    ll_dist <- dfr_loglogistic()
    ll_ll <- loglik(ll_dist)
    s_ll <- score(ll_dist)
    par_ll <- c(3, 2)

    grad_numerical <- numDeriv::grad(function(p) ll_ll(df, p), par_ll)
    grad_analytical <- s_ll(df, par_ll)
    expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})
