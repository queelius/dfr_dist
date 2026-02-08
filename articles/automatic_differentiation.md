# Custom Derivatives for Maximum Likelihood Estimation

## Introduction

Maximum likelihood estimation (MLE) for survival models requires
computing derivatives of the log-likelihood function. The **score
function** (gradient) is needed for optimization, and the **Hessian
matrix** (second derivatives) is needed for computing standard errors
and confidence intervals via the observed Fisher information.

The `dfr.dist` package uses a simple 2-tier fallback for each:

| Derivative | If provided              | Otherwise                                                              |
|------------|--------------------------|------------------------------------------------------------------------|
| Score      | `score_fn(df, par, ...)` | [`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html)       |
| Hessian    | `hess_fn(df, par, ...)`  | [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html) |

**You decide how to compute derivatives.** Hand-derive them, use an AD
library, or let the package fall back to numerical methods. The package
accepts functions; it doesn’t care how they were produced.

This vignette shows how to provide custom derivatives using the Weibull
distribution as a case study.

## Setup

``` r
library(dfr.dist)
```

## The Weibull Distribution

The Weibull distribution is widely used in survival analysis. Its hazard
function is:

$$h(t;k,\sigma) = \frac{k}{\sigma}\left( \frac{t}{\sigma} \right)^{k - 1}$$

The cumulative hazard is:

$$H(t;k,\sigma) = \left( \frac{t}{\sigma} \right)^{k}$$

For uncensored data, the log-likelihood is:

$$\ell(k,\sigma) = n\log k + (k - 1)\sum\limits_{i}\log t_{i} - nk\log\sigma - \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k}$$

## Simulating Data

``` r
set.seed(42)
true_k <- 2
true_sigma <- 3
n <- 100

u <- runif(n)
times <- true_sigma * (-log(u))^(1/true_k)
df <- data.frame(t = times, delta = rep(1, n))

cat("Sample size:", n, "\n")
#> Sample size: 100
cat("True parameters: k =", true_k, ", sigma =", true_sigma, "\n")
#> True parameters: k = 2 , sigma = 3
```

## Approach 1: Minimal — Just the Hazard Rate

The simplest approach provides only the hazard function. Everything else
is computed numerically:

``` r
weibull_minimal <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    }
)

ll <- loglik(weibull_minimal)
s <- score(weibull_minimal)
H <- hess_loglik(weibull_minimal)

test_par <- c(1.8, 2.8)

cat("Log-likelihood:", ll(df, par = test_par), "\n")
#> Log-likelihood: -179.2143
cat("Score (numerical):", s(df, par = test_par), "\n")
#> Score (numerical): -10.49419 8.878177
cat("Hessian (numerical):\n")
#> Hessian (numerical):
print(round(H(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4602  30.3417
#> [2,]  30.3417 -50.2047
```

**Pros:** No math required beyond the hazard function.

**Cons:** Slower (numerical integration for H(t), finite differences for
derivatives).

## Approach 2: Analytical Score + Analytical Hessian

For the Weibull, we can derive exact gradient and Hessian analytically:

**Score:**

$$\frac{\partial\ell}{\partial k} = \frac{n}{k} + \sum\limits_{i}\log t_{i} - n\log\sigma - \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k}\log\left( \frac{t_{i}}{\sigma} \right)$$

$$\frac{\partial\ell}{\partial\sigma} = \frac{k}{\sigma}\left\lbrack \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k} - n \right\rbrack$$

``` r
weibull_full <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)

        dk <- n_events / k + sum(delta * log_t_ratio) -
            sum(t_ratio^k * log_t_ratio)
        dsigma <- -n_events * k / sigma +
            (k / sigma) * sum(t_ratio^k)

        c(dk, dsigma)
    },
    hess_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)
        t_ratio_k <- t_ratio^k

        d2k <- -n_events / k^2 - sum(t_ratio_k * log_t_ratio^2)
        d2k_sigma <- -n_events / sigma +
            (1 / sigma) * sum(t_ratio_k) +
            (k / sigma) * sum(t_ratio_k * log_t_ratio)
        d2sigma <- n_events * k / sigma^2 -
            k * (k + 1) / sigma^2 * sum(t_ratio_k)

        matrix(c(d2k, d2k_sigma, d2k_sigma, d2sigma), nrow = 2)
    }
)

s_full <- score(weibull_full)
H_full <- hess_loglik(weibull_full)

cat("Score (analytical):", s_full(df, par = test_par), "\n")
#> Score (analytical): -10.4947 8.878147
cat("Hessian (analytical):\n")
#> Hessian (analytical):
print(round(H_full(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
```

## Approach 3: Provide Score Only, Let numDeriv Handle the Hessian

If deriving the Hessian is tedious but the score is straightforward,
provide just `score_fn`:

``` r
weibull_score_only <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)

        dk <- n_events / k + sum(delta * log_t_ratio) -
            sum(t_ratio^k * log_t_ratio)
        dsigma <- -n_events * k / sigma +
            (k / sigma) * sum(t_ratio^k)

        c(dk, dsigma)
    }
    # No hess_fn -> numDeriv::hessian fallback
)

H_mixed <- hess_loglik(weibull_score_only)
cat("Hessian (numerical fallback):\n")
#> Hessian (numerical fallback):
print(round(H_mixed(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
```

This gives you fast optimization (analytical gradient) with correct but
slower Hessian computation at the end.

## Using Pre-Built Constructors

The package provides helper constructors with analytical derivatives
already implemented:

``` r
# These include score_fn and hess_fn (where available)
weib <- dfr_weibull(shape = 2, scale = 3)
exp_d <- dfr_exponential(lambda = 0.5)

# Fit the model
solver <- fit(weib)
result <- solver(df, par = c(1.5, 2.5))

cat("Fitted parameters:\n")
#> Fitted parameters:
print(coef(result))
#> [1] 1.702321 2.963588
cat("\nTrue parameters: k =", true_k, ", sigma =", true_sigma, "\n")
#> 
#> True parameters: k = 2 , sigma = 3
cat("\n95% Confidence intervals:\n")
#> 
#> 95% Confidence intervals:
print(confint(result))
#>          2.5%    97.5%
#> [1,] 1.450448 1.954194
#> [2,] 2.603077 3.324099
```

## Computing Standard Errors

With the Hessian, standard errors come from the observed Fisher
information:

``` r
# Hessian at MLE
mle_par <- coef(result)
hess_mle <- H_full(df, par = mle_par)

# Observed Fisher information = -Hessian
obs_info <- -hess_mle

# Standard errors = sqrt(diag(inverse Fisher info))
se <- sqrt(diag(solve(obs_info)))

cat("Standard errors:\n")
#> Standard errors:
cat("  SE(k) =", round(se[1], 4), "\n")
#>   SE(k) = 0.1285
cat("  SE(sigma) =", round(se[2], 4), "\n")
#>   SE(sigma) = 0.1839
```

## Distribution Reference

| Distribution | Constructor                                                                             | Score      | Hessian    |
|--------------|-----------------------------------------------------------------------------------------|------------|------------|
| Exponential  | [`dfr_exponential()`](https://queelius.github.io/dfr.dist/reference/dfr_exponential.md) | Analytical | Analytical |
| Weibull      | [`dfr_weibull()`](https://queelius.github.io/dfr.dist/reference/dfr_weibull.md)         | Analytical | Analytical |
| Gompertz     | [`dfr_gompertz()`](https://queelius.github.io/dfr.dist/reference/dfr_gompertz.md)       | Analytical | numDeriv   |
| Log-logistic | [`dfr_loglogistic()`](https://queelius.github.io/dfr.dist/reference/dfr_loglogistic.md) | Analytical | numDeriv   |

For Gompertz and log-logistic, you can supply your own `hess_fn` if
needed.

## Next Steps

- **[`vignette("getting_started")`](https://queelius.github.io/dfr.dist/articles/getting_started.md)** -
  Quick 5-minute introduction
- **[`vignette("reliability_engineering")`](https://queelius.github.io/dfr.dist/articles/reliability_engineering.md)** -
  Real-world applications
- **[`vignette("failure_rate")`](https://queelius.github.io/dfr.dist/articles/failure_rate.md)** -
  Mathematical foundations
- **[`vignette("custom_distributions")`](https://queelius.github.io/dfr.dist/articles/custom_distributions.md)** -
  The three-level optimization paradigm
