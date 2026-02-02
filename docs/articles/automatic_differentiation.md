# Automatic Differentiation for Maximum Likelihood Estimation

## Introduction

Maximum likelihood estimation (MLE) for survival models requires
computing derivatives of the log-likelihood function. The **score
function** (gradient) is needed for optimization, and the **Hessian
matrix** (second derivatives) is needed for computing standard errors
and confidence intervals via the observed Fisher information.

The `dfr.dist` package supports three approaches for computing these
derivatives:

1.  **Pure numerical differentiation** - finite differences via
    `numDeriv`
2.  **Full automatic differentiation** - AD computes both score and
    Hessian
3.  **Hybrid approach** - analytical score + AD Jacobian for Hessian

This vignette explores these approaches using the Weibull distribution
as a case study, demonstrating that the **hybrid approach** offers the
best combination of speed and precision.

**Note:** The package provides helper constructors
([`dfr_exponential()`](https://queelius.github.io/dfr_dist/reference/dfr_exponential.md),
[`dfr_weibull()`](https://queelius.github.io/dfr_dist/reference/dfr_weibull.md),
[`dfr_gompertz()`](https://queelius.github.io/dfr_dist/reference/dfr_gompertz.md),
[`dfr_loglogistic()`](https://queelius.github.io/dfr_dist/reference/dfr_loglogistic.md))
that include optimized implementations with analytical cumulative
hazards and score functions. This vignette shows how to build these from
scratch to understand the underlying principles.

## Setup

``` r
library(dfr.dist)

# Check if femtograd is available for AD
# Note: we use requireNamespace() rather than library() to avoid
# femtograd::fit() masking the generics::fit() S3 generic
has_ad <- requireNamespace("femtograd", quietly = TRUE)
if (has_ad) {
  cat("femtograd available - AD examples will run\n")
} else {
  cat("femtograd not installed - AD examples will be skipped\n
       Install with: devtools::install_github('queelius/femtograd')\n")
}
#> femtograd available - AD examples will run
```

## How femtograd Works

The `femtograd` package implements two AD modes that `dfr.dist` uses for
different computations:

- **Reverse-mode AD** (for gradients): femtograd builds a computational
  graph of `value` objects during the forward pass, then calls
  `backward()` in reverse topological order to compute all partial
  derivatives in a single backward pass. This is what
  [`ad_gradient()`](https://queelius.github.io/dfr_dist/reference/ad_gradient.md)
  uses — efficient for scalar outputs (like log-likelihood) with many
  parameters.

- **Forward-over-reverse** (for Hessians): To compute second
  derivatives, femtograd wraps its `value` graph nodes in `dual_num`
  objects, combining forward-mode dual numbers with the reverse-mode
  graph. Our
  [`ad_jacobian()`](https://queelius.github.io/dfr_dist/reference/ad_jacobian.md)
  uses femtograd’s dual numbers directly in a forward-mode pass to
  compute the Jacobian of the score function column-by-column — one pass
  per parameter.

The **hybrid approach** (analytical score + AD Jacobian) leverages this:
you provide an analytical gradient, and
[`ad_jacobian()`](https://queelius.github.io/dfr_dist/reference/ad_jacobian.md)
differentiates it using forward-mode AD to obtain the Hessian.

## The Weibull Distribution

The Weibull distribution is widely used in survival analysis. Its hazard
function is:

$$h(t;k,\sigma) = \frac{k}{\sigma}\left( \frac{t}{\sigma} \right)^{k - 1}$$

The cumulative hazard is:

$$H(t;k,\sigma) = \left( \frac{t}{\sigma} \right)^{k}$$

For uncensored data, the log-likelihood is:

$$\ell(k,\sigma) = n\log k + (k - 1)\sum\limits_{i}\log t_{i} - nk\log\sigma - \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k}$$

## Simulating Weibull Data

``` r
set.seed(42)
true_k <- 2
true_sigma <- 3
n <- 100

# Inverse CDF sampling for Weibull
u <- runif(n)
times <- true_sigma * (-log(u))^(1/true_k)
df <- data.frame(t = times, delta = rep(1, n))

cat("Sample size:", n, "\n")
#> Sample size: 100
cat("True parameters: k =", true_k, ", sigma =", true_sigma, "\n")
#> True parameters: k = 2 , sigma = 3
cat("Sample mean:", round(mean(times), 3), "\n")
#> Sample mean: 2.634
```

## Approach 1: Pure Numerical Differentiation

The simplest approach uses finite differences to approximate
derivatives. This requires only the hazard rate function:

``` r
weibull_numerical <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    }
)

# Get likelihood functions
ll <- loglik(weibull_numerical)
score_num <- score(weibull_numerical)
hess_num <- hess_loglik(weibull_numerical)

# Evaluate at test point
test_par <- c(1.8, 2.8)

cat("Log-likelihood:", ll(df, par = test_par), "\n\n")
#> Log-likelihood: -179.2143

cat("Score (numerical gradient):\n")
#> Score (numerical gradient):
print(score_num(df, par = test_par))
#> [1] -10.494194   8.878177

cat("\nHessian (numerical):\n")
#> 
#> Hessian (numerical):
print(round(hess_num(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4602  30.3417
#> [2,]  30.3417 -50.2047
```

**Pros:** Simple, no analytical derivatives needed.

**Cons:** Slow (many function evaluations), potential numerical
precision issues.

## Approach 2: Hybrid - Analytical Score + AD Jacobian

The score function for Weibull can be derived analytically:

$$\frac{\partial\ell}{\partial k} = \frac{n}{k} + \sum\limits_{i}\log t_{i} - n\log\sigma - \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k}\log\left( \frac{t_{i}}{\sigma} \right)$$

$$\frac{\partial\ell}{\partial\sigma} = \frac{k}{\sigma}\left\lbrack \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k} - n \right\rbrack$$

We provide this analytical score, and let AD compute the Jacobian (which
gives us the Hessian):

``` r
weibull_hybrid <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[[1]]
        sigma <- par[[2]]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[[1]]
        sigma <- par[[2]]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k <- par[[1]]
        sigma <- par[[2]]
        t <- df$t
        n <- length(t)

        # Vectorized computation
        t_over_sigma <- t / sigma
        t_over_sigma_k <- t_over_sigma^k

        score_k <- n/k + sum(log(t)) - n*log(sigma) -
                   sum(t_over_sigma_k * log(t_over_sigma))
        score_sigma <- k/sigma * (sum(t_over_sigma_k) - n)

        c(score_k, score_sigma)
    }
)

score_hybrid <- score(weibull_hybrid)
hess_hybrid <- hess_loglik(weibull_hybrid)

cat("Score (analytical):\n")
#> Score (analytical):
print(score_hybrid(df, par = test_par))
#> [1] -10.494704   8.878147

cat("\nHessian (AD Jacobian of analytical score):\n")
#> 
#> Hessian (AD Jacobian of analytical score):
print(round(hess_hybrid(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
```

**Key point:** The `score_fn` uses `par[[1]]` indexing (double brackets)
for AD compatibility. This allows femtograd’s dual number objects to
flow through the Jacobian computation (which uses forward-mode AD to
differentiate the score).

## Approach 3: Full AD from Log-Likelihood

With AD, we can also compute the score directly from the log-likelihood,
then get the Hessian as the Jacobian of that score:

``` r
weibull_full_ad <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[[1]]
        sigma <- par[[2]]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[[1]]
        sigma <- par[[2]]
        (t / sigma)^k
    }
    # No score_fn - AD will compute it from log-likelihood
)

score_full_ad <- score(weibull_full_ad)
hess_full_ad <- hess_loglik(weibull_full_ad)

cat("Score (AD gradient of log-likelihood):\n")
#> Score (AD gradient of log-likelihood):
print(score_full_ad(df, par = test_par))
#> [1] -10.494704   8.878147

cat("\nHessian (AD Jacobian of AD score):\n")
#> 
#> Hessian (AD Jacobian of AD score):
print(round(hess_full_ad(df, par = test_par), 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
```

## Comparing Precision

All three approaches should give similar results. Let’s compare:

``` r
h_num <- hess_num(df, par = test_par)
h_hybrid <- hess_hybrid(df, par = test_par)
h_full <- hess_full_ad(df, par = test_par)

cat("Max absolute differences:\n")
#> Max absolute differences:
cat("  Hybrid vs Numerical:", max(abs(h_hybrid - h_num)), "\n")
#>   Hybrid vs Numerical: 0.0006991277
cat("  Full AD vs Numerical:", max(abs(h_full - h_num)), "\n")
#>   Full AD vs Numerical: 0.0006991277
cat("  Hybrid vs Full AD:", max(abs(h_hybrid - h_full)), "\n")
#>   Hybrid vs Full AD: 5.558789e-10
```

All methods agree to high precision (differences \< 0.01).

## Speed Benchmark

The real difference is in computational speed:

``` r
# Warmup
invisible(hess_num(df, par = test_par))
invisible(hess_hybrid(df, par = test_par))
invisible(hess_full_ad(df, par = test_par))

n_iter <- 20

# Numerical
t0 <- Sys.time()
for (i in 1:n_iter) invisible(hess_num(df, par = test_par))
time_num <- as.numeric(Sys.time() - t0) / n_iter * 1000

# Hybrid
t0 <- Sys.time()
for (i in 1:n_iter) invisible(hess_hybrid(df, par = test_par))
time_hybrid <- as.numeric(Sys.time() - t0) / n_iter * 1000

# Full AD
t0 <- Sys.time()
for (i in 1:n_iter) invisible(hess_full_ad(df, par = test_par))
time_full <- as.numeric(Sys.time() - t0) / n_iter * 1000

cat("Average time per Hessian computation:\n")
#> Average time per Hessian computation:
cat(sprintf("  Numerical:      %7.2f ms\n", time_num))
#>   Numerical:       238.21 ms
cat(sprintf("  Hybrid (recommended): %7.2f ms\n", time_hybrid))
#>   Hybrid (recommended):    3.89 ms
cat(sprintf("  Full AD:        %7.2f ms\n", time_full))
#>   Full AD:          46.82 ms
cat(sprintf("\nSpeedup (Numerical / Hybrid): %.1fx\n", time_num / time_hybrid))
#> 
#> Speedup (Numerical / Hybrid): 61.3x
```

## Scaling with Sample Size

The advantage of the hybrid approach grows with sample size:

``` r
cat("Scaling with sample size:\n\n")
#> Scaling with sample size:
cat(sprintf("%8s %10s %10s %8s\n", "n", "Hybrid", "Numerical", "Speedup"))
#>        n     Hybrid  Numerical  Speedup
cat(sprintf("%8s %10s %10s %8s\n", "---", "------", "---------", "-------"))
#>      ---     ------  ---------  -------

for (sample_n in c(50, 100, 500, 1000)) {
    set.seed(42)
    times_bench <- true_sigma * (-log(runif(sample_n)))^(1/true_k)
    df_bench <- data.frame(t = times_bench, delta = rep(1, sample_n))

    # Time hybrid
    t0 <- Sys.time()
    for (i in 1:10) invisible(hess_hybrid(df_bench, par = test_par))
    t_hybrid <- as.numeric(Sys.time() - t0) / 10 * 1000

    # Time numerical
    t0 <- Sys.time()
    for (i in 1:5) invisible(hess_num(df_bench, par = test_par))
    t_num <- as.numeric(Sys.time() - t0) / 5 * 1000

    cat(sprintf("%8d %9.1f ms %9.1f ms %7.0fx\n",
                sample_n, t_hybrid, t_num, t_num/t_hybrid))
}
#>       50       3.3 ms     117.6 ms      36x
#>      100       3.4 ms     239.8 ms      70x
#>      500       2.7 ms    1101.2 ms     411x
#>     1000       3.2 ms    2240.6 ms     698x
```

The hybrid approach stays nearly constant (~2ms) while numerical scales
linearly with n.

## Why is the Hybrid Approach Fastest?

1.  **Numerical differentiation** requires many log-likelihood
    evaluations, each involving numerical integration for the cumulative
    hazard.

2.  **Full AD** uses femtograd’s reverse-mode mechanism, building
    computation graphs for every observation, which has significant
    overhead.

3.  **Hybrid approach** uses an analytical score that directly computes
    the gradient without numerical integration. The AD Jacobian then
    requires only 2 forward passes (one per parameter) on this efficient
    score function.

## Writing AD-Compatible Score Functions

For the hybrid approach, your `score_fn` must be compatible with
femtograd’s dual number objects (used in the forward-mode Jacobian
pass):

### Do’s:

- Use `par[[1]]` (double brackets) to access parameters
- Use vectorized operations: `t / sigma`, `t^k`
- Use [`sum()`](https://rdrr.io/r/base/sum.html) to reduce vectors to
  scalars
- Use standard math functions:
  [`log()`](https://rdrr.io/r/base/Log.html),
  [`exp()`](https://rdrr.io/r/base/Log.html), `^`

### Don’ts:

- Don’t use `par[1]` (single bracket) - it won’t extract dual numbers
  correctly
- Avoid control flow that depends on parameter values
- Don’t use functions that don’t have AD implementations

### Example Pattern:

``` r
score_fn = function(df, par, ...) {
    # Extract parameters with [[ ]]
    theta1 <- par[[1]]
    theta2 <- par[[2]]

    # Get data (numeric, not dual)
    t <- df$t
    n <- length(t)

    # Vectorized operations work with duals
    z <- t / theta2

    # Sum reduces to scalar dual
    score1 <- n / theta1 - sum(z^theta1)
    score2 <- -n * theta1 / theta2 + ...

    # Return as vector
    c(score1, score2)
}
```

## Computing Standard Errors and Confidence Intervals

With the Hessian, we can compute standard errors via the observed Fisher
information:

``` r
# Find MLE
mle_result <- optim(c(1.5, 2.5),
                    function(p) -ll(df, par = p),
                    method = "L-BFGS-B",
                    lower = c(0.1, 0.1))
mle_par <- mle_result$par

cat("MLE estimates:\n")
#> MLE estimates:
cat("  k =", round(mle_par[1], 4), "(true:", true_k, ")\n")
#>   k = 1.7023 (true: 2 )
cat("  sigma =", round(mle_par[2], 4), "(true:", true_sigma, ")\n\n")
#>   sigma = 2.9636 (true: 3 )

# Hessian at MLE
hess_mle <- hess_hybrid(df, par = mle_par)

# Observed Fisher information = -Hessian
obs_info <- -hess_mle

# Standard errors = sqrt(diag(inverse Fisher info))
se <- sqrt(diag(solve(obs_info)))

cat("Standard errors:\n")
#> Standard errors:
cat("  SE(k) =", round(se[1], 4), "\n")
#>   SE(k) = 0.1285
cat("  SE(sigma) =", round(se[2], 4), "\n\n")
#>   SE(sigma) = 0.1839

# 95% confidence intervals
cat("95% Confidence Intervals:\n")
#> 95% Confidence Intervals:
cat("  k: [", round(mle_par[1] - 1.96*se[1], 3), ",",
    round(mle_par[1] + 1.96*se[1], 3), "]\n")
#>   k: [ 1.45 , 1.954 ]
cat("  sigma: [", round(mle_par[2] - 1.96*se[2], 3), ",",
    round(mle_par[2] + 1.96*se[2], 3), "]\n")
#>   sigma: [ 2.603 , 3.324 ]
```

## Summary

| Approach   | Speed    | Precision     | Effort       |
|------------|----------|---------------|--------------|
| Numerical  | Slow     | Good          | Minimal      |
| Full AD    | Slow     | Excellent     | Moderate     |
| **Hybrid** | **Fast** | **Excellent** | **Moderate** |

**Recommendation:** For production use, implement analytical score
functions and let AD compute the Hessian via Jacobian. This gives you:

- 70-2000x speedup over numerical (depending on sample size)
- Exact derivatives (no finite difference approximation)
- Automatic Hessian computation (no manual second derivatives)

The initial effort to derive the score function pays off substantially
in computational efficiency.

## Additional Distribution Examples

Below are several common survival distributions with their analytical
score functions, ready for use with the hybrid AD approach.

### Exponential Distribution

The simplest case - constant hazard rate:

$$h(t;\lambda) = \lambda,\quad H(t;\lambda) = \lambda t$$

``` r
exponential_dist <- dfr_dist(
    rate = function(t, par, ...) rep(par[[1]], length(t)),
    cum_haz_rate = function(t, par, ...) par[[1]] * t,
    score_fn = function(df, par, ...) {
        lambda <- par[[1]]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, length(t))
        # Score: sum(delta)/lambda - sum(t)
        c(sum(delta) / lambda - sum(t))
    }
)

# Test
set.seed(123)
exp_times <- rexp(50, rate = 0.5)
exp_df <- data.frame(t = exp_times, delta = rep(1, 50))

ll_exp <- loglik(exponential_dist)
score_exp <- score(exponential_dist)

cat("Exponential at lambda=0.5:\n")
#> Exponential at lambda=0.5:
cat("  Log-likelihood:", round(ll_exp(exp_df, par = c(0.5)), 3), "\n")
#>   Log-likelihood: -91.176
cat("  Score:", round(score_exp(exp_df, par = c(0.5)), 6), "(should be ~0 at MLE)\n")
#>   Score: -13.03708 (should be ~0 at MLE)
```

### Gompertz Distribution

Used in actuarial science for modeling mortality. The hazard rate
increases exponentially:

$$h(t;a,b) = ae^{bt},\quad H(t;a,b) = \frac{a}{b}\left( e^{bt} - 1 \right)$$

Score functions:
$$\frac{\partial\ell}{\partial a} = \frac{n}{a} - \frac{1}{b}\sum\limits_{i}\left( e^{bt_{i}} - 1 \right)$$$$\frac{\partial\ell}{\partial b} = \sum\limits_{i}t_{i} - \frac{a}{b}\sum\limits_{i}\left( \frac{t_{i}e^{bt_{i}}}{1} - \frac{e^{bt_{i}} - 1}{b} \right)$$

``` r
gompertz_dist <- dfr_dist(
    rate = function(t, par, ...) {
        a <- par[[1]]
        b <- par[[2]]
        a * exp(b * t)
    },
    cum_haz_rate = function(t, par, ...) {
        a <- par[[1]]
        b <- par[[2]]
        (a / b) * (exp(b * t) - 1)
    },
    score_fn = function(df, par, ...) {
        a <- par[[1]]
        b <- par[[2]]
        t <- df$t
        n <- length(t)
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, n)

        exp_bt <- exp(b * t)

        # Score w.r.t. a: sum(delta)/a - (1/b) * sum(exp_bt - 1)
        score_a <- sum(delta) / a - sum(exp_bt - 1) / b

        # Score w.r.t. b: sum(delta * t) - (a/b) * sum(t * exp_bt) + (a/b^2) * sum(exp_bt - 1)
        score_b <- sum(delta * t) - (a / b) * sum(t * exp_bt) +
                   (a / b^2) * sum(exp_bt - 1)

        c(score_a, score_b)
    }
)

# Test with simulated Gompertz data
# Using inverse CDF: T = (1/b) * log(1 - (b/a) * log(U))
set.seed(456)
a_true <- 0.01
b_true <- 0.5
u <- runif(80)
# Inverse CDF for Gompertz (when b > 0)
gompertz_times <- (1/b_true) * log(1 - (b_true/a_true) * log(1-u))
gompertz_times <- gompertz_times[is.finite(gompertz_times) & gompertz_times > 0]
gompertz_df <- data.frame(t = gompertz_times, delta = rep(1, length(gompertz_times)))

score_gomp <- score(gompertz_dist)
hess_gomp <- hess_loglik(gompertz_dist)

cat("Gompertz distribution (a=0.01, b=0.5):\n")
#> Gompertz distribution (a=0.01, b=0.5):
cat("  Score at true params:", round(score_gomp(gompertz_df, par = c(a_true, b_true)), 4), "\n")
#>   Score at true params: -276.3175 -13.9818
cat("  Hessian:\n")
#>   Hessian:
print(round(hess_gomp(gompertz_df, par = c(a_true, b_true)), 4))
#>            [,1]       [,2]
#> [1,] -800000.00 -56937.280
#> [2,]  -56937.28  -4327.902
```

### Log-Logistic Distribution

Common in survival analysis, especially for accelerated failure time
models:

$$h(t;\alpha,\beta) = \frac{(\beta/\alpha)(t/\alpha)^{\beta - 1}}{1 + (t/\alpha)^{\beta}},\quad H(t;\alpha,\beta) = \log\left( 1 + (t/\alpha)^{\beta} \right)$$

``` r
# Log-logistic: h(t) = (beta/alpha)(t/alpha)^(beta-1) / (1 + (t/alpha)^beta)
# H(t) = log(1 + (t/alpha)^beta)
# For uncensored data:
# loglik = n*log(beta) - n*beta*log(alpha) + (beta-1)*sum(log(t)) - 2*sum(log(1 + (t/alpha)^beta))

loglogistic_dist <- dfr_dist(
    rate = function(t, par, ...) {
        alpha <- par[[1]]
        beta <- par[[2]]
        z <- (t / alpha)^beta
        (beta / alpha) * (t / alpha)^(beta - 1) / (1 + z)
    },
    cum_haz_rate = function(t, par, ...) {
        alpha <- par[[1]]
        beta <- par[[2]]
        log(1 + (t / alpha)^beta)
    },
    score_fn = function(df, par, ...) {
        alpha <- par[[1]]
        beta <- par[[2]]
        t <- df$t
        n <- length(t)

        z <- (t / alpha)^beta
        log_t_alpha <- log(t / alpha)

        # d(loglik)/d(alpha) = -n*beta/alpha + 2*beta/alpha * sum(z/(1+z))
        score_alpha <- -n * beta / alpha + 2 * (beta / alpha) * sum(z / (1 + z))

        # d(loglik)/d(beta) = n/beta - n*log(alpha) + sum(log(t)) - 2*sum(z*log(t/alpha)/(1+z))
        score_beta <- n / beta - n * log(alpha) + sum(log(t)) -
                      2 * sum(z * log_t_alpha / (1 + z))

        c(score_alpha, score_beta)
    }
)

# Simulate log-logistic data using inverse CDF
# F(t) = 1 - 1/(1 + (t/alpha)^beta), so T = alpha * ((1-U)/U)^(1/beta)
set.seed(789)
alpha_true <- 2
beta_true <- 3
u <- runif(100)
loglog_times <- alpha_true * ((1 - u) / u)^(1 / beta_true)
loglog_df <- data.frame(t = loglog_times, delta = rep(1, length(loglog_times)))

score_ll <- score(loglogistic_dist)
hess_ll <- hess_loglik(loglogistic_dist)

cat("Log-logistic distribution (alpha=2, beta=3):\n")
#> Log-logistic distribution (alpha=2, beta=3):
cat("  Score at true params:", round(score_ll(loglog_df, par = c(alpha_true, beta_true)), 4), "\n")
#>   Score at true params: 9.7502 3.0433
cat("  Hessian:\n")
#>   Hessian:
print(round(hess_ll(loglog_df, par = c(alpha_true, beta_true)), 4))
#>          [,1]     [,2]
#> [1,] -83.5178   4.6166
#> [2,]   4.6166 -15.5408
```

### Log-Normal Distribution

Another common AFT model with a non-monotone hazard:

$$h(t;\mu,\sigma) = \frac{\phi(z)}{t\sigma\Phi( - z)},\quad z = \frac{\log t - \mu}{\sigma}$$

where $\phi$ and $\Phi$ are the standard normal PDF and CDF. The
cumulative hazard is $H(t) = - \log\left( \Phi( - z) \right)$.

``` r
# Note: dnorm/pnorm don't work with femtograd's AD objects,
# so log-normal uses pure numerical differentiation
lognormal_dist <- dfr_dist(
    rate = function(t, par, ...) {
        mu <- par[1]
        sigma <- par[2]
        z <- (log(t) - mu) / sigma
        dnorm(z) / (t * sigma * pnorm(-z))
    }
    # No cum_haz_rate or score_fn - uses numerical fallback
)

# Test with simulated log-normal data
set.seed(321)
mu_true <- 1
sigma_true <- 0.5
lognorm_times <- rlnorm(100, meanlog = mu_true, sdlog = sigma_true)
lognorm_df <- data.frame(t = lognorm_times, delta = rep(1, 100))

ll_ln <- loglik(lognormal_dist)
cat("Log-normal distribution (mu=1, sigma=0.5):\n")
#> Log-normal distribution (mu=1, sigma=0.5):
cat("  Log-likelihood at true params:", round(ll_ln(lognorm_df, par = c(mu_true, sigma_true)), 3), "\n")
#>   Log-likelihood at true params: -168.255
```

## Complete Distribution Reference

| Distribution | Parameters              | Hazard h(t)                   | Analytical Score |
|--------------|-------------------------|-------------------------------|------------------|
| Exponential  | λ (rate)                | λ                             | ✓ Simple         |
| Weibull      | k (shape), σ (scale)    | (k/σ)(t/σ)^(k-1)              | ✓                |
| Gompertz     | a (scale), b (shape)    | a·exp(bt)                     | ✓                |
| Log-logistic | α (scale), β (shape)    | (β/α)(t/α)^((β-1)/(1+(t/α))β) | ✓                |
| Log-normal   | μ (location), σ (scale) | φ(z)/(tσΦ(-z))                | Complex          |

The distributions with analytical score functions will use the efficient
hybrid AD approach automatically when `femtograd` is installed.

## Next Steps

- **[`vignette("failure_rate")`](https://queelius.github.io/dfr_dist/articles/failure_rate.md)** -
  Introduction to DFR distributions and hazard-based thinking
- **[`vignette("custom_distributions")`](https://queelius.github.io/dfr_dist/articles/custom_distributions.md)** -
  Complete guide to creating your own distributions
- **[`vignette("reliability_engineering")`](https://queelius.github.io/dfr_dist/articles/reliability_engineering.md)** -
  Real-world applications in reliability analysis
- **[`vignette("getting_started")`](https://queelius.github.io/dfr_dist/articles/getting_started.md)** -
  Quick 5-minute introduction
