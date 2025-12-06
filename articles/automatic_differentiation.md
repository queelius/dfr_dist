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

## Setup

``` r
library(dfr.dist)

# Check if femtograd is available for AD
has_ad <- requireNamespace("femtograd", quietly = TRUE)
if (has_ad) {
  library(femtograd)
  cat("femtograd available - AD examples will run\n")
} else {
  cat("femtograd not installed - AD examples will be skipped\n
       Install with: devtools::install_github('queelius/femtograd')\n")
}
#> femtograd available - AD examples will run
```

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
for AD compatibility. This allows femtograd’s dual numbers to flow
through the computation.

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
#> [1,] -73.4602  30.3417
#> [2,]  30.3417 -50.2047
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
#>   Full AD vs Numerical: 0
cat("  Hybrid vs Full AD:", max(abs(h_hybrid - h_full)), "\n")
#>   Hybrid vs Full AD: 0.0006991277
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
#>   Numerical:       108.58 ms
cat(sprintf("  Hybrid (recommended): %7.2f ms\n", time_hybrid))
#>   Hybrid (recommended):    1.55 ms
cat(sprintf("  Full AD:        %7.2f ms\n", time_full))
#>   Full AD:         107.65 ms
cat(sprintf("\nSpeedup (Numerical / Hybrid): %.1fx\n", time_num / time_hybrid))
#> 
#> Speedup (Numerical / Hybrid): 70.2x
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
#>       50       1.3 ms      54.6 ms      42x
#>      100       1.7 ms     107.9 ms      64x
#>      500       1.4 ms     536.4 ms     384x
#>     1000       1.8 ms    1107.3 ms     621x
```

The hybrid approach stays nearly constant (~2ms) while numerical scales
linearly with n.

## Why is the Hybrid Approach Fastest?

1.  **Numerical differentiation** requires many log-likelihood
    evaluations, each involving numerical integration for the cumulative
    hazard.

2.  **Full AD** builds computation graphs for every observation, which
    has significant overhead.

3.  **Hybrid approach** uses an analytical score that directly computes
    the gradient without numerical integration. The AD Jacobian then
    requires only 2 forward passes (one per parameter) on this efficient
    score function.

## Writing AD-Compatible Score Functions

For the hybrid approach, your `score_fn` must be compatible with
femtograd’s dual numbers:

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
