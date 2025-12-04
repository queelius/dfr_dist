---
output:
  github_document:
    toc: true
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# dfr.dist

<!-- badges: start -->
[![R-CMD-check](https://github.com/queelius/dfr_dist/workflows/R-CMD-check/badge.svg)](https://github.com/queelius/dfr_dist/actions)
<!-- badges: end -->

**Dynamic Failure Rate Distributions for Survival Analysis**

The `dfr.dist` package provides a flexible framework for specifying survival
distributions through their **hazard (failure rate) functions**. Instead of
choosing from a fixed catalog of distributions, you directly specify the hazard
function, giving complete control over time-varying failure patterns.

## Features

- **Flexible hazard specification**: Define any hazard function h(t, par, ...)
- **Complete distribution interface**: Automatic computation of survival, CDF, PDF, quantiles
- **Likelihood model support**: Log-likelihood, score, Hessian for MLE
- **Censoring support**: Handle exact and right-censored survival data
- **Integration with ecosystem**: Works with `algebraic.dist`, `likelihood.model`, `algebraic.mle`

## Installation

Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("queelius/dfr_dist")
```

## Quick Start


``` r
library(dfr.dist)
library(algebraic.dist)
```

### Define a distribution via its hazard function


``` r
# Exponential distribution: constant hazard h(t) = lambda
exp_dist <- dfr_dist(
  rate = function(t, par, ...) rep(par[1], length(t)),
  par = c(lambda = 0.5)
)

# All distribution functions are automatically available
S <- surv(exp_dist)
S(2)  # Survival probability at t=2
#> [1] 0.3678794
```

### Maximum Likelihood Estimation


``` r
# Simulate failure times
set.seed(42)
times <- rexp(50, rate = 1)
df <- data.frame(t = times, delta = 1)

# Create distribution without fixed parameters
dist <- dfr_dist(
  rate = function(t, par, ...) rep(par[1], length(t))
)

# Fit via MLE
solver <- fit(dist)
result <- solver(df, par = c(0.5), method = "BFGS")
params(result)
#> [1] 0.8808457
```

### Custom hazard functions

Model complex failure patterns like bathtub curves (infant mortality + useful life + wear-out):


``` r
# h(t) = a*exp(-b*t) + c + d*t^k
bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    par[1] * exp(-par[2] * t) + par[3] + par[4] * t^par[5]
  },
  par = c(a = 1, b = 2, c = 0.02, d = 0.001, k = 2)
)

h <- hazard(bathtub)
curve(sapply(x, h), 0, 15, xlab = "Time", ylab = "Hazard rate",
      main = "Bathtub hazard curve")
```

<div class="figure">
<img src="man/figures/README-bathtub-1.png" alt="plot of chunk bathtub" width="100%" />
<p class="caption">plot of chunk bathtub</p>
</div>

## Mathematical Background

For a lifetime $T$, the hazard function is:
$$h(t) = \frac{f(t)}{S(t)}$$

From the hazard, all other quantities follow:

| Function | Formula | Method |
|----------|---------|--------|
| Cumulative hazard | $H(t) = \int_0^t h(u) du$ | `cum_haz()` |
| Survival | $S(t) = e^{-H(t)}$ | `surv()` |
| CDF | $F(t) = 1 - S(t)$ | `cdf()` |
| PDF | $f(t) = h(t) \cdot S(t)$ | `pdf()` |

## Likelihood for Survival Data

For exact observations: $\log L = \log h(t) - H(t)$

For right-censored: $\log L = -H(t)$


``` r
# Mixed data with censoring
df <- data.frame(
  t = c(1, 2, 3, 4, 5),
  delta = c(1, 1, 0, 1, 0)  # 1 = exact, 0 = censored
)

ll <- loglik(dist)
ll(df, par = c(0.5))
#> [1] -9.579442
```

## Documentation

- [Vignette: Dynamic Failure Rate Distributions](https://queelius.github.io/dfr_dist/articles/failure_rate.html)
- [Function Reference](https://queelius.github.io/dfr_dist/reference/)

## Related Packages

- [`algebraic.dist`](https://github.com/queelius/algebraic.dist): Generic distribution interface
- [`likelihood.model`](https://github.com/queelius/likelihood.model): Likelihood model framework
- [`algebraic.mle`](https://github.com/queelius/algebraic.mle): MLE utilities
