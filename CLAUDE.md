# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

`flexhaz` is an R package for working with dynamic failure rate (DFR)
distributions in survival analysis. The package parameterizes
distributions using flexible failure rate (hazard) functions that can
depend on time and any set of predictors/covariates.

## Development Commands

``` r

devtools::document()      # Update NAMESPACE and .Rd files from roxygen2
devtools::build()         # Build package tarball
devtools::check()         # Run R CMD check
devtools::install()       # Install locally
devtools::test()          # Run all tests
devtools::test(filter = "derivatives")  # Run specific test file
devtools::build_vignettes()    # Build vignettes
pkgdown::build_site()          # Build package website
pkgdown::deploy_to_branch()    # Deploy docs to GitHub Pages (gh-pages branch)
covr::package_coverage()       # Check test coverage
```

Always run `devtools::document()` after modifying roxygen2 comments.
NAMESPACE is auto-generated - never edit manually.

## Dependencies

**Required (Imports):** - `stats`, `numDeriv` - numerical methods -
`algebraic.dist` - distribution interface generics -
`likelihood.model` - likelihood model interface - `generics` - generic
function infrastructure

**Optional (Suggests):** - `testthat`, `knitr`, `rmarkdown` - testing
and documentation

## Core Architecture

### The `dfr_dist` Object

The central abstraction in `R/dfr_dist.R`:

``` r

dfr_dist(rate, par = NULL, eps = 0.01,
         ob_col = "t", delta_col = "delta",
         cum_haz_rate = NULL, score_fn = NULL,
         hess_fn = NULL)
```

- `rate`: Function computing hazard rate h(t, par, …)
- `par`: Distribution parameters (can be NULL if unknown)
- `eps`: Epsilon for numerical integration/sampling
- `ob_col`, `delta_col`: Column names for survival data
- `cum_haz_rate`: Optional analytical H(t, par) for performance
- `score_fn`: Optional score function (gradient of log-likelihood)
- `hess_fn`: Optional Hessian function (second derivatives of
  log-likelihood)

Inherits from `"univariate_dist"` and `"dist"` classes, implementing the
`algebraic.dist` interface.

### Key Constraints on Hazard Functions

All failure rate functions must satisfy: 1. **Non-negativity**: h(t, …)
≥ 0 for all t 2. **Infinite cumulative hazard**: lim(t→∞) H(t) = ∞,
where H(t) = ∫₀ᵗ h(u) du

These are not enforced in code but are required for mathematical
correctness.

### Closure-Returning Pattern

All distribution methods follow a two-step pattern - understanding this
is critical for working with the package:

1.  Method is called on `dfr_dist` object, returning a closure
2.  Closure accepts time `t`, parameters `par`, and additional args
    `...`

``` r

# Create distribution
d <- dfr_dist(rate = function(t, par, ...) par, par = 1)

# Get hazard function (returns closure)
h <- hazard(d)

# Use the closure
h(t = 5)              # Uses default par = 1
h(t = 5, par = 2)     # Overrides with par = 2
```

This enables flexible parameter handling and late binding of parameters
(useful for MLE).

### Distribution Methods

| Method | Returns | Formula |
|----|----|----|
| [`hazard()`](https://queelius.github.io/algebraic.dist/reference/hazard.html) | h(t, par, …) | Direct rate function |
| [`cum_haz()`](https://queelius.github.io/flexhaz/reference/cum_haz.md) | H(t, par, …) | ∫₀ᵗ h(u) du (numerical) |
| [`surv()`](https://queelius.github.io/algebraic.dist/reference/surv.html) | S(t, par, …) | exp(-H(t)) |
| [`cdf()`](https://queelius.github.io/algebraic.dist/reference/cdf.html) | F(t, par, …) | 1 - S(t) |
| [`density()`](https://rdrr.io/r/stats/density.html) | f(t, par, …) | h(t) × S(t) |
| [`inv_cdf()`](https://queelius.github.io/algebraic.dist/reference/inv_cdf.html) | t given F(t) | Uses [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) |
| [`sampler()`](https://queelius.github.io/algebraic.dist/reference/sampler.html) | Generates samples | Inverse CDF sampling |
| [`params()`](https://queelius.github.io/algebraic.dist/reference/params.html) | Resolved parameters | Handles NA/NULL defaults |

### Likelihood Model Interface

The `dfr_dist` class implements `likelihood_model` from the
likelihood.model package:

| Method | Returns | Usage |
|----|----|----|
| [`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html) | Log-likelihood function | `ll <- loglik(dist); ll(df, par)` |
| [`score()`](https://queelius.github.io/likelihood.model/reference/score.html) | Score function (gradient) | `score_fn` → [`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) |
| [`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html) | Hessian of log-likelihood | `hess_fn` → [`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html) |
| [`fit()`](https://generics.r-lib.org/reference/fit.html) | MLE solver | `solver <- fit(dist); result <- solver(df, par)` → `fisher_mle` |
| [`assumptions()`](https://queelius.github.io/likelihood.model/reference/assumptions.html) | Model assumptions | Returns character vector of assumptions |

### Derivative Computation

The package uses simple 2-tier fallbacks for derivatives:

**Score:** `score_fn(df, par)` → `numDeriv::grad(loglik, par)`
**Hessian:** `hess_fn(df, par)` → `numDeriv::hessian(loglik, par)`

The user decides how to compute derivatives (hand-derive, use AD, etc.).
The package just accepts functions or falls back to numerical methods.

**Example with analytical score and Hessian:**

``` r

exp_dist <- dfr_dist(
    rate = function(t, par, ...) rep(par[1], length(t)),
    cum_haz_rate = function(t, par, ...) par[1] * t,
    score_fn = function(df, par, ...) {
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
        c(sum(delta) / par[1] - sum(df$t))
    },
    hess_fn = function(df, par, ...) {
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))
        matrix(-sum(delta == 1) / par[1]^2, 1, 1)
    }
)
```

### Survival Data Conventions

Data frames passed to likelihood methods use: - `t` column (or
`ob_col`): Observation times - `delta` column (or `delta_col`): Event
indicators - `1` = exact observation (failure) - `0` = right-censored
(survived past time t) - `-1` = left-censored (failed before time t)

Log-likelihood contributions: - Exact (δ=1): log L_i = log h(t_i) -
H(t_i) - Right-censored (δ=0): log L_i = -H(t_i) - Left-censored (δ=-1):
log L_i = log(1 - exp(-H(t_i)))

### Numerical Methods

**Cumulative hazard** uses
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) with
defaults: - subdivisions = 1000L, abs_tol = 1e-3 (overridable)

**Sampling** uses inverse CDF method: 1. Generate u ~ Uniform(0, 1) 2.
Return inv_cdf(u, par) via
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html)

## Package Ecosystem Integration

This package is designed to work with a family of related packages:

    algebraic.dist          # Generic distribution interface (imported)
        ↓
    flexhaz               # This package - DFR distributions
        ↓
    likelihood.model       # Likelihood model interface + fisher_mle class

### Usage with MLE

The `dfr_dist` object implements the `likelihood_model` interface and
has its own [`fit()`](https://generics.r-lib.org/reference/fit.html)
method:

``` r

# Create distribution
dist <- dfr_dist(rate = function(t, par, ...) rep(par[1], length(t)))

# Fit to data
solver <- fit(dist)
result <- solver(df, par = c(1), method = "BFGS")

# Result is a fisher_mle object with standard methods
coef(result)      # Parameter estimates
vcov(result)      # Variance-covariance matrix
confint(result)   # Confidence intervals
summary(result)   # Full summary
```

### Series System Packages

Two packages handle **series systems** (systems that fail when any
component fails) with **masked failure data** (where the causing
component is uncertain):

- **maskedcauses**: Series system likelihood models with masked
  component cause data. Currently implements exponential and Weibull
  components. Its README explicitly mentions dfr_dist as a future
  integration for general series systems with arbitrary component hazard
  functions.

- **mdrelax**: Relaxed masking conditions for series systems
  (Weibull/exponential).

**Integration potential**: dfr_dist could serve as a component
distribution engine for series systems, enabling time-dependent and
covariate-dependent component hazards beyond the current
exponential/Weibull implementations. The architectural foundation exists
in maskedcauses’s `utils.R`
([`cum_haz()`](https://queelius.github.io/flexhaz/reference/cum_haz.md),
`qcomp()`, `rcomp()` for arbitrary hazard functions).

### Core Packages

- **algebraic.dist**: Provides generics (hazard, cdf, pdf, surv, etc.)
- **likelihood.model**: Defines likelihood interface (loglik, score,
  hess_loglik) and fisher_mle class
- **numerical.mle**: Numerical optimization methods (optional, for
  advanced use)
- **md.tools**: Masked data utilities (Boolean matrix encoding/decoding)

## Helper Distribution Constructors

The package provides ready-to-use constructors for classic survival
distributions:

| Constructor | Hazard Pattern | Parameters |
|----|----|----|
| `dfr_exponential(lambda)` | Constant | lambda = failure rate |
| `dfr_weibull(shape, scale)` | Power-law | shape (k), scale (σ) |
| `dfr_gompertz(a, b)` | Exponential growth | a = initial hazard, b = growth rate |
| `dfr_loglogistic(alpha, beta)` | Non-monotonic | alpha = scale, beta = shape |

Each provides analytical `rate`, `cum_haz_rate`, and `score_fn`.
Exponential and Weibull also include analytical `hess_fn`:

``` r

# Quick creation with parameters
weib <- dfr_weibull(shape = 2, scale = 100)

# Template for fitting (parameters supplied later)
solver <- fit(dfr_weibull())
result <- solver(df, par = c(2, 100))
```

## Diagnostic Methods

| Method                         | Purpose                               |
|--------------------------------|---------------------------------------|
| `residuals(dist, data, type)`  | Cox-Snell or Martingale residuals     |
| `plot(dist, what)`             | Visualize survival, hazard, or cumhaz |
| `qqplot_residuals(dist, data)` | Q-Q plot for model validation         |

Cox-Snell residuals should follow Exp(1) if model is correct.

## File Structure

- `R/dfr_dist.R`: Main S3 class and all methods
- `R/distributions.R`: Helper constructors (dfr_exponential,
  dfr_weibull, etc.)
- `R/diagnostics.R`: Residuals, plot, and Q-Q plot methods
- `R/generic_functions.R`: Generic
  [`cum_haz()`](https://queelius.github.io/flexhaz/reference/cum_haz.md)
  declaration
- `R/utils.R`: Helper `get_params()` for parameter resolution
- `R/reexports.R`: Re-exports from likelihood.model and algebraic.dist
- `tests/testthat/test-dfr_dist.R`: Core distribution tests
- `tests/testthat/test-distributions.R`: Tests for helper constructors
- `tests/testthat/test-diagnostics.R`: Tests for diagnostic methods
- `tests/testthat/test-derivatives.R`: Tests for score/Hessian
  computation
- `tests/testthat/test-likelihood_model.R`: Likelihood interface tests
- `vignettes/failure_rate.Rmd`: Hazard-based modeling deep dive
- `vignettes/custom_distributions.Rmd`: How to create custom
  distributions
- `vignettes/reliability_engineering.Rmd`: Real-world applications
- `vignettes/custom_derivatives.Rmd`: Custom derivatives tutorial

## Common Pitfalls

1.  **Parameter handling**: Always use `params(x, par)` to handle
    NA/NULL correctly
2.  **Integration warnings**: Cumulative hazard uses numerical
    integration - check for convergence
3.  **Sampling efficiency**: Inverse CDF sampling uses uniroot - can be
    slow for complex rate functions
4.  **Closure pattern**: Remember methods return functions, not values
    directly
5.  **Derivative functions**: `score_fn(df, par, ...)` must return a
    numeric vector; `hess_fn(df, par, ...)` must return a matrix
