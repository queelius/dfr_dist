# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

`dfr.dist` is an R package for working with dynamic failure rate (DFR)
distributions in survival analysis. The package parameterizes
distributions using flexible failure rate (hazard) functions that can
depend on time and any set of predictors/covariates.

## Development Commands

``` r
devtools::document()      # Update NAMESPACE and .Rd files from roxygen2
devtools::build()         # Build package tarball
devtools::check()         # Run R CMD check
devtools::install()       # Install locally
devtools::test()          # Run tests
devtools::build_vignettes()    # Build vignettes
pkgdown::build_site()          # Build package website
```

Always run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
after modifying roxygen2 comments. NAMESPACE is auto-generated - never
edit manually.

## Core Architecture

### The `dfr_dist` Object

The central abstraction in `R/dfr_dist.R`:

``` r
dfr_dist(rate, par = NULL, eps = 0.01,
         ob_col = "t", delta_col = "delta",
         cum_haz_rate = NULL, score_fn = NULL)
```

- `rate`: Function computing hazard rate h(t, par, …)
- `par`: Distribution parameters (can be NULL if unknown)
- `eps`: Epsilon for numerical integration/sampling
- `ob_col`, `delta_col`: Column names for survival data
- `cum_haz_rate`: Optional analytical H(t, par) for AD gradient
  computation
- `score_fn`: Optional analytical score for AD Hessian via Jacobian

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

| Method                                                                  | Returns             | Formula                                                  |
|-------------------------------------------------------------------------|---------------------|----------------------------------------------------------|
| `hazard()`                                                              | h(t, par, …)        | Direct rate function                                     |
| [`cum_haz()`](https://queelius.github.io/dfr_dist/reference/cum_haz.md) | H(t, par, …)        | ∫₀ᵗ h(u) du (numerical)                                  |
| `surv()`                                                                | S(t, par, …)        | exp(-H(t))                                               |
| `cdf()`                                                                 | F(t, par, …)        | 1 - S(t)                                                 |
| [`pdf()`](https://rdrr.io/r/grDevices/pdf.html)                         | f(t, par, …)        | h(t) × S(t)                                              |
| `inv_cdf()`                                                             | t given F(t)        | Uses [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) |
| `sampler()`                                                             | Generates samples   | Inverse CDF sampling                                     |
| `params()`                                                              | Resolved parameters | Handles NA/NULL defaults                                 |

### Likelihood Model Interface

The `dfr_dist` class implements `likelihood_model` from the
likelihood.model package:

| Method                                                                                    | Returns                   | Usage                                  |
|-------------------------------------------------------------------------------------------|---------------------------|----------------------------------------|
| [`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html)           | Log-likelihood function   | `ll <- loglik(dist); ll(df, par)`      |
| [`score()`](https://queelius.github.io/likelihood.model/reference/score.html)             | Score function (gradient) | See fallback chain below               |
| [`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html) | Hessian of log-likelihood | See fallback chain below               |
| [`fit()`](https://generics.r-lib.org/reference/fit.html)                                  | MLE solver                | `solver <- fit(dist); solver(df, par)` |

### Automatic Differentiation Integration

The package optionally integrates with `femtograd` for exact
gradient/Hessian computation:

**Score function fallback chain:** 1. Analytical `score_fn` if provided
→ exact 2. AD gradient via `femtograd` if `cum_haz_rate` provided →
exact 3. Numerical gradient via
[`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html) →
approximate

**Hessian fallback chain (hybrid approach):** 1. AD Jacobian of
analytical `score_fn` via `femtograd` → exact 2. Numerical Hessian via
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html) →
approximate

The hybrid approach (analytical gradient + AD Jacobian) is often more
practical than full AD because: - Gradients are easier to derive
analytically than Hessians - Forward-mode AD Jacobian of a p-dimensional
gradient costs O(p) passes - Works even when full AD through
log-likelihood is impractical

**Example with analytical score:**

``` r
exp_dist <- dfr_dist(
    rate = function(t, par, ...) rep(par[1], length(t)),
    cum_haz_rate = function(t, par, ...) par[1] * t,
    score_fn = function(df, par, ...) {
        c(sum(df$delta) / par[1] - sum(df$t))
    }
)
# Hessian computed via AD Jacobian of score_fn
H <- hess_loglik(exp_dist)
hess <- H(df, par = c(1))  # Uses femtograd if available
```

### Survival Data Conventions

Data frames passed to likelihood methods use: - `t` column (or
`ob_col`): Observation times - `delta` column (or `delta_col`): Event
indicators - `1` = exact observation (failure) - `0` = right-censored
(survived past time t)

Log-likelihood contributions: - Exact: log L_i = log h(t_i) - H(t_i) -
Censored: log L_i = -H(t_i)

### Numerical Methods

**Cumulative hazard** uses
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) with
defaults: - subdivisions = 1000L, abs_tol = 1e-3 (overridable)

**Sampling** uses iterative rejection: 1. Start at t = 0, sample from
Exp(rate(t, par)) 2. Accept if runif(1) ≤ S(t, par), else increment t by
eps

## Package Ecosystem Integration

This package is designed to work with a family of related packages:

    algebraic.dist          # Generic distribution interface (imported)
        ↓
    dfr.dist               # This package - DFR distributions
        ↓
    likelihood.model       # Likelihood model interface (dfr_dist satisfies this)
        ↓
    algebraic.mle          # MLE algebra for hypothesis testing, model selection
        ↓
    numerical.mle          # Optimization: gradient ascent, Newton-Raphson, etc.

### Usage with MLE

The `dfr_dist` object satisfies the `likelihood_model` interface,
enabling:

``` r
# Fit DFR distribution to data using algebraic.mle
library(algebraic.mle)
mle_result <- mle(dfr_model, data, start_params)
```

### Series System Packages

Two packages handle **series systems** (systems that fail when any
component fails) with **masked failure data** (where the causing
component is uncertain):

- **likelihood.model.series.md**: Series system likelihood models with
  masked component cause data. Currently implements exponential and
  Weibull components. Its README explicitly mentions dfr_dist as a
  future integration for general series systems with arbitrary component
  hazard functions.

- **mdrelax**: Relaxed masking conditions for series systems
  (Weibull/exponential).

**Integration potential**: dfr_dist could serve as a component
distribution engine for series systems, enabling time-dependent and
covariate-dependent component hazards beyond the current
exponential/Weibull implementations. The architectural foundation exists
in likelihood.model.series.md’s `utils.R`
([`cum_haz()`](https://queelius.github.io/dfr_dist/reference/cum_haz.md),
`qcomp()`, `rcomp()` for arbitrary hazard functions).

### Core Packages

- **algebraic.dist**: Provides generics (hazard, cdf, pdf, surv, etc.)
- **likelihood.model**: Defines likelihood interface (loglik, score,
  hess_loglik)
- **algebraic.mle**: MLE algebra for inference
- **numerical.mle**: Numerical optimization methods
- **md.tools**: Masked data utilities (Boolean matrix encoding/decoding)

## File Structure

- `R/dfr_dist.R`: Main S3 class and all methods (~440 lines)
- `R/ad_utils.R`: AD helpers for femtograd integration
- `R/generic_functions.R`: Generic
  [`cum_haz()`](https://queelius.github.io/dfr_dist/reference/cum_haz.md)
  declaration
- `R/utils.R`: Helper `get_params()` for parameter resolution
- `R/reexports.R`: Re-exports from likelihood.model (loglik, score,
  etc.)
- `tests/testthat/test-ad_integration.R`: Tests for AD integration
- `tests/testthat/test-likelihood_model.R`: Comprehensive tests for
  likelihood interface
- `vignettes/failure_rate.Rmd`: Main tutorial

## Common Pitfalls

1.  **Parameter handling**: Always use `params(x, par)` to handle
    NA/NULL correctly
2.  **Integration warnings**: Cumulative hazard uses numerical
    integration - check for convergence
3.  **Sampling efficiency**: Rejection sampler can be slow for complex
    rate functions
4.  **Closure pattern**: Remember methods return functions, not values
    directly
5.  **AD compatibility**: When writing `score_fn` for AD Hessian, use
    `[[` indexing for parameters (e.g., `par[[1]]` not `par[1]`) to work
    with femtograd dual numbers
