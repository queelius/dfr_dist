# MLE solver for dfr_dist objects

Returns a solver function that fits a dfr_dist model to survival data
using maximum likelihood estimation. The solver uses gradient-based
optimization (BFGS by default) with analytical or numerical gradients.

## Usage

``` r
# S3 method for class 'dfr_dist'
fit(object, ...)
```

## Arguments

- object:

  A `dfr_dist` object

- ...:

  Additional arguments passed to the log-likelihood, score, and Hessian
  constructors (e.g., integration parameters for cum_haz)

## Value

A solver function that accepts:

- `df`: Data frame with observation times and censoring indicators

- `par`: Initial parameter values (uses object's params if NULL)

- `method`: Optimization method (default "BFGS")

- `control`: Control parameters for optim()

- `...`: Additional arguments passed to likelihood functions

## Details

The solver returns a `fisher_mle` object (from likelihood.model)
containing:

- Parameter estimates

- Log-likelihood value at MLE

- Variance-covariance matrix (from Hessian)

- Convergence status

Use methods like [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`summary()`](https://rdrr.io/r/base/summary.html) on the result.

## Examples

``` r
if (FALSE) { # \dontrun{
# Exponential distribution
exp_dist <- dfr_dist(
  rate = function(t, par, ...) rep(par[1], length(t)),
  par = c(lambda = 1)
)

# Simulate data
set.seed(42)
df <- data.frame(t = rexp(100, rate = 2), delta = 1)

# Fit model
solver <- fit(exp_dist)
result <- solver(df, par = c(1))
summary(result)
confint(result)
} # }
```
