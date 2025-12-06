# Hessian of log-likelihood for dfr_dist

Returns a function that computes the Hessian matrix of the
log-likelihood. If analytical score is provided, uses AD Jacobian for
exact computation. Otherwise falls back to numerical differentiation.

## Usage

``` r
# S3 method for class 'dfr_dist'
hess_loglik(model, ...)
```

## Arguments

- model:

  A dfr_dist object

- ...:

  Additional arguments passed to score

## Value

A function that computes the Hessian matrix
