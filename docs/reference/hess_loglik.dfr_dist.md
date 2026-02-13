# Hessian of log-likelihood for dfr_dist

Returns a function that computes the Hessian matrix of the
log-likelihood. Uses user-provided Hessian function if available,
otherwise falls back to numerical differentiation via numDeriv::hessian.

## Usage

``` r
# S3 method for class 'dfr_dist'
hess_loglik(model, ...)
```

## Arguments

- model:

  A dfr_dist object

- ...:

  Additional arguments passed to loglik

## Value

A function that computes the Hessian matrix
