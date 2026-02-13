# Score function (gradient of log-likelihood) for dfr_dist

Returns a function that computes the score (gradient of log-likelihood)
with respect to parameters. Uses user-provided score function if
available, otherwise falls back to numerical differentiation via
numDeriv::grad.

## Usage

``` r
# S3 method for class 'dfr_dist'
score(model, ...)
```

## Arguments

- model:

  A dfr_dist object

- ...:

  Additional arguments passed to loglik

## Value

A function that computes the score vector
