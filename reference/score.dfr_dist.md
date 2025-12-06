# Score function (gradient of log-likelihood) for dfr_dist

Returns a function that computes the score (gradient of log-likelihood)
with respect to parameters. Uses analytical score if provided, otherwise
falls back to numerical differentiation.

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
