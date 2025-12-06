# Log-likelihood method for `dfr_dist` objects

Returns a function that computes the log-likelihood of the data given
the distribution parameters. The log-likelihood for survival data is:

## Usage

``` r
# S3 method for class 'dfr_dist'
loglik(model, ...)
```

## Arguments

- model:

  The `dfr_dist` object

- ...:

  Additional arguments to pass to the hazard and cumulative hazard

## Value

A function that computes the log-likelihood. It accepts: - `df`: A data
frame with observation times and censoring indicators - `par`: The
parameters of the distribution - `...`: Additional arguments passed to
internal functions

## Details

For exact observations (uncensored): log(f(t)) = log(h(t)) - H(t) For
right-censored observations: log(S(t)) = -H(t)

where h(t) is the hazard function, H(t) is the cumulative hazard, f(t) =
h(t)\*S(t) is the pdf, and S(t) = exp(-H(t)) is the survival function.
