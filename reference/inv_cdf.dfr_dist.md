# Method for obtaining the quantile (inverse cdf) of an object.

Method for obtaining the quantile (inverse cdf) of an object.

## Usage

``` r
# S3 method for class 'dfr_dist'
inv_cdf(x, ...)
```

## Arguments

- x:

  The object to obtain the inverse cdf of.

- ...:

  Additional arguments to pass into `cdf` constructor.

## Value

A function that computes the quantile of the distribution. It accepts
`p`, the probability at which to compute the quantile, `par`, the
parameters of the distribution, and `...`, any additional arguments to
pass into the constructed cdf.
