# Method for obtaining the cdf of a `dfr_dist` object.

Method for obtaining the cdf of a `dfr_dist` object.

## Usage

``` r
# S3 method for class 'dfr_dist'
cdf(x, ...)
```

## Arguments

- x:

  The object to obtain the cdf of.

- ...:

  Additional arguments to pass into the `cum_haz` constructor.

## Value

A function that computes the cdf of the distribution. It accepts `t`,
the time at which to compute the cdf, `par`, the parameters of the
distribution, `log.p` argument that determines whether to compute the
log of the cdf, `lower.limit`, whether to compute the lower limit (F(t))
or upper limit (S(t) = 1-F(t)). Finally, it passes any additional
arguments `...` to the `rate` function of the `dfr_dist` object `x`.
