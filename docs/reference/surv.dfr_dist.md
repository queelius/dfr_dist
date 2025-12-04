# Method for obtaining the survival function of a `dfr_dist` object.

Method for obtaining the survival function of a `dfr_dist` object.

## Usage

``` r
# S3 method for class 'dfr_dist'
surv(x, ...)
```

## Arguments

- x:

  The object to obtain the survival function of.

- ...:

  Additional arguments to pass into the `cum_haz` constructor.

## Value

A function that computes the survival function of the distribution. It
accepts `t`, the time at which to compute the survival, `par`, the
parameters of the distribution, `log.p` argument that determines whether
to compute the log of the survival, and it passes any additional
arguments into the `rate` function of the `dfr_dist` object `x`.
