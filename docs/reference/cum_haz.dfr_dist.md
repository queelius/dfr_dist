# Method for obtaining the hazard function of a `dfr_dist` object.

Method for obtaining the hazard function of a `dfr_dist` object.

## Usage

``` r
# S3 method for class 'dfr_dist'
cum_haz(x, ...)
```

## Arguments

- x:

  The object to obtain the hazard function of.

- ...:

  Additional arguments to pass into the `integrate` function.

## Value

A function that computes the hazard function of the distribution. It
accepts `t`, the time at which to compute the hazard function, and
`par`, the parameters of the distribution. If `par` is `NULL`, then the
parameters of the `dfr_dist` object `x` are used. Finally, it passes any
additional arguments `...` to the `rate` function of the `dfr_dist`
object `x`.
