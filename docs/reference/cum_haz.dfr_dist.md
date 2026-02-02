# Method for obtaining the cumulative hazard function of a `dfr_dist` object.

Method for obtaining the cumulative hazard function of a `dfr_dist`
object.

## Usage

``` r
# S3 method for class 'dfr_dist'
cum_haz(x, ...)
```

## Arguments

- x:

  The object to obtain the cumulative hazard function of.

- ...:

  Additional arguments to pass into the `integrate` function (only used
  when no analytical cum_haz_rate is provided).

## Value

A function that computes the cumulative hazard H(t) of the distribution.
It accepts `t`, the time at which to compute the cumulative hazard, and
`par`, the parameters of the distribution. If `par` is `NULL`, then the
parameters of the `dfr_dist` object `x` are used. Finally, it passes any
additional arguments `...` to the rate function.

## Details

If the `dfr_dist` object has an analytical `cum_haz_rate` function, that
is used directly for fast, exact computation. Otherwise, numerical
integration of the hazard function is performed.
