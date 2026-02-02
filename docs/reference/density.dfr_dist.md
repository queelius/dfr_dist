# Method for obtaining the density (pdf) of a `dfr_dist` object.

Method for obtaining the density (pdf) of a `dfr_dist` object.

## Usage

``` r
# S3 method for class 'dfr_dist'
density(x, ...)
```

## Arguments

- x:

  The object to obtain the density of.

- ...:

  Additional arguments to pass.

## Value

A function that computes the density of the distribution. It accepts
`t`, the time at which to compute the density, `par` is the parameters
of the distribution, and `log` determines whether to compute the log of
the density. Finally, it passes any additional arguments `...` to the
`rate` function of the `dfr_dist` object `x`.
