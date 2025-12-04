# Method for obtaining the pdf of a `dfr_dist` object.

Method for obtaining the pdf of a `dfr_dist` object.

## Usage

``` r
# S3 method for class 'dfr_dist'
pdf(x, ...)
```

## Arguments

- x:

  The object to obtain the pdf of.

- ...:

  Additional arguments to pass.

## Value

A function that computes the pdf of the distribution. It accepts `y`,
the value at which to compute the pdf, `t`, the time at which to compute
the pdf, `par` is the parameters of the distribution, and `log`
determines whether to compute the log of the pdf. Finally, it passes any
additional arguments `...` to the `rate` function of the `dfr_dist`
object `x`.
