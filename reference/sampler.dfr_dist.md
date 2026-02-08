# Sampling function for `dfr_dist` objects.

Uses inverse CDF sampling: generates uniform random values and
transforms them through the quantile function (inverse CDF).

## Usage

``` r
# S3 method for class 'dfr_dist'
sampler(x, ...)
```

## Arguments

- x:

  The object to obtain the sampler of.

- ...:

  Additional arguments to pass into the inverse CDF constructor.

## Value

A function that samples from the distribution. It accepts `n`, the
number of samples to take, `par`, the parameters of the distribution,
and `...`, additional arguments passed to the quantile function.
