# Sampling function for `dfr_dist` objects.

Since S(t,par) = exp(-cum_hz(t,par)), we can sample from the
distribution by letting t = 0 (or some other positive number if we want
to condition on T \> t_min), sampling from an exponential distribution
with `lambda = rate(t, par)`, and then rejecting the sample if
`runif(1) > S(t, par)`. If accepted, add that observation to the sample,
otherwise reject it, let `t = t + eps` where `eps` is some small number,
and repeat. We continue this process until we have `n` observations for
the sample.

## Usage

``` r
# S3 method for class 'dfr_dist'
sampler(x, ...)
```

## Arguments

- x:

  The object to obtain the sampler of.

- ...:

  Additional arguments to pass into the survival function

## Value

A function that samples from the distribution. It accepts `n`, the
number of samples to take, `t` is the time at which to start sampling,
`par` are the parameters of the distribution, and `eps` is the update
for numerical integration. Finally, we pass additional arguments `...`
into the hazard function.
