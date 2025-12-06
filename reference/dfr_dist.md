# Constructor for `dfr_dist` objects

We assume that the hazard rate is a function of time and any other
predictors. We also assume that integrate(rate(t), 0, Inf) = infinity
and that the support is (0, Inf).

## Usage

``` r
dfr_dist(
  rate,
  par = NULL,
  eps = 0.01,
  ob_col = "t",
  delta_col = "delta",
  cum_haz_rate = NULL,
  score_fn = NULL
)
```

## Arguments

- rate:

  A function that computes the hazard rate at time `t`.

- par:

  The parameters of the distribution. Defaults to `NULL`, which means
  that the parameters are unknown.

- eps:

  The epsilon update for numerical integration. Defaults to 0.01.

- ob_col:

  The column name for observation times in data frames. Defaults to "t".

- delta_col:

  The column name for event indicators in data frames. Uses standard
  survival analysis convention: 1 = event observed (exact), 0 =
  right-censored. Defaults to "delta".

- cum_haz_rate:

  Optional analytical cumulative hazard function H(t, par). If provided,
  enables exact AD-based gradient computation. Should return the
  integral of rate from 0 to t.

- score_fn:

  Optional analytical score function score(df, par). If provided,
  enables exact AD-based Hessian computation via Jacobian of the score.
  Should return gradient vector.

## Value

A `dfr_dist` object that inherits from `likelihood_model`.
