# Constructor for `dfr_dist` objects

We assume that the hazard rate is a function of time and any other
predictors. We also assume that integrate(rate(t), 0, Inf) = infinity
and that the support is (0, Inf).

## Usage

``` r
dfr_dist(
  rate,
  par = NULL,
  ob_col = "t",
  delta_col = "delta",
  cum_haz_rate = NULL,
  score_fn = NULL,
  hess_fn = NULL
)
```

## Arguments

- rate:

  A function that computes the hazard rate at time `t`.

- par:

  The parameters of the distribution. Defaults to `NULL`, which means
  that the parameters are unknown.

- ob_col:

  The column name for observation times in data frames. Defaults to "t".

- delta_col:

  The column name for event indicators in data frames. Uses standard
  survival analysis convention: 1 = event observed (exact), 0 =
  right-censored, -1 = left-censored. Defaults to "delta".

- cum_haz_rate:

  Optional analytical cumulative hazard function H(t, par). If provided,
  used for faster exact cumulative hazard computation instead of
  numerical integration. Should return the integral of rate from 0 to t.

- score_fn:

  Optional score function (gradient of log-likelihood). Signature:
  score_fn(df, par, ob_col, delta_col, ...) returning a numeric vector.
  The ob_col and delta_col arguments indicate which columns in df
  contain observation times and event indicators. If NULL, falls back to
  numerical gradient via numDeriv::grad. Analytical score functions that
  only handle delta in {0, 1} are automatically bypassed when
  left-censored data (delta = -1) is present.

- hess_fn:

  Optional Hessian function (second derivatives of log-likelihood).
  Signature: hess_fn(df, par, ob_col, delta_col, ...) returning a
  matrix. The ob_col and delta_col arguments indicate which columns in
  df contain observation times and event indicators. If NULL, falls back
  to numerical Hessian via numDeriv::hessian. Analytical Hessian
  functions that only handle delta in {0, 1} are automatically bypassed
  when left-censored data (delta = -1) is present.

## Value

A `dfr_dist` object that inherits from `likelihood_model`.
