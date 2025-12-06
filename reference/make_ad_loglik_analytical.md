# Create AD-compatible log-likelihood using analytical cumulative hazard

Create AD-compatible log-likelihood using analytical cumulative hazard

## Usage

``` r
make_ad_loglik_analytical(model, df)
```

## Arguments

- model:

  A dfr_dist object with cum_haz_rate defined

- df:

  Data frame with observations

## Value

A function f(par) compatible with femtograd
