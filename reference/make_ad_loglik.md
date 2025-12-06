# Create AD-aware log-likelihood function for dfr_dist

Returns a function that computes log-likelihood using femtograd objects,
enabling automatic differentiation.

## Usage

``` r
make_ad_loglik(model, df)
```

## Arguments

- model:

  A dfr_dist object

- df:

  Data frame with observations

## Value

A function f(par) that works with femtograd value/dual objects
