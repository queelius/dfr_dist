# Compute gradient using femtograd reverse-mode AD

Compute gradient using femtograd reverse-mode AD

## Usage

``` r
ad_gradient(f, par)
```

## Arguments

- f:

  A function f(par) returning a scalar. The function should work with
  parameters accessed via `par[i]` indexing.

- par:

  Numeric vector of parameters

## Value

Numeric vector of partial derivatives
