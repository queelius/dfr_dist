# Compute Hessian using forward-over-reverse AD

For when you have only the objective function (not analytical gradient).
Uses dual numbers wrapping value objects.

## Usage

``` r
ad_hessian_fwd_rev(f, par)
```

## Arguments

- f:

  A function f(par) returning a scalar

- par:

  Numeric vector of parameters

## Value

Hessian matrix (p x p)
