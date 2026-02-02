# Compute Hessian via AD Jacobian of gradient

Uses the hybrid approach: takes an analytical (or AD) gradient function
and computes its Jacobian using forward-mode AD to get the Hessian.

## Usage

``` r
ad_hessian(score_fn, par)
```

## Arguments

- score_fn:

  A function score_fn(par) returning the gradient vector. Should accept
  either numeric vector or list of value/dual objects.

- par:

  Numeric vector of parameters

## Value

Hessian matrix (p x p)
