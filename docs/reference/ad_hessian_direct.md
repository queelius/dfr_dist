# Compute Hessian directly using femtograd's forward-over-reverse

For when you have the objective function directly. Uses femtograd's
built-in hessian() function.

## Usage

``` r
ad_hessian_direct(f, par)
```

## Arguments

- f:

  A function f(par) returning a scalar. Should work with a list of value
  objects accessed via `[[` indexing.

- par:

  Numeric vector of parameters

## Value

Hessian matrix (p x p)
