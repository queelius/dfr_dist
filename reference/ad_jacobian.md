# Compute Jacobian using femtograd forward-mode AD

For a function f: R^n -\> R^m, computes the m x n Jacobian matrix. Each
column is computed by a forward pass with tangent = 1 for that input.

## Usage

``` r
ad_jacobian(f, par)
```

## Arguments

- f:

  A function f(par) returning a numeric vector. The function should
  access parameters using `[[` indexing (e.g., `par[[1]]`) for AD
  compatibility, or accept a vector where femtograd ops are overloaded.

- par:

  Numeric vector of parameters

## Value

Jacobian matrix (output_dim x input_dim)
