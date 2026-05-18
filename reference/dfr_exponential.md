# Exponential Distribution (Constant Hazard)

Creates a DFR distribution with constant failure rate (exponential). The
exponential distribution is "memoryless" - the hazard does not depend on
time, making it appropriate for random failures unrelated to age.

## Usage

``` r
dfr_exponential(lambda = NULL)
```

## Arguments

- lambda:

  Rate parameter (failure rate). If NULL, must be provided when calling
  methods or fitting. Must be positive.

## Value

A `dfr_dist` object of class `c("dfr_exponential", "dfr_dist", ...)`
with analytical rate, cumulative hazard, and score function.

## Details

The exponential distribution has:

- Hazard: \\h(t) = \lambda\\

- Cumulative hazard: \\H(t) = \lambda t\\

- Survival: \\S(t) = e^{-\lambda t}\\

- Mean time to failure: \\1/\lambda\\

## Reliability Interpretation

Use exponential for:

- Electronic components during useful life (random failures)

- Systems with redundancy where failures are independent

- As a baseline model to test against more complex alternatives

## Examples

``` r
# Component with MTBF of 1000 hours (lambda = 0.001)
comp <- dfr_exponential(lambda = 0.001)

# Survival probability at 500 hours
S <- surv(comp)
S(500)  # ~60.6%
#> [1] 0.6065307

# Fit to failure data
set.seed(42)
failures <- data.frame(t = rexp(50, rate = 0.001), delta = 1)
solver <- fit(comp)
result <- solver(failures, par = c(0.002))
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
coef(result)  # Should be close to 0.001
#> [1] 0.0008808454
```
