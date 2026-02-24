# Log-Logistic Distribution (Non-Monotonic Hazard)

Creates a DFR distribution with log-logistic hazard function. The
log-logistic has a non-monotonic hazard that increases then decreases,
useful for modeling processes with an initial risk that diminishes.

## Usage

``` r
dfr_loglogistic(alpha = NULL, beta = NULL)
```

## Arguments

- alpha:

  Scale parameter. Median lifetime when beta \> 1.

- beta:

  Shape parameter. Controls hazard shape: beta \<= 1: monotonically
  decreasing hazard beta \> 1: hazard increases to a peak then decreases

## Value

A `dfr_dist` object with analytical rate, cumulative hazard, and score
function.

## Details

The log-logistic distribution has:

- Hazard: \\h(t) = \frac{(\beta/\alpha)(t/\alpha)^{\beta-1}}{1 +
  (t/\alpha)^\beta}\\

- Cumulative hazard: \\H(t) = \log(1 + (t/\alpha)^\beta)\\

- Survival: \\S(t) = \frac{1}{1 + (t/\alpha)^\beta}\\

- Median: \\\alpha\\ (when beta \> 1)

## Reliability Interpretation

The log-logistic is useful when:

- Initial failures decrease after screening period

- Risk peaks early then declines (therapy response)

- Hazard is not monotonic throughout lifetime

The cumulative hazard has a closed form and is provided analytically.

## Examples

``` r
# Component with peak hazard around t = alpha
comp <- dfr_loglogistic(alpha = 1000, beta = 2)

# Non-monotonic hazard
h <- hazard(comp)
h(500)   # increasing phase
#> [1] 8e-04
h(1000)  # near peak
#> [1] 0.001
h(2000)  # decreasing phase
#> [1] 8e-04

# Survival function
S <- surv(comp)
S(1000)  # 50% survival at median (alpha)
#> [1] 0.5
```
