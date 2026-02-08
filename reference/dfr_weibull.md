# Weibull Distribution (Power-Law Hazard)

Creates a DFR distribution with Weibull hazard function. The Weibull is
extremely versatile: it can model increasing (wear-out), decreasing
(infant mortality), or constant (exponential) failure rates.

## Usage

``` r
dfr_weibull(shape = NULL, scale = NULL)
```

## Arguments

- shape:

  Shape parameter (k). Controls hazard behavior: k \< 1: decreasing
  hazard (infant mortality) k = 1: constant hazard (exponential) k \> 1:
  increasing hazard (wear-out)

- scale:

  Scale parameter (sigma). Controls time scale.

## Value

A `dfr_dist` object with analytical rate, cumulative hazard, and score
function.

## Details

The Weibull distribution has:

- Hazard: \\h(t) = (k/\sigma)(t/\sigma)^{k-1}\\

- Cumulative hazard: \\H(t) = (t/\sigma)^k\\

- Survival: \\S(t) = e^{-(t/\sigma)^k}\\

- Characteristic life (63.2% failure): \\\sigma\\

## Reliability Interpretation

- Shape \< 1: Infant mortality (burn-in failures, defects)

- Shape = 1: Random failures (reduces to exponential)

- Shape = 2: Rayleigh distribution (linear hazard increase)

- Shape \> 2: Accelerating wear-out (fatigue, corrosion)

## B-Life Calculation

The B10 life (10% failure quantile) is commonly used in reliability:
\\B10 = \sigma \cdot (-\log(0.9))^{1/k}\\

## Examples

``` r
# Bearing with wear-out failure (shape > 1)
bearing <- dfr_weibull(shape = 2.5, scale = 50000)

# Hazard increases with time
h <- hazard(bearing)
h(10000)  # hazard at 10k hours
#> [1] 4.472136e-06
h(40000)  # much higher at 40k hours
#> [1] 3.577709e-05

# B10 life calculation
Q <- inv_cdf(bearing)
B10 <- Q(0.10)  # 10% failure quantile

# Fit to test data with right-censoring
set.seed(123)
test_data <- data.frame(
  t = pmin(rweibull(100, shape = 2.5, scale = 50000), 30000),
  delta = as.integer(rweibull(100, shape = 2.5, scale = 50000) <= 30000)
)
solver <- fit(dfr_weibull())
result <- solver(test_data, par = c(2, 40000))
coef(result)
#> [1]     6.091528 39999.995178
```
