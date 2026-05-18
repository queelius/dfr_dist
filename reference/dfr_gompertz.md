# Gompertz Distribution (Exponential Growth Hazard)

Creates a DFR distribution with Gompertz hazard function. The Gompertz
models exponentially increasing failure rate, often used for biological
aging and wear-out processes that accelerate over time.

## Usage

``` r
dfr_gompertz(a = NULL, b = NULL)
```

## Arguments

- a:

  Initial hazard rate at t=0. Must be positive.

- b:

  Growth rate of the hazard. Must be positive.

## Value

A `dfr_dist` object of class `c("dfr_gompertz", "dfr_dist", ...)` with
analytical rate, cumulative hazard, and score function.

## Details

The Gompertz distribution has:

- Hazard: \\h(t) = a \cdot e^{bt}\\

- Cumulative hazard: \\H(t) = (a/b)(e^{bt} - 1)\\

- Survival: \\S(t) = \exp(-(a/b)(e^{bt} - 1))\\

## Reliability Interpretation

Use Gompertz for:

- Aging systems where failure rate grows exponentially

- Biological mortality (human lifespans)

- Corrosion/degradation with accelerating kinetics

When b is small, Gompertz approximates exponential early in life. As b
increases, wear-out acceleration becomes more pronounced.

## Examples

``` r
# Aging system: initial hazard 0.001, doubling every 1000 hours
# b = log(2)/1000 gives doubling time of 1000
system <- dfr_gompertz(a = 0.001, b = log(2)/1000)

# Hazard at various ages
h <- hazard(system)
h(0)      # 0.001 (initial)
#> [1] 0.001
h(1000)   # 0.002 (doubled)
#> [1] 0.002
h(2000)   # 0.004 (quadrupled)
#> [1] 0.004

# Survival probability
S <- surv(system)
S(5000)   # probability of surviving 5000 hours
#> [1] 3.774076e-20
```
