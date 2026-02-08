# Compute Kaplan-Meier and Nelson-Aalen Estimates

Internal function to compute non-parametric survival estimates.

## Usage

``` r
kaplan_meier(time, delta)
```

## Arguments

- time:

  Observation times

- delta:

  Event indicators (1 = event, 0 = censored)

## Value

List with time, surv, cumhaz
