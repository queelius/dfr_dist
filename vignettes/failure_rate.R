## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(dfr.dist)
library(algebraic.dist)

## -----------------------------------------------------------------------------
# Exponential distribution: constant hazard h(t) = lambda
exp_dist <- dfr_dist(
  rate = function(t, par, ...) rep(par[1], length(t)),
  par = c(lambda = 0.5)
)

print(exp_dist)
is_dfr_dist(exp_dist)

## -----------------------------------------------------------------------------
weibull_dist <- dfr_dist(
  rate = function(t, par, ...) {
    k <- par[1]      # shape
    sigma <- par[2]  # scale
    (k / sigma) * (t / sigma)^(k - 1)
  },
  par = c(shape = 2, scale = 3)
)

## -----------------------------------------------------------------------------
h <- hazard(exp_dist)
H <- cum_haz(exp_dist)

# Evaluate at specific times
h(1)      # Hazard at t=1
h(c(1, 2, 3))  # Vectorized

H(2)      # Cumulative hazard at t=2

## -----------------------------------------------------------------------------
S <- surv(exp_dist)
F <- cdf(exp_dist)

# Verify S(t) + F(t) = 1
t <- 2
c(survival = S(t), cdf = F(t), sum = S(t) + F(t))

