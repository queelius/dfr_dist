# dfr.dist 0.5.0

## Changes

* Added helper distribution constructors: `dfr_exponential()`, `dfr_weibull()`,
  `dfr_gompertz()`, `dfr_loglogistic()` with analytical hazard, cumulative
  hazard, score, and Hessian functions where available.
* Added diagnostic methods: `residuals()` (Cox-Snell and Martingale),
  `plot()` (survival, hazard, cumulative hazard), `qqplot_residuals()`.
* Added `density()` method (alias for pdf).
* Added `assumptions()` method for listing model assumptions.
* Added `kaplan_meier()` internal utility for empirical survival estimation.
* Added support for left-censored observations (delta = -1).
* Improved numerical stability in log-likelihood computation.
* Removed femtograd dependency — users supply their own derivative functions
  via `score_fn` and `hess_fn`, or the package falls back to numDeriv.

# dfr.dist 0.1.0

Initial release.

## Features

* `dfr_dist()` constructor for creating distributions from hazard functions
* Complete distribution interface: `hazard()`, `cum_haz()`, `surv()`, `cdf()`, `pdf()`, `inv_cdf()`, `sampler()`
* Likelihood model interface: `loglik()`, `score()`, `hess_loglik()`, `fit()`
* Support for right-censored survival data
* Numerical integration for cumulative hazard
* Inverse CDF sampling

## Examples

Supports modeling complex hazard patterns:

* Exponential (constant hazard)
* Weibull (monotonic hazard)
* Bathtub curves (infant mortality + useful life + wear-out)
* Proportional hazards with covariates
