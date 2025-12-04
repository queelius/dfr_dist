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
