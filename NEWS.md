# flexhaz 0.5.2

* Prototype constructors (`dfr_exponential()`, `dfr_weibull()`,
  `dfr_gompertz()`, `dfr_loglogistic()`) now prepend a family-specific
  S3 subclass (e.g. `"dfr_exponential"`) to the class chain. This
  enables downstream packages to dispatch on component type via
  `inherits()` rather than maintaining parallel lookup tables. Existing
  `inherits(x, "dfr_dist")` checks are unaffected (additive change).

# flexhaz 0.5.1

* Initial CRAN release.
* Added JOSS paper (`paper.md`, `paper.bib`).
* Added community guidelines (`CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`).
* Added Zenodo archive (DOI: 10.5281/zenodo.18757479).

# flexhaz 0.5.0

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

# flexhaz 0.1.0

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
