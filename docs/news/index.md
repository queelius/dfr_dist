# Changelog

## dfr.dist 0.1.0

Initial release.

### Features

- [`dfr_dist()`](https://queelius.github.io/dfr_dist/reference/dfr_dist.md)
  constructor for creating distributions from hazard functions
- Complete distribution interface:
  [`hazard()`](https://queelius.github.io/algebraic.dist/reference/hazard.html),
  [`cum_haz()`](https://queelius.github.io/dfr_dist/reference/cum_haz.md),
  [`surv()`](https://queelius.github.io/algebraic.dist/reference/surv.html),
  [`cdf()`](https://queelius.github.io/algebraic.dist/reference/cdf.html),
  [`pdf()`](https://queelius.github.io/algebraic.dist/reference/pdf.html),
  [`inv_cdf()`](https://queelius.github.io/algebraic.dist/reference/inv_cdf.html),
  [`sampler()`](https://queelius.github.io/algebraic.dist/reference/sampler.html)
- Likelihood model interface:
  [`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html),
  [`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
  [`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html),
  [`fit()`](https://generics.r-lib.org/reference/fit.html)
- Support for right-censored survival data
- Numerical integration for cumulative hazard
- Inverse CDF sampling

### Examples

Supports modeling complex hazard patterns:

- Exponential (constant hazard)
- Weibull (monotonic hazard)
- Bathtub curves (infant mortality + useful life + wear-out)
- Proportional hazards with covariates
