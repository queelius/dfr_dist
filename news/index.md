# Changelog

## flexhaz 0.5.0

### Changes

- Added helper distribution constructors:
  [`dfr_exponential()`](https://queelius.github.io/flexhaz/reference/dfr_exponential.md),
  [`dfr_weibull()`](https://queelius.github.io/flexhaz/reference/dfr_weibull.md),
  [`dfr_gompertz()`](https://queelius.github.io/flexhaz/reference/dfr_gompertz.md),
  [`dfr_loglogistic()`](https://queelius.github.io/flexhaz/reference/dfr_loglogistic.md)
  with analytical hazard, cumulative hazard, score, and Hessian
  functions where available.
- Added diagnostic methods:
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) (Cox-Snell and
  Martingale), [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  (survival, hazard, cumulative hazard),
  [`qqplot_residuals()`](https://queelius.github.io/flexhaz/reference/qqplot_residuals.md).
- Added [`density()`](https://rdrr.io/r/stats/density.html) method
  (alias for pdf).
- Added
  [`assumptions()`](https://queelius.github.io/likelihood.model/reference/assumptions.html)
  method for listing model assumptions.
- Added
  [`kaplan_meier()`](https://queelius.github.io/flexhaz/reference/kaplan_meier.md)
  internal utility for empirical survival estimation.
- Added support for left-censored observations (delta = -1).
- Improved numerical stability in log-likelihood computation.
- Removed femtograd dependency — users supply their own derivative
  functions via `score_fn` and `hess_fn`, or the package falls back to
  numDeriv.

## flexhaz 0.1.0

Initial release.

### Features

- [`dfr_dist()`](https://queelius.github.io/flexhaz/reference/dfr_dist.md)
  constructor for creating distributions from hazard functions
- Complete distribution interface:
  [`hazard()`](https://queelius.github.io/algebraic.dist/reference/hazard.html),
  [`cum_haz()`](https://queelius.github.io/flexhaz/reference/cum_haz.md),
  [`surv()`](https://queelius.github.io/algebraic.dist/reference/surv.html),
  [`cdf()`](https://queelius.github.io/algebraic.dist/reference/cdf.html),
  [`pdf()`](https://rdrr.io/r/grDevices/pdf.html),
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
