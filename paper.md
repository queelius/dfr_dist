# Summary

`flexhaz` is an R package for defining survival distributions through
their hazard (failure rate) function. The user writes the hazard, any R
function of time and parameters, and the package derives everything
else. Given a hazard $`h(t)`$, the full distribution follows by
numerical or analytical integration: cumulative hazard
$`H(t) = \int_0^t h(u)\,du`$, survival $`S(t) = \exp(-H(t))`$, density
$`f(t) = h(t) S(t)`$, quantiles, and random sampling. The package
provides a complete likelihood model interface (log-likelihood, score,
Hessian, and maximum likelihood estimation via `optim`) supporting
exact, right-censored, and left-censored observations. Model diagnostics
include Cox-Snell and Martingale residuals with Q-Q plots for
goodness-of-fit assessment.

The central abstraction is the `dfr_dist` S3 class. All distribution
methods return closures that accept time values and optional parameter
overrides, enabling late binding of parameters for MLE:

``` r

library(flexhaz)
bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    par[1] * exp(-par[2] * t) + par[3] + par[4] * t^par[5]
  },
  par = c(a = 0.5, b = 0.5, c = 0.01, d = 5e-5, k = 2.5)
)
S <- surv(bathtub)
S(10)  # survival probability at t = 10

# Fit to censored data (delta: 1=exact, 0=right-censored)
df <- data.frame(t = rexp(50, 0.02), delta = sample(0:1, 50, TRUE))
solver <- fit(bathtub)
result <- solver(df, par = c(0.3, 0.3, 0.02, 1e-4, 2))
```

# Statement of Need

Reliability engineers and biostatisticians who work with time-to-event
data frequently encounter failure patterns that do not conform to
standard parametric families. Bathtub curves, competing failure modes
with additive hazards, and time-covariate interactions all require
custom hazard specifications. Existing tools either restrict the user to
a catalog of built-in distributions or require specifying the density
and cumulative distribution function directly, quantities that are often
harder to derive or reason about than the hazard itself. Researchers in
reliability engineering and survival analysis need a framework that
accepts a hazard function as the sole specification and provides the
full distributional and inferential machinery automatically, including
MLE with censored data and residual diagnostics.

# State of the Field

Several R packages address parametric and semi-parametric survival
modeling, but none accept a user-defined hazard function as the sole
specification and derive the full distributional and inferential
machinery.

**Parametric frameworks.** The `survival` package \[@survival\] provides
Kaplan-Meier estimation, Cox proportional hazards, and parametric
accelerated failure time models, but does not support user-defined
hazard functions. `flexsurv` \[@flexsurv\] allows user-defined
distributions but requires the density or hazard together with its
cumulative form, and focuses on regression modeling rather than
standalone distribution objects with late-binding parameters. The `eha`
package \[@eha\] supports proportional hazards and AFT models but is
limited to built-in families.

**Non-parametric and spline methods.** `muhaz` \[@muhaz\] and `bshazard`
\[@bshazard\] estimate hazard functions non-parametrically via kernel
smoothing and B-splines, respectively, but produce exploratory estimates
rather than parametric distribution objects suitable for simulation or
MLE. `survPen` \[@survPen\] fits penalized spline hazard models with
automatic smoothing, and `rstpm2` \[@rstpm2\] implements Royston-Parmar
models using restricted cubic splines on the log cumulative hazard. Both
target data-driven flexible fitting rather than user-specified
mechanistic hazard functions.

`flexhaz` fills the gap between these approaches. It accepts a single R
function (the hazard) and derives the complete distribution, likelihood
model, and diagnostics. The package integrates with the `algebraic.dist`
\[@algebraicdist\] distribution interface and the `likelihood.model`
\[@likelihoodmodel\] framework, yielding `fisher_mle` result objects
with standard `coef`, `vcov`, `confint`, and `summary` methods.

# Software Design

The `dfr_dist` constructor accepts a hazard rate function and optional
analytical cumulative hazard, score, and Hessian functions. The
resulting S3 object inherits from `likelihood_model`, `univariate_dist`,
and `dist` classes. All distribution methods (`hazard`, `surv`, `cdf`,
`density`, `inv_cdf`, `sampler`) return closures that accept time `t`,
an optional parameter vector `par`, and additional arguments. This
closure-returning pattern separates the distribution specification from
parameter values, enabling the same object to serve as both a
parameterized distribution for simulation and a template for MLE.
Because additional arguments propagate through `...` in every closure,
the same architecture supports covariate-dependent hazards:
`h(t, par, x, ...)` naturally yields covariate-conditional survival,
density, and likelihood functions without special-case machinery.

The design follows a three-level optimization paradigm. At Level 1, the
user supplies only the hazard function; cumulative hazard, score, and
Hessian are all computed numerically. At Level 2, adding an analytical
cumulative hazard eliminates numerical integration during likelihood
evaluation. At Level 3, supplying analytical score and Hessian functions
yields production-quality MLE with fast, accurate derivatives. Built-in
constructors for exponential, Weibull, Gompertz, and log-logistic
distributions provide Level 3 implementations. When analytical
derivatives handle only exact and right-censored observations, the
package automatically falls back to
[`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) and
[`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
\[@numDeriv\] for left-censored data, ensuring correctness without user
intervention.

The `fit` method returns a solver closure:
`solver <- fit(dist); result <- solver(df, par)`. The result is a
`fisher_mle` object from the `likelihood.model` package, providing
`coef`, `vcov`, `confint`, `logLik`, and `summary` methods.

# Research Impact Statement

`flexhaz` is the foundation layer for an in-development family of
packages addressing masked failure data in series systems. `serieshaz`
composes `dfr_dist` components into series system distributions via
additive hazards, and `maskedhaz` builds likelihood models for masked
component-cause data with arbitrary hazard functions. A companion
package `maskedcauses`, currently under CRAN review, provides
closed-form solutions for exponential and Weibull components (798 tests,
99.45% coverage); `maskedhaz` results are validated against these
analytical solutions to verify the numerical approach. The full
ecosystem is available on r-universe at
<https://queelius.r-universe.dev>.

# AI Usage Disclosure

Claude Code (Anthropic) was used to assist with code refactoring, test
generation, and drafting of documentation including this manuscript. All
generated content was reviewed, edited, and validated by the author. The
package design, mathematical formulations, and research decisions are
entirely the author’s own work.

# Acknowledgements

The author thanks the R Core Team for the R statistical computing
environment \[@R\] and the developers of `numDeriv`, `algebraic.dist`,
and `likelihood.model` for foundational infrastructure.

# References
