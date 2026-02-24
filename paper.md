---
title: 'flexhaz: Hazard-First Survival Distributions for R'
tags:
  - R
  - survival analysis
  - hazard function
  - reliability engineering
  - maximum likelihood estimation
authors:
  - name: Alexander Towell
    orcid: 0000-0001-6443-9897
    affiliation: 1
affiliations:
  - name: Southern Illinois University Edwardsville
    index: 1
date: "24 February 2026"
bibliography: paper.bib
---

# Summary

Standard parametric survival models force analysts to choose from a short list of hazard shapes: constant (exponential), power-law (Weibull), exponential growth (Gompertz). But real failure patterns rarely cooperate. Ceramic capacitors exhibit wear-out that accelerates faster than any power law. Software systems follow bathtub-shaped crash rates. Post-surgical patients face risk that drops sharply then slowly climbs again. Each demands a hazard function that standard families cannot express.

`flexhaz` is an R package that takes a different approach: the user writes the hazard function --- any R function of time and parameters --- and the package derives everything else. Given a hazard $h(t)$, the full distribution follows by numerical or analytical integration: cumulative hazard $H(t) = \int_0^t h(u)\,du$, survival $S(t) = \exp(-H(t))$, density $f(t) = h(t) S(t)$, quantiles, and random sampling. The package also provides a complete likelihood model interface --- log-likelihood, score, Hessian, and maximum likelihood estimation (MLE) via `optim` --- supporting exact, right-censored, and left-censored observations. Model diagnostics include Cox-Snell and Martingale residuals with Q-Q plots for goodness-of-fit assessment.

The central abstraction is the `dfr_dist` S3 class. All distribution methods return closures that accept time values and optional parameter overrides, enabling late binding of parameters for MLE:

```r
library(flexhaz)
bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    par[1] * exp(-par[2] * t) + par[3] + par[4] * t^par[5]
  },
  par = c(a = 0.5, b = 0.5, c = 0.01, d = 5e-5, k = 2.5)
)
S <- surv(bathtub)
S(10)  # survival probability at t = 10
solver <- fit(bathtub)
result <- solver(data, par = c(0.3, 0.3, 0.02, 1e-4, 2))
```

# Statement of Need

Reliability engineers and biostatisticians who work with time-to-event data frequently encounter failure patterns that do not conform to standard parametric families. Bathtub curves, competing failure modes with additive hazards, and time-covariate interactions all require custom hazard specifications. Existing tools either restrict the user to a catalog of built-in distributions or require specifying the density and cumulative distribution function directly --- quantities that are often harder to derive or reason about than the hazard itself. Researchers in reliability engineering and survival analysis need a framework that accepts a hazard function as the sole specification and provides the full distributional and inferential machinery automatically, including MLE with censored data and residual diagnostics.

# State of the Field

Several R packages address parametric and semi-parametric survival modeling. The `survival` package [@survival] provides the foundational infrastructure for survival analysis in R, including Kaplan-Meier estimation, Cox proportional hazards regression, and parametric accelerated failure time models, but does not support user-defined hazard functions as a distribution specification mechanism. `flexsurv` [@flexsurv] extends parametric survival modeling by allowing user-defined distributions, but requires the user to supply a density or hazard function together with its cumulative form; it focuses on regression modeling with covariate effects on distribution parameters rather than providing a standalone distribution object with late-binding parameters and a likelihood model interface. The `eha` package [@eha] supports parametric proportional hazards and accelerated failure time models for demographic and epidemiological applications but is limited to built-in distribution families. Non-parametric hazard estimation is available through `muhaz` [@muhaz], which uses kernel smoothing, and `bshazard` [@bshazard], which uses B-splines within a generalized linear mixed model framework; these produce exploratory hazard estimates but do not yield parametric distribution objects suitable for simulation, extrapolation, or MLE. The `survPen` package [@survPen] fits penalized spline hazard models with automatic smoothing parameter selection, targeting flexible regression rather than user-specified mechanistic hazard functions. The `rstpm2` package [@rstpm2] implements Royston-Parmar flexible parametric survival models using restricted cubic splines on the log cumulative hazard, with support for time-varying effects; like `flexsurv`, it targets data-driven flexible fitting rather than mechanistic hazard specification.

`flexhaz` fills the gap between these approaches. It accepts a single R function --- the hazard --- and derives the complete distribution, likelihood model, and diagnostics. A three-level optimization paradigm lets users start with just the hazard (Level 1), add an analytical cumulative hazard for speed (Level 2), and supply analytical score and Hessian functions for production-quality MLE (Level 3). Built-in constructors for exponential, Weibull, Gompertz, and log-logistic distributions provide Level 3 implementations out of the box. The package integrates with the `algebraic.dist` [@algebraicdist] distribution interface and the `likelihood.model` [@likelihoodmodel] framework, yielding `fisher_mle` result objects with standard `coef`, `vcov`, `confint`, and `summary` methods.

# Software Design

The `dfr_dist` constructor accepts a hazard rate function and optional analytical cumulative hazard, score, and Hessian functions. The resulting S3 object inherits from `likelihood_model`, `univariate_dist`, and `dist` classes. All distribution methods (`hazard`, `surv`, `cdf`, `density`, `inv_cdf`, `sampler`) return closures that accept time `t`, an optional parameter vector `par`, and additional arguments. This closure-returning pattern separates the distribution specification from parameter values, enabling the same object to serve as both a parameterized distribution for simulation and a template for MLE.

Derivative computation follows a two-tier fallback: if the user supplies `score_fn` or `hess_fn`, those are used directly; otherwise, the package falls back to `numDeriv::grad` and `numDeriv::hessian` [@numDeriv]. Analytical derivatives that handle only exact and right-censored observations (delta in {0, 1}) are automatically bypassed when left-censored data (delta = -1) is present, ensuring correct numerical fallback without user intervention.

# Research Impact Statement

`flexhaz` serves as the foundation layer for a family of packages addressing masked failure data in series systems: `serieshaz` composes `dfr_dist` components into series system distributions via additive hazards, and `maskedhaz` builds likelihood models for masked component-cause data with arbitrary hazard functions. These packages support ongoing research on model selection for masked series systems [@towell2024masked], where the ability to define arbitrary component hazard functions --- rather than being restricted to exponential or Weibull families --- enables cross-validation between analytical and numerical approaches. The masked-data research program includes a companion package `maskedcauses` with closed-form solutions for exponential and Weibull components (798 tests, 99.45% coverage), against which `maskedhaz` results are validated. The full ecosystem is available on r-universe at <https://queelius.r-universe.dev>.

# AI Usage Disclosure

Claude Code (Anthropic) was used to assist with code refactoring, test generation, and drafting of documentation including this manuscript. All generated content was reviewed, edited, and validated by the author. The package design, mathematical formulations, and research decisions are entirely the author's own work.

# Acknowledgements

The author thanks the R Core Team for the R statistical computing environment [@R] and the developers of `numDeriv`, `algebraic.dist`, and `likelihood.model` for foundational infrastructure.

# References
