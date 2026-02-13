# Classic Survival Distribution Constructors

Convenient constructors for commonly-used survival distributions. Each
provides the complete specification (rate, cum_haz_rate, score_fn, and
where practical, hess_fn) for optimal performance.

## Left-Censoring Note

The analytical score and Hessian functions provided by these
constructors assume event indicators in {0, 1} (right-censored and exact
observations). For left-censored data (delta = -1), these functions are
not applicable and the package automatically falls back to numerical
differentiation via
[`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) and
[`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
through the log-likelihood, which handles all censoring types correctly.
