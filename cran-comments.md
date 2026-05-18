## R CMD check results

0 errors | 0 warnings | 1 note

* One NOTE: "unable to verify current time" (system clock check, benign).

## Patch (v0.5.1 -> v0.5.2)

Additive backwards-compatible change. The prototype constructors
dfr_exponential(), dfr_weibull(), dfr_gompertz(), and dfr_loglogistic()
now prepend a family-specific S3 subclass (e.g., "dfr_exponential",
"dfr_weibull") to the class chain. Downstream packages can dispatch on
component type via inherits() rather than maintaining parallel lookup
tables.

The change is purely additive: existing inherits(x, "dfr_dist") checks
and all distribution methods (hazard, surv, cdf, density, sampler, etc.)
continue to dispatch and behave identically.

## Coordinated submission

This is part of a coordinated 3-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- flexhaz 0.5.2 (this package)
- serieshaz 0.2.0
- maskedcauses 0.10.0

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (R-devel and R-release)
