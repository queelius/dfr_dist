
  - [R package `dfr.dist`: dynamic failure rate (DFR)
    distributions](#r-package-dfrdist-dynamic-failure-rate-dfr-distributions)
      - [Installation](#installation)

<!-- README.md is generated from README.Rmd. Please edit that file -->

## R package `dfr.dist`: dynamic failure rate (DFR) distributions

<!-- badges: start -->

<!-- badges: end -->

An R package for working with models in survival analysis in which the
distribution is parameterized by a very flexible failure rate function
(any function that satisfies properties like being non-negative,
integrating to infinity over the domain, and having a support of `(0,
Inf)`.

### Installation

You can install the development version of `dfr.dist` from
[GitHub](https://github.com/queelius/dfr_dist) with:

``` r
# install.packages("devtools")
devtools::install_github("queelius/dfr_dist")
```
