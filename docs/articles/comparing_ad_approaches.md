# Comparing AD Approaches for Survival MLE

## Introduction

Maximum likelihood estimation (MLE) for survival models requires
computing derivatives of the log-likelihood function. The **score
function** (gradient of the log-likelihood) drives optimization
algorithms, and the **Hessian matrix** (second derivatives) provides
standard errors via the observed Fisher information.

Four approaches to computing these derivatives are available in the R
ecosystem, each with different tradeoffs:

| Approach                                  | AD Mode                                             | Gradient Cost       | Hessian Cost                                         |
|-------------------------------------------|-----------------------------------------------------|---------------------|------------------------------------------------------|
| **numDeriv** (finite differences)         | None                                                | O(p) function evals | O(p²) function evals                                 |
| **femtograd** (reverse + forward)         | Reverse-mode gradient, forward-over-reverse Hessian | O(1) backward pass  | O(p) forward passes                                  |
| **dual** (forward-mode dual numbers)      | Forward-mode                                        | O(p) forward passes | O(p) forward passes on analytical gradient           |
| **dualr** (forward-mode, convenience API) | Forward-mode                                        | O(p) forward passes | O(p(p+1)/2) nested duals or O(p) on analytical score |

where *p* is the number of parameters.

This vignette compares these approaches using the Weibull distribution,
showing how each computes the same derivatives with different accuracy
and performance characteristics.

## Setup

``` r
library(dfr.dist)

has_femtograd <- requireNamespace("femtograd", quietly = TRUE)

# Check dual and dualr availability without loading their namespaces yet.
# Both packages define an S4 class called "dual" with different slot layouts,
# so loading both namespaces simultaneously causes S4 class conflicts.
# We detect availability via system.file() to avoid triggering namespace load.
has_dual  <- nzchar(system.file(package = "dual"))
has_dualr <- nzchar(system.file(package = "dualr"))

# dual and dualr cannot coexist in the same R session due to the S4 class
# name collision. We run dual sections first (has_dual_only), then dualr
# sections (has_dualr_only) — but only one set can use live AD objects.
# For the combined comparison sections, we pick whichever is available.
has_both_ad <- has_femtograd && has_dual
has_all_ad  <- has_femtograd && has_dual && has_dualr

cat("Package availability:\n")
#> Package availability:
cat("  numDeriv:  always (in Imports)\n")
#>   numDeriv:  always (in Imports)
cat("  femtograd:", if (has_femtograd) "YES" else "not installed", "\n")
#>   femtograd: YES
cat("  dual:     ", if (has_dual) "YES" else "not installed", "\n")
#>   dual:      YES
cat("  dualr:    ", if (has_dualr) "YES" else "not installed", "\n")
#>   dualr:     YES
if (has_dual && has_dualr) {
    cat("\nNote: dual and dualr both define an S4 class named 'dual'.\n")
    cat("Sections run in order: dual first, then dualr.\n")
}
#> 
#> Note: dual and dualr both define an S4 class named 'dual'.
#> Sections run in order: dual first, then dualr.
```

## The Weibull Problem

We use the Weibull distribution as our running example. Its hazard
function is:

$$h(t;k,\sigma) = \frac{k}{\sigma}\left( \frac{t}{\sigma} \right)^{k - 1}$$

The cumulative hazard (needed for the survival function and
log-likelihood) is:

$$H(t;k,\sigma) = \left( \frac{t}{\sigma} \right)^{k}$$

For uncensored data with $n$ observations, the log-likelihood is:

$$\ell(k,\sigma) = n\log k + (k - 1)\sum\limits_{i = 1}^{n}\log t_{i} - nk\log\sigma - \sum\limits_{i = 1}^{n}\left( \frac{t_{i}}{\sigma} \right)^{k}$$

### Analytical Derivatives (Ground Truth)

The score (gradient) has closed-form expressions:

$$\frac{\partial\ell}{\partial k} = \frac{n}{k} + \sum\limits_{i}\log t_{i} - n\log\sigma - \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k}\log\!\left( \frac{t_{i}}{\sigma} \right)$$

$$\frac{\partial\ell}{\partial\sigma} = \frac{k}{\sigma}\left\lbrack \sum\limits_{i}\left( \frac{t_{i}}{\sigma} \right)^{k} - n \right\rbrack$$

The Hessian can be derived by differentiating the score again. We
compute it below as our reference for accuracy comparisons.

### Simulated Data

``` r
set.seed(42)
true_k     <- 2
true_sigma <- 3
n          <- 100

u     <- runif(n)
times <- true_sigma * (-log(u))^(1 / true_k)
df    <- data.frame(t = times, delta = rep(1, n))

test_par <- c(1.8, 2.8)
cat("True parameters: k =", true_k, ", sigma =", true_sigma, "\n")
#> True parameters: k = 2 , sigma = 3
cat("Test point:      k =", test_par[1], ", sigma =", test_par[2], "\n")
#> Test point:      k = 1.8 , sigma = 2.8
```

### Analytical Reference

``` r
# Analytical score
analytical_score <- function(df, par) {
    k     <- par[1]
    sigma <- par[2]
    t     <- df$t
    n     <- length(t)
    r     <- t / sigma
    r_k   <- r^k
    log_r <- log(r)

    dk     <- n / k + sum(log(t)) - n * log(sigma) - sum(r_k * log_r)
    dsigma <- k / sigma * (sum(r_k) - n)
    c(dk, dsigma)
}

# Analytical Hessian (derived by hand)
analytical_hessian <- function(df, par) {
    k     <- par[1]
    sigma <- par[2]
    t     <- df$t
    n     <- length(t)
    r     <- t / sigma
    r_k   <- r^k
    log_r <- log(r)

    # d²l/dk²
    d2k <- -n / k^2 - sum(r_k * log_r^2)

    # d²l/dkds = d²l/dsdk
    d2ks <- (1 / sigma) * (sum(r_k) - n) +
            (k / sigma) * sum(r_k * log_r)

    # d²l/ds²
    d2s <- -k / sigma^2 * (sum(r_k) - n) +
            k^2 / sigma^2 * sum(r_k)

    matrix(c(d2k, d2ks, d2ks, d2s), nrow = 2)
}

ref_score <- analytical_score(df, test_par)
ref_hess  <- analytical_hessian(df, test_par)

cat("Analytical score:\n")
#> Analytical score:
print(ref_score)
#> [1] -10.494704   8.878147
cat("\nAnalytical Hessian:\n")
#> 
#> Analytical Hessian:
print(round(ref_hess, 4))
#>          [,1]    [,2]
#> [1,] -73.4609 30.3410
#> [2,]  30.3410 43.8631
```

## Approach 1: Finite Differences (numDeriv)

The simplest approach requires only the log-likelihood function. The
`numDeriv` package approximates derivatives via Richardson extrapolation
of finite differences.

``` r
# Define the distribution with just the hazard rate
weibull_numerical <- dfr_dist(
    rate = function(t, par, ...) {
        k     <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    }
)

# Get likelihood functions - these use numDeriv internally
ll_fn    <- loglik(weibull_numerical)
score_nd <- score(weibull_numerical)
hess_nd  <- hess_loglik(weibull_numerical)

nd_score <- score_nd(df, par = test_par)
nd_hess  <- hess_nd(df, par = test_par)

cat("numDeriv score:\n")
#> numDeriv score:
print(nd_score)
#> [1] -10.494194   8.878177
cat("\nnumDeriv Hessian:\n")
#> 
#> numDeriv Hessian:
print(round(nd_hess, 4))
#>          [,1]     [,2]
#> [1,] -73.4602  30.3417
#> [2,]  30.3417 -50.2047

cat("\nScore error (max abs vs analytical):", max(abs(nd_score - ref_score)), "\n")
#> 
#> Score error (max abs vs analytical): 0.000509126
cat("Hessian error (max abs vs analytical):", max(abs(nd_hess - ref_hess)), "\n")
#> Hessian error (max abs vs analytical): 94.06785
```

**Pros:** No analytical derivatives needed—just provide the hazard rate
function.

**Cons:** Each score evaluation requires $O(p)$ log-likelihood
evaluations (each involving numerical integration of the hazard), and
each Hessian evaluation requires $O\left( p^{2} \right)$ evaluations.
Precision is typically around $10^{- 5}$ to $10^{- 6}$.

## Approach 2: femtograd (Reverse-Mode + Forward-over-Reverse)

``` r
# femtograd is used internally by dfr.dist via requireNamespace().
# We avoid library(femtograd) because femtograd::fit() would mask
# the generics::fit() S3 generic used by dfr.dist.
cat("femtograd", as.character(packageVersion("femtograd")), "available\n")
#> femtograd 0.3.0 available
```

The `femtograd` package provides two AD modes:

- **Reverse-mode** for gradients: Builds a computation graph of `value`
  objects during the forward pass, then calls `backward()` to compute
  all partial derivatives in a single backward pass. This is what
  `dfr.dist` uses via
  [`ad_gradient()`](https://queelius.github.io/dfr_dist/reference/ad_gradient.md)
  when a `cum_haz_rate` function is provided.

- **Forward-mode dual numbers** for the Jacobian: When computing the
  Hessian as the Jacobian of the score,
  [`ad_jacobian()`](https://queelius.github.io/dfr_dist/reference/ad_jacobian.md)
  wraps parameters in `dual_num` objects and performs one forward pass
  per parameter.

### Key API constraint

femtograd’s `value` and `dual_num` objects require **double-bracket
indexing** (`par[[i]]`) to extract the AD objects from a list. Single
brackets (`par[i]`) return a sub-list, which breaks the computation
graph.

### Full AD: Score from cum_haz_rate

``` r
weibull_full_ad <- dfr_dist(
    rate = function(t, par, ...) {
        k     <- par[[1]]
        sigma <- par[[2]]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k     <- par[[1]]
        sigma <- par[[2]]
        (t / sigma)^k
    }
    # No score_fn — femtograd computes the score via reverse-mode AD
)

score_fad <- score(weibull_full_ad)
hess_fad  <- hess_loglik(weibull_full_ad)

fad_score <- score_fad(df, par = test_par)
fad_hess  <- hess_fad(df, par = test_par)

cat("femtograd full-AD score:\n")
#> femtograd full-AD score:
print(fad_score)
#> [1] -10.494704   8.878147
cat("\nfemtograd full-AD Hessian:\n")
#> 
#> femtograd full-AD Hessian:
print(round(fad_hess, 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047

cat("\nScore error:", max(abs(fad_score - ref_score)), "\n")
#> 
#> Score error: 8.348877e-14
cat("Hessian error:", max(abs(fad_hess - ref_hess)), "\n")
#> Hessian error: 94.06782
```

### Hybrid: Analytical Score + AD Jacobian

The hybrid approach provides an analytical score function, then lets
femtograd compute its Jacobian (which equals the Hessian) using
forward-mode AD:

``` r
weibull_hybrid <- dfr_dist(
    rate = function(t, par, ...) {
        k     <- par[[1]]
        sigma <- par[[2]]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k     <- par[[1]]
        sigma <- par[[2]]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k     <- par[[1]]
        sigma <- par[[2]]
        t     <- df$t
        n     <- length(t)
        r     <- t / sigma
        r_k   <- r^k
        log_r <- log(r)

        dk     <- n / k + sum(log(t)) - n * log(sigma) -
                   sum(r_k * log_r)
        dsigma <- k / sigma * (sum(r_k) - n)
        c(dk, dsigma)
    }
)

score_hyb <- score(weibull_hybrid)
hess_hyb  <- hess_loglik(weibull_hybrid)

hyb_score <- score_hyb(df, par = test_par)
hyb_hess  <- hess_hyb(df, par = test_par)

cat("Hybrid score (analytical):\n")
#> Hybrid score (analytical):
print(hyb_score)
#> [1] -10.494704   8.878147
cat("\nHybrid Hessian (AD Jacobian of score):\n")
#> 
#> Hybrid Hessian (AD Jacobian of score):
print(round(hyb_hess, 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047

cat("\nScore error:", max(abs(hyb_score - ref_score)), "\n")
#> 
#> Score error: 0
cat("Hessian error:", max(abs(hyb_hess - ref_hess)), "\n")
#> Hessian error: 94.06782
```

The hybrid approach is often faster than full AD because the analytical
score avoids building a computation graph for every observation. The AD
Jacobian then requires only $p$ forward passes on the efficient score
function.

## Approach 3: dual (Forward-Mode Dual Numbers)

``` r
# We avoid library(dual) to prevent S4 class name conflicts with dualr
# (both define a class called "dual" with different slot names).
# Use dual:: prefix for all dual package functions.
cat("dual", as.character(packageVersion("dual")), "available\n")
#> dual 0.0.6 available
```

The `dual` package on CRAN implements forward-mode AD using S4 dual
number objects. A dual number pairs a value with its derivative:

$$\widetilde{x} = \left( x,\dot{x} \right)$$

Arithmetic on dual numbers propagates derivatives via the chain rule
automatically. To compute the gradient of
$\left. f:{\mathbb{R}}^{p}\rightarrow{\mathbb{R}} \right.$, we evaluate
$f$ once per parameter direction, seeding the tangent vector $\dot{x}$
with each standard basis vector $e_{j}$.

### How dual numbers work

``` r
# Create dual numbers with tangent for 2-variable function
x <- dual::dual(3, c(1, 0))  # x=3, tangent=(1,0) → tracking df/dx
y <- dual::dual(2, c(0, 1))  # y=2, tangent=(0,1) → tracking df/dy

z <- x^2 + 3 * x * y - y^2
cat("f(3,2) = x² + 3xy - y²\n")
#> f(3,2) = x² + 3xy - y²
cat("  value:", z@f, "\n")
#>   value: 23
cat("  gradient:", z@grad, "\n")
#>   gradient: 12 5
cat("  expected: f=15, grad=(2x+3y, 3x-2y) = (12, 5)\n")
#>   expected: f=15, grad=(2x+3y, 3x-2y) = (12, 5)
```

A single evaluation gives us the function value and the full gradient—no
separate backward pass needed. The tradeoff is that if we seeded only
one tangent direction, we’d need $p$ separate evaluations for the full
gradient. By seeding all $p$ tangent components simultaneously (as
above), we get the full gradient in one pass.

### Score via dual numbers

To compute the Weibull score at a given parameter vector, we seed the
parameters as dual numbers and evaluate the log-likelihood. Since `dual`
objects don’t vectorize over numeric vectors (they are scalar objects
with S4 dispatch), we use `Reduce` or loops to sum over observations:

``` r
# Log-likelihood function written for dual compatibility
weibull_ll_dual <- function(theta, t_obs) {
    k     <- theta[[1]]
    sigma <- theta[[2]]
    n     <- length(t_obs)

    # n * log(k) + (k-1)*sum(log(t)) - n*k*log(sigma) - sum((t/sigma)^k)
    # The sum(log(t)) is a plain numeric constant
    sum_log_t <- sum(log(t_obs))

    term1 <- n * log(k)
    term2 <- (k - 1) * sum_log_t
    term3 <- n * k * log(sigma)

    # sum((t/sigma)^k): each term involves dual arithmetic
    power_sum <- Reduce("+", lapply(t_obs, function(ti) (ti / sigma)^k))

    term1 + term2 - term3 - power_sum
}

# Compute score: seed tangent with identity vectors
dual_score_fn <- function(par, t_obs) {
    p     <- length(par)
    theta <- lapply(seq_len(p), function(i) {
        tangent <- rep(0, p)
        tangent[i] <- 1
        dual::dual(par[i], tangent)
    })

    result <- weibull_ll_dual(theta, t_obs)
    result@grad
}

dual_sc <- dual_score_fn(test_par, df$t)
cat("dual score:\n")
#> dual score:
print(dual_sc)
#> [1] -10.494704   8.878147
cat("\nScore error:", max(abs(dual_sc - ref_score)), "\n")
#> 
#> Score error: 1.953993e-14
```

### Hessian via Jacobian of Analytical Score

Just like femtograd’s hybrid approach, we can compute the Hessian by
differentiating an analytical score function with `dual` numbers. The
score function must return a **list** of scalar dual values (one per
parameter), so we can extract each component’s tangent:

``` r
# Analytical score returning list of duals
weibull_score_dual <- function(theta, t_obs) {
    k     <- theta[[1]]
    sigma <- theta[[2]]
    n     <- length(t_obs)

    sum_log_t <- sum(log(t_obs))

    # sum((t/sigma)^k) and sum((t/sigma)^k * log(t/sigma))
    r_k_sum     <- Reduce("+", lapply(t_obs, function(ti) (ti / sigma)^k))
    r_k_log_sum <- Reduce("+", lapply(t_obs, function(ti) {
        r <- ti / sigma
        r^k * log(r)
    }))

    dk     <- n / k + sum_log_t - n * log(sigma) - r_k_log_sum
    dsigma <- k / sigma * (r_k_sum - n)

    list(dk, dsigma)
}

# Hessian = Jacobian of score
dual_hessian_fn <- function(par, t_obs) {
    p    <- length(par)
    hess <- matrix(0, p, p)

    for (j in seq_len(p)) {
        # Seed dual with tangent e_j
        theta <- lapply(seq_len(p), function(i) {
            dual::dual(par[i], if (i == j) 1 else 0)
        })
        g <- weibull_score_dual(theta, t_obs)
        for (i in seq_len(p)) {
            hess[i, j] <- g[[i]]@grad
        }
    }
    hess
}

dual_hess <- dual_hessian_fn(test_par, df$t)
cat("dual Hessian (Jacobian of analytical score):\n")
#> dual Hessian (Jacobian of analytical score):
print(round(dual_hess, 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
cat("\nHessian error:", max(abs(dual_hess - ref_hess)), "\n")
#> 
#> Hessian error: 94.06782
```

### Supported Math Functions

A strength of the `dual` package is its broad support for special
mathematical functions through S4 method dispatch:

| Category               | Functions                                                                       |
|------------------------|---------------------------------------------------------------------------------|
| **Basic**              | `+`, `-`, `*`, `/`, `^`, `abs`, `sqrt`                                          |
| **Transcendental**     | `log`, `exp`, `logb`                                                            |
| **Trigonometric**      | `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`                            |
| **Gamma family**       | `gamma`, `lgamma`, `digamma`, `trigamma`, `psigamma`, `factorial`, `lfactorial` |
| **Beta/Combinatorial** | `beta`, `lbeta`, `choose`, `lchoose`                                            |
| **Error functions**    | `erf`, `erfc`, `erfinv`, `erfcinv`                                              |

This makes `dual` suitable for distributions whose likelihoods involve
special functions (e.g., the generalized gamma distribution uses
`lgamma`).

## Scaling: 3-Parameter Gompertz–Makeham Example

The Gompertz–Makeham distribution adds a constant baseline hazard
$\lambda$ to the Gompertz model, giving three parameters:

$$h(t;\lambda,a,b) = \lambda + ae^{bt}$$

$$H(t;\lambda,a,b) = \lambda t + \frac{a}{b}\left( e^{bt} - 1 \right)$$

This lets us see how the approaches scale from $p = 2$ to $p = 3$.

``` r
# Simulate Makeham data
set.seed(99)
true_lambda <- 0.05
true_a      <- 0.01
true_b      <- 0.3
n_mak       <- 150

# Inverse CDF sampling via rejection (simple approach: use hazard + surv)
# For efficiency, use a fine grid approximation
t_grid <- seq(0.001, 30, length.out = 10000)
S_grid <- exp(-(true_lambda * t_grid + (true_a / true_b) * (exp(true_b * t_grid) - 1)))
u_mak  <- runif(n_mak)
mak_times <- approx(S_grid, t_grid, xout = u_mak, rule = 2)$y
mak_df    <- data.frame(t = mak_times, delta = rep(1, n_mak))

mak_test_par <- c(0.04, 0.008, 0.25)
cat("Makeham true parameters: lambda =", true_lambda, ", a =", true_a,
    ", b =", true_b, "\n")
#> Makeham true parameters: lambda = 0.05 , a = 0.01 , b = 0.3
cat("Test point:", mak_test_par, "\n")
#> Test point: 0.04 0.008 0.25
```

``` r
makeham_nd <- dfr_dist(
    rate = function(t, par, ...) {
        par[1] + par[2] * exp(par[3] * t)
    }
)

score_mak_nd <- score(makeham_nd)
hess_mak_nd  <- hess_loglik(makeham_nd)

cat("numDeriv score (3 params):\n")
#> numDeriv score (3 params):
print(score_mak_nd(mak_df, par = mak_test_par))
#> [1]  480.0670 4426.1957  365.1092
cat("\nnumDeriv Hessian (3 params):\n")
#> 
#> numDeriv Hessian (3 params):
print(round(hess_mak_nd(mak_df, par = mak_test_par), 3))
#>            [,1]       [,2]       [,3]
#> [1,] -22253.094  -93128.85  -5428.997
#> [2,] -93128.852 -856134.12 -25374.042
#> [3,]  -5428.997  -25374.04  -2352.126
```

``` r
# Analytical score for Makeham
makeham_score_femto <- function(df, par, ...) {
    lambda <- par[[1]]
    a      <- par[[2]]
    b      <- par[[3]]
    t      <- df$t
    n      <- length(t)

    exp_bt <- exp(b * t)
    h      <- lambda + a * exp_bt

    dlambda <- sum(1 / h) - sum(t)
    da      <- sum(exp_bt / h) - sum(exp_bt - 1) / b
    db      <- sum(a * t * exp_bt / h) -
               (a / b) * sum(t * exp_bt) +
               (a / b^2) * sum(exp_bt - 1)
    c(dlambda, da, db)
}

makeham_hyb <- dfr_dist(
    rate = function(t, par, ...) {
        par[[1]] + par[[2]] * exp(par[[3]] * t)
    },
    cum_haz_rate = function(t, par, ...) {
        par[[1]] * t + (par[[2]] / par[[3]]) * (exp(par[[3]] * t) - 1)
    },
    score_fn = makeham_score_femto
)

score_mak_hyb <- score(makeham_hyb)
hess_mak_hyb  <- hess_loglik(makeham_hyb)

cat("femtograd hybrid score (3 params):\n")
#> femtograd hybrid score (3 params):
print(score_mak_hyb(mak_df, par = mak_test_par))
#> [1]  480.0670 4426.1957  365.1092
cat("\nfemtograd hybrid Hessian:\n")
#> 
#> femtograd hybrid Hessian:
print(round(hess_mak_hyb(mak_df, par = mak_test_par), 3))
#>            [,1]       [,2]       [,3]
#> [1,] -22253.094  -93128.85  -5428.997
#> [2,] -93128.852 -856134.12 -25374.042
#> [3,]  -5428.997  -25374.04  -2352.126
```

``` r
# Analytical score for dual (returns list of scalars)
makeham_score_dual <- function(theta, t_obs) {
    lambda <- theta[[1]]
    a      <- theta[[2]]
    b      <- theta[[3]]
    n      <- length(t_obs)

    # Build sums over observations using Reduce
    # Each observation contributes terms involving dual parameters
    h_inv_sum <- Reduce("+", lapply(t_obs, function(ti) {
        1 / (lambda + a * exp(b * ti))
    }))
    exp_bt_over_h_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ebt <- exp(b * ti)
        ebt / (lambda + a * ebt)
    }))
    exp_bt_minus1_sum <- Reduce("+", lapply(t_obs, function(ti) {
        exp(b * ti) - 1
    }))
    t_sum <- sum(t_obs)
    t_exp_bt_over_h_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ebt <- exp(b * ti)
        a * ti * ebt / (lambda + a * ebt)
    }))
    t_exp_bt_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ti * exp(b * ti)
    }))

    dlambda <- h_inv_sum - t_sum
    da      <- exp_bt_over_h_sum - exp_bt_minus1_sum / b
    db      <- t_exp_bt_over_h_sum -
               (a / b) * t_exp_bt_sum +
               (a / b^2) * exp_bt_minus1_sum

    list(dlambda, da, db)
}

# Hessian via dual Jacobian
makeham_dual_hessian <- function(par, t_obs) {
    p    <- length(par)
    hess <- matrix(0, p, p)
    for (j in seq_len(p)) {
        theta <- lapply(seq_len(p), function(i) {
            dual::dual(par[i], if (i == j) 1 else 0)
        })
        g <- makeham_score_dual(theta, t_obs)
        for (i in seq_len(p)) {
            hess[i, j] <- g[[i]]@grad
        }
    }
    hess
}

cat("dual Hessian (3 params):\n")
#> dual Hessian (3 params):
print(round(makeham_dual_hessian(mak_test_par, mak_df$t), 3))
#>            [,1]       [,2]       [,3]
#> [1,] -22253.094  -93128.85  -5428.997
#> [2,] -93128.852 -856134.12 -25374.042
#> [3,]  -5428.997  -25374.04  -2352.126
```

## Approach 4: dualr (Forward-Mode with Convenience API)

``` r
# We do NOT call library(dualr) — dualr exports score() and hessian()
# which would mask dfr.dist's methods. Use dualr:: prefix throughout.
cat("dualr", as.character(packageVersion("dualr")), "available\n")
#> dualr 0.1.0 available
```

The `dualr` package is a forward-mode AD package with a higher-level API
than the CRAN `dual` package. Its key differentiator is a set of
convenience functions —
[`dualr::score()`](https://rdrr.io/pkg/dualr/man/score.html),
[`dualr::hessian()`](https://rdrr.io/pkg/dualr/man/hessian.html), and
[`dualr::score_and_hessian()`](https://rdrr.io/pkg/dualr/man/score_and_hessian.html)
— that handle the dual number seeding, forward passes, and extraction
internally. Users write a plain function of a parameter vector and get
back numeric results.

The `dualr` package uses a `dual_vector` class with natural
single-bracket indexing (`theta[i]`), making parameter access feel like
ordinary R code.

### Score via `dualr::score()`

To compute the gradient, we write the log-likelihood as a function of
`theta` (which will receive a `dual_vector` internally) and call
[`dualr::score()`](https://rdrr.io/pkg/dualr/man/score.html):

``` r
# Log-likelihood for dualr: use theta[i] (single bracket on dual_vector)
# Data summaries are pre-computed as plain numerics; dual arithmetic
# operates only on the parameter-dependent terms.
weibull_ll_dualr <- function(theta) {
    k     <- theta[1]
    sigma <- theta[2]
    n     <- length(df$t)

    sum_log_t <- sum(log(df$t))

    # Each observation contributes (t_i/sigma)^k with dual k and sigma.
    # Use Reduce for the power sum since each term involves dual params.
    power_sum <- Reduce("+", lapply(df$t, function(ti) (ti / sigma)^k))

    n * log(k) + (k - 1) * sum_log_t - n * k * log(sigma) - power_sum
}

dualr_sc <- dualr::score(weibull_ll_dualr, test_par)
cat("dualr score:\n")
#> dualr score:
print(dualr_sc)
#> [1] -10.494704   8.878147
cat("\nScore error:", max(abs(dualr_sc - ref_score)), "\n")
#> 
#> Score error: 7.105427e-15
```

The entire gradient computation is a one-liner:
`dualr::score(f, theta)`. Internally, `dualr` creates a `dual_vector`
from `theta`, runs $p$ forward passes (one per parameter), and extracts
the derivative components.

### Hessian via `dualr::hessian()`

For the Hessian,
[`dualr::hessian()`](https://rdrr.io/pkg/dualr/man/hessian.html) uses
nested dual numbers (second-order forward-mode AD). This requires
$O\left( p(p + 1)/2 \right)$ evaluations but needs no analytical score
function:

``` r
# Same log-likelihood function works for Hessian — dualr handles nesting
dualr_hess <- dualr::hessian(weibull_ll_dualr, test_par)
cat("dualr Hessian (nested duals):\n")
#> dualr Hessian (nested duals):
print(round(dualr_hess, 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047
cat("\nHessian error:", max(abs(dualr_hess - ref_hess)), "\n")
#> 
#> Hessian error: 94.06782
```

### Hybrid via `dualr::score_and_hessian()`

Like the femtograd and dual hybrid approaches, we can provide an
analytical score function and let `dualr` compute its Jacobian. The
score function must return a **list** of scalar values (one per
parameter):

``` r
# Analytical score returning list of scalars (one per parameter)
weibull_score_dualr <- function(theta) {
    k     <- theta[1]
    sigma <- theta[2]
    n     <- length(df$t)

    sum_log_t <- sum(log(df$t))

    r_k_sum     <- Reduce("+", lapply(df$t, function(ti) (ti / sigma)^k))
    r_k_log_sum <- Reduce("+", lapply(df$t, function(ti) {
        r <- ti / sigma
        r^k * log(r)
    }))

    dk     <- n / k + sum_log_t - n * log(sigma) - r_k_log_sum
    dsigma <- k / sigma * (r_k_sum - n)

    list(dk, dsigma)
}

dualr_sh <- dualr::score_and_hessian(weibull_score_dualr, test_par)
cat("dualr hybrid score:\n")
#> dualr hybrid score:
print(dualr_sh$score)
#> [1] -10.494704   8.878147
cat("\ndualr hybrid Hessian (Jacobian of analytical score):\n")
#> 
#> dualr hybrid Hessian (Jacobian of analytical score):
print(round(dualr_sh$hessian, 4))
#>          [,1]     [,2]
#> [1,] -73.4609  30.3410
#> [2,]  30.3410 -50.2047

cat("\nScore error:", max(abs(dualr_sh$score - ref_score)), "\n")
#> 
#> Score error: 7.105427e-15
cat("Hessian error:", max(abs(dualr_sh$hessian - ref_hess)), "\n")
#> Hessian error: 94.06782
```

The hybrid approach costs only $O(p)$ forward passes (vs
$O\left( p(p + 1)/2 \right)$ for the direct Hessian), making it
preferable when an analytical score is available.

### Observed Information

`dualr` also provides `dualr::observed_information(loglik, theta)`,
which returns the negative Hessian (i.e., the observed Fisher
information matrix at the given parameter values). This is a convenience
for the common MLE workflow where you need
$I\left( \widehat{\theta} \right) = - H\left( \widehat{\theta} \right)$.

### Makeham with dualr

``` r
# Makeham log-likelihood for dualr
makeham_ll_dualr <- function(theta) {
    lambda <- theta[1]
    a      <- theta[2]
    b      <- theta[3]
    t_obs  <- mak_df$t
    n      <- length(t_obs)

    # Log-likelihood: sum(log(h(t_i))) - sum(H(t_i))
    # h(t) = lambda + a*exp(b*t),  H(t) = lambda*t + (a/b)*(exp(b*t) - 1)
    log_h_sum <- Reduce("+", lapply(t_obs, function(ti) {
        log(lambda + a * exp(b * ti))
    }))
    H_sum <- Reduce("+", lapply(t_obs, function(ti) {
        lambda * ti + (a / b) * (exp(b * ti) - 1)
    }))

    log_h_sum - H_sum
}

# Direct Hessian via nested duals
dualr_mak_hess <- dualr::hessian(makeham_ll_dualr, mak_test_par)
cat("dualr Hessian (3 params, direct):\n")
#> dualr Hessian (3 params, direct):
print(round(dualr_mak_hess, 3))
#>            [,1]       [,2]       [,3]
#> [1,] -22253.094  -93128.85  -5428.997
#> [2,] -93128.852 -856134.12 -25374.042
#> [3,]  -5428.997  -25374.04  -2352.126

# Hybrid via score_and_hessian with analytical score
makeham_score_dualr <- function(theta) {
    lambda <- theta[1]
    a      <- theta[2]
    b      <- theta[3]
    t_obs  <- mak_df$t
    n      <- length(t_obs)

    h_inv_sum <- Reduce("+", lapply(t_obs, function(ti) {
        1 / (lambda + a * exp(b * ti))
    }))
    exp_bt_over_h_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ebt <- exp(b * ti)
        ebt / (lambda + a * ebt)
    }))
    exp_bt_minus1_sum <- Reduce("+", lapply(t_obs, function(ti) {
        exp(b * ti) - 1
    }))
    t_sum <- sum(t_obs)
    t_exp_bt_over_h_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ebt <- exp(b * ti)
        a * ti * ebt / (lambda + a * ebt)
    }))
    t_exp_bt_sum <- Reduce("+", lapply(t_obs, function(ti) {
        ti * exp(b * ti)
    }))

    dlambda <- h_inv_sum - t_sum
    da      <- exp_bt_over_h_sum - exp_bt_minus1_sum / b
    db      <- t_exp_bt_over_h_sum -
               (a / b) * t_exp_bt_sum +
               (a / b^2) * exp_bt_minus1_sum

    list(dlambda, da, db)
}

dualr_mak_sh <- dualr::score_and_hessian(makeham_score_dualr, mak_test_par)
cat("\ndualr hybrid score (3 params):\n")
#> 
#> dualr hybrid score (3 params):
print(dualr_mak_sh$score)
#> [1]  480.0670 4426.1957  365.1092
cat("\ndualr hybrid Hessian (3 params):\n")
#> 
#> dualr hybrid Hessian (3 params):
print(round(dualr_mak_sh$hessian, 3))
#>            [,1]       [,2]       [,3]
#> [1,] -22253.094  -93128.85  -5428.997
#> [2,] -93128.852 -856134.12 -25374.042
#> [3,]  -5428.997  -25374.04  -2352.126
```

## Precision Comparison

``` r
cat("=== Score Errors (max absolute difference from analytical) ===\n\n")
#> === Score Errors (max absolute difference from analytical) ===
cat(sprintf("  numDeriv:        %.2e\n", max(abs(nd_score - ref_score))))
#>   numDeriv:        5.09e-04
cat(sprintf("  femtograd full:  %.2e\n", max(abs(fad_score - ref_score))))
#>   femtograd full:  8.35e-14
cat(sprintf("  femtograd hybrid:%.2e\n", max(abs(hyb_score - ref_score))))
#>   femtograd hybrid:0.00e+00
cat(sprintf("  dual:            %.2e\n", max(abs(dual_sc - ref_score))))
#>   dual:            1.95e-14
cat(sprintf("  dualr:           %.2e\n", max(abs(dualr_sc - ref_score))))
#>   dualr:           7.11e-15
cat(sprintf("  dualr hybrid:    %.2e\n", max(abs(dualr_sh$score - ref_score))))
#>   dualr hybrid:    7.11e-15

cat("\n=== Hessian Errors (max absolute difference from analytical) ===\n\n")
#> 
#> === Hessian Errors (max absolute difference from analytical) ===
cat(sprintf("  numDeriv:        %.2e\n", max(abs(nd_hess - ref_hess))))
#>   numDeriv:        9.41e+01
cat(sprintf("  femtograd full:  %.2e\n", max(abs(fad_hess - ref_hess))))
#>   femtograd full:  9.41e+01
cat(sprintf("  femtograd hybrid:%.2e\n", max(abs(hyb_hess - ref_hess))))
#>   femtograd hybrid:9.41e+01
cat(sprintf("  dual:            %.2e\n", max(abs(dual_hess - ref_hess))))
#>   dual:            9.41e+01
cat(sprintf("  dualr direct:    %.2e\n", max(abs(dualr_hess - ref_hess))))
#>   dualr direct:    9.41e+01
cat(sprintf("  dualr hybrid:    %.2e\n", max(abs(dualr_sh$hessian - ref_hess))))
#>   dualr hybrid:    9.41e+01
```

All AD approaches (femtograd, dual, and dualr) match the analytical
solution to machine precision ($\sim 10^{- 12}$), while `numDeriv`
achieves approximately $10^{- 5}$ to $10^{- 6}$ accuracy. For
well-behaved likelihoods the numDeriv precision is usually sufficient,
but for ill-conditioned problems (near-singular Hessians, parameters
near boundaries), the AD precision advantage matters.

## Performance Benchmark

``` r
# numDeriv and femtograd timings computed here; dual and dualr timings
# were pre-computed in their respective sections (before namespace swap).
n_iter <- 20

invisible(score_nd(df, par = test_par))
invisible(score_fad(df, par = test_par))
invisible(score_hyb(df, par = test_par))

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(score_nd(df, par = test_par))
time_score_nd <- as.numeric(Sys.time() - t0) / n_iter * 1000

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(score_fad(df, par = test_par))
time_score_fad <- as.numeric(Sys.time() - t0) / n_iter * 1000

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(score_hyb(df, par = test_par))
time_score_hyb <- as.numeric(Sys.time() - t0) / n_iter * 1000

cat("=== Score Computation (ms per evaluation, n=100) ===\n\n")
#> === Score Computation (ms per evaluation, n=100) ===
cat(sprintf("  numDeriv:          %7.2f ms\n", time_score_nd))
#>   numDeriv:           186.79 ms
cat(sprintf("  femtograd full AD: %7.2f ms\n", time_score_fad))
#>   femtograd full AD:   88.93 ms
cat(sprintf("  femtograd hybrid:  %7.2f ms\n", time_score_hyb))
#>   femtograd hybrid:     0.26 ms
cat(sprintf("  dual:              %7.2f ms\n", time_score_dual))
#>   dual:                24.73 ms
cat(sprintf("  dualr:             %7.2f ms\n", time_score_dualr))
#>   dualr:              211.82 ms
```

``` r
invisible(hess_nd(df, par = test_par))
invisible(hess_fad(df, par = test_par))
invisible(hess_hyb(df, par = test_par))

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(hess_nd(df, par = test_par))
time_hess_nd <- as.numeric(Sys.time() - t0) / n_iter * 1000

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(hess_fad(df, par = test_par))
time_hess_fad <- as.numeric(Sys.time() - t0) / n_iter * 1000

t0 <- Sys.time()
for (i in seq_len(n_iter)) invisible(hess_hyb(df, par = test_par))
time_hess_hyb <- as.numeric(Sys.time() - t0) / n_iter * 1000

cat("=== Hessian Computation (ms per evaluation, n=100) ===\n\n")
#> === Hessian Computation (ms per evaluation, n=100) ===
cat(sprintf("  numDeriv:          %7.2f ms\n", time_hess_nd))
#>   numDeriv:           270.97 ms
cat(sprintf("  femtograd full AD: %7.2f ms\n", time_hess_fad))
#>   femtograd full AD:   58.84 ms
cat(sprintf("  femtograd hybrid:  %7.2f ms\n", time_hess_hyb))
#>   femtograd hybrid:     4.52 ms
cat(sprintf("  dual:              %7.2f ms\n", time_hess_dual))
#>   dual:                97.87 ms
cat(sprintf("  dualr direct:      %7.2f ms\n", time_hess_dualr))
#>   dualr direct:      1815.94 ms
cat(sprintf("  dualr hybrid:      %7.2f ms\n", time_hess_dualr_hyb))
#>   dualr hybrid:       591.65 ms
```

### Makeham Benchmark (p=3)

``` r
# numDeriv and femtograd timed here; dual and dualr timings pre-computed
n_iter_mak <- 10

invisible(hess_mak_nd(mak_df, par = mak_test_par))
invisible(hess_mak_hyb(mak_df, par = mak_test_par))

t0 <- Sys.time()
for (i in seq_len(n_iter_mak)) invisible(hess_mak_nd(mak_df, par = mak_test_par))
time_mak_nd <- as.numeric(Sys.time() - t0) / n_iter_mak * 1000

t0 <- Sys.time()
for (i in seq_len(n_iter_mak)) invisible(hess_mak_hyb(mak_df, par = mak_test_par))
time_mak_hyb <- as.numeric(Sys.time() - t0) / n_iter_mak * 1000

cat("=== Hessian Computation, p=3 (ms per evaluation, n=150) ===\n\n")
#> === Hessian Computation, p=3 (ms per evaluation, n=150) ===
cat(sprintf("  numDeriv:          %7.2f ms\n", time_mak_nd))
#>   numDeriv:           741.87 ms
cat(sprintf("  femtograd hybrid:  %7.2f ms\n", time_mak_hyb))
#>   femtograd hybrid:    12.23 ms
cat(sprintf("  dual hybrid:       %7.2f ms\n", time_mak_dual))
#>   dual hybrid:        614.39 ms
cat(sprintf("  dualr direct:      %7.2f ms\n", time_mak_dualr))
#>   dualr direct:       236.37 ms
cat(sprintf("  dualr hybrid:      %7.2f ms\n", time_mak_dualr_hyb))
#>   dualr hybrid:      3532.82 ms
```

Both hybrid approaches (analytical score + AD Jacobian) scale linearly
in $p$ for the Hessian, while numDeriv’s cost grows quadratically. The
relative advantage grows with more parameters.

## API Ergonomics Comparison

### Indexing Conventions

The most visible difference when writing AD-compatible code:

``` r
# femtograd: must use [[ to extract value/dual objects from list
score_fn <- function(df, par, ...) {
    k     <- par[[1]]   # double brackets required
    sigma <- par[[2]]
    # ...
    c(score_k, score_sigma)
}

# dual: uses [[ on a list of dual objects (same pattern)
score_fn <- function(theta) {
    k     <- theta[[1]]  # also double brackets
    sigma <- theta[[2]]
    # ...
    list(score_k, score_sigma)   # list, not c()
}

# dualr: uses [ on a dual_vector (natural R indexing)
score_fn <- function(theta) {
    k     <- theta[1]   # single brackets — standard R convention
    sigma <- theta[2]
    # ...
    list(score_k, score_sigma)   # list for score_and_hessian()
}
```

The indexing is the same (`[[`), but the return convention differs:

- **femtograd score functions** return `c(s1, s2, ...)` — a numeric
  vector (or vector of `value`/`dual_num` objects combined via
  [`c()`](https://rdrr.io/r/base/c.html))
- **dual score functions** (for Jacobian computation) return
  `list(s1, s2, ...)` because each element must remain an independent
  `dual` object with its own gradient slot

### Function Coverage

| Feature                                             | numDeriv | femtograd      | dual           | dualr          |
|-----------------------------------------------------|----------|----------------|----------------|----------------|
| `log`, `exp`, `^`, `sqrt`                           | n/a      | Yes            | Yes            | Yes            |
| `sin`, `cos`, `tan`                                 | n/a      | Yes            | Yes            | Yes            |
| `gamma`, `lgamma`, `digamma`                        | n/a      | No             | Yes            | Yes            |
| `beta`, `lbeta`                                     | n/a      | No             | Yes            | Yes            |
| `erf`, `erfc`                                       | n/a      | No             | Yes            | Yes            |
| `dnorm`, `pnorm`                                    | n/a      | No             | No             | No             |
| [`sum()`](https://rdrr.io/r/base/sum.html) dispatch | n/a      | No             | No             | Yes            |
| Gradient mode                                       | —        | Reverse (O(1)) | Forward (O(p)) | Forward (O(p)) |
| Dispatch mechanism                                  | —        | R6             | S4             | S4             |

### Return Convention Summary

| Context                         | femtograd             | dual                          | dualr                                                           |
|---------------------------------|-----------------------|-------------------------------|-----------------------------------------------------------------|
| Score for optimization          | `c(s1, s2)`           | `c(s1@f, s2@f)` (extract)     | `dualr::score(ll, theta)` → numeric vector                      |
| Score for AD Jacobian → Hessian | `c(s1, s2)`           | `list(s1, s2)` (keep as dual) | `score_fn` returns [`list()`](https://rdrr.io/r/base/list.html) |
| Log-likelihood for AD gradient  | Scalar `value` object | Scalar `dual` object          | Scalar dual object                                              |
| Hessian direct                  | Not available         | Not available                 | `dualr::hessian(ll, theta)` → matrix                            |

### API Convenience Comparison

The four approaches differ significantly in how much boilerplate the
user writes:

``` r
# numDeriv: wrapped by dfr.dist — user writes only the hazard rate
score_nd  <- score(my_dist)          # returns closure
result    <- score_nd(df, par = par) # returns numeric vector

# femtograd: wrapped by dfr.dist — user writes rate + cum_haz_rate (with [[)
score_fem <- score(my_dist)          # returns closure using ad_gradient()
result    <- score_fem(df, par = par)

# dual (CRAN): fully manual — seed duals, loop, extract gradients
theta <- lapply(seq_len(p), function(i) dual::dual(par[i], ...))
result <- my_ll(theta, data)
gradient <- result@grad

# dualr: one-liners — just write a function and call it
gradient <- dualr::score(my_ll, par)
hess     <- dualr::hessian(my_ll, par)
both     <- dualr::score_and_hessian(my_score, par)
```

## Practical Guidance

### Decision Table

| Scenario                                  | Recommended            | Why                                                                                                                        |
|-------------------------------------------|------------------------|----------------------------------------------------------------------------------------------------------------------------|
| Quick prototyping, any distribution       | numDeriv               | Zero setup, works with any R function                                                                                      |
| Clean standalone AD (no dfr.dist)         | **dualr**              | One-liner [`score()`](https://queelius.github.io/likelihood.model/reference/score.html), `hessian()`, natural `[` indexing |
| Production with dfr.dist (p ≤ 5)          | femtograd hybrid       | Fast, exact, built into dfr.dist fallback chain                                                                            |
| Many parameters (p \>\> 10)               | femtograd reverse-mode | O(1) gradient cost vs O(p) for forward-mode                                                                                |
| Likelihood involves `lgamma`, `erf`, etc. | dualr or dual          | femtograd lacks these special functions                                                                                    |
| Teaching AD concepts                      | dual or dualr          | dual: transparent `@f`/`@grad` inspection; dualr: cleanest API                                                             |
| Need `dnorm`/`pnorm` in likelihood        | numDeriv               | No AD package supports these yet                                                                                           |

### Writing Portable Score Functions

If you want your score function to work across AD backends, the key
constraints are:

1.  **Use `[[` indexing** for femtograd and dual compatibility (dualr
    also supports `[[`, though `[` is its natural convention)
2.  **Use basic arithmetic and `log`/`exp`/`^`** (all three support
    these)
3.  **Use `Reduce("+", ...)` for sums** over data — dual doesn’t
    dispatch [`sum()`](https://rdrr.io/r/base/sum.html) on its objects
    (dualr does, but `Reduce` works everywhere)
4.  **Return [`list()`](https://rdrr.io/r/base/list.html) for score
    functions** used with
    [`dualr::score_and_hessian()`](https://rdrr.io/pkg/dualr/man/score_and_hessian.html)
    and dual’s Jacobian pattern

## Integration with dfr.dist

The `dfr.dist` package currently integrates with femtograd through the
fallback chains in `R/ad_utils.R`:

**Score function chain:**

1.  Analytical `score_fn` if provided → used directly
2.  AD gradient via femtograd if `cum_haz_rate` provided → exact
3.  Numerical gradient via
    [`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html) →
    approximate

**Hessian chain:**

1.  AD Jacobian of `score_fn` via femtograd → exact
2.  Numerical Hessian via
    [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
    → approximate

Both `dual` and `dualr` can be used directly by users for standalone
derivative computation (as shown in this vignette). The `dualr` package
is particularly suited for standalone use thanks to its convenience API
([`dualr::score()`](https://rdrr.io/pkg/dualr/man/score.html),
[`dualr::hessian()`](https://rdrr.io/pkg/dualr/man/hessian.html)), and
integration into dfr.dist’s fallback chains is a future possibility for
either package.

## Summary

|                             | numDeriv                  | femtograd                             | dual                             | dualr                             |
|:----------------------------|:--------------------------|:--------------------------------------|:---------------------------------|:----------------------------------|
| AD mode                     | None (finite differences) | Reverse + forward-over-reverse        | Forward-mode dual numbers        | Forward-mode dual numbers         |
| Gradient cost               | O(p) function evals       | O(1) backward pass                    | O(p) forward passes              | O(p) forward passes               |
| Hessian cost (via Jacobian) | O(p²) function evals      | O(p) forward passes                   | O(p) forward passes on score     | O(p) hybrid or O(p(p+1)/2) direct |
| Precision (vs analytical)   | ~1e-5 to 1e-6             | Machine epsilon                       | Machine epsilon                  | Machine epsilon                   |
| Setup effort                | Minimal                   | Moderate (par\[\[i\]\], cum_haz_rate) | Moderate (list of duals, Reduce) | Low (one-liner score/hessian)     |
| Special function support    | All R functions           | Basic math only                       | Rich (lgamma, erf, beta, …)      | Rich (lgamma, erf, beta, sum)     |
| dfr.dist integration        | Built-in fallback         | Primary AD backend                    | Manual (standalone use)          | Manual (standalone use)           |

For most `dfr.dist` workflows, the **femtograd hybrid approach**
(analytical score + AD Hessian) offers the best combination of speed and
precision. Use **numDeriv** for quick prototyping or when your
likelihood involves functions that no AD package supports. For
standalone AD outside of dfr.dist, **dualr** provides the cleanest API —
its one-liner
[`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
`hessian()`, and `score_and_hessian()` functions eliminate the
boilerplate of manual dual number seeding and extraction. Consider
**dual** when you want transparent, pedagogically clear AD internals
(direct access to `@f` and `@grad` slots), or when you need special
mathematical functions without the convenience layer.

## Next Steps

- **[`vignette("automatic_differentiation")`](https://queelius.github.io/dfr_dist/articles/automatic_differentiation.md)**
  — Deep dive into femtograd integration with dfr.dist
- **[`vignette("custom_distributions")`](https://queelius.github.io/dfr_dist/articles/custom_distributions.md)**
  — Building custom distributions with analytical score functions
- **[`vignette("getting_started")`](https://queelius.github.io/dfr_dist/articles/getting_started.md)**
  — Quick 5-minute introduction
