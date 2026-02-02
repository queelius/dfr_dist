#' Constructor for `dfr_dist` objects
#'
#' We assume that the hazard rate is a function of time and any other
#' predictors. We also assume that integrate(rate(t), 0, Inf) = infinity
#' and that the support is (0, Inf).
#'
#' @param rate A function that computes the hazard rate at time `t`.
#' @param par The parameters of the distribution. Defaults to `NULL`,
#'            which means that the parameters are unknown.
#' @param eps The epsilon update for numerical integration. Defaults to 0.01.
#' @param ob_col The column name for observation times in data frames.
#'               Defaults to "t".
#' @param delta_col The column name for event indicators in data frames.
#'                  Uses standard survival analysis convention: 1 = event
#'                  observed (exact), 0 = right-censored. Defaults to "delta".
#' @param cum_haz_rate Optional analytical cumulative hazard function H(t, par).
#'                     If provided, used for faster exact cumulative hazard
#'                     computation instead of numerical integration.
#'                     Should return the integral of rate from 0 to t.
#' @param score_fn Optional score function (gradient of log-likelihood).
#'                 Signature: score_fn(df, par, ...) returning a numeric vector.
#'                 If NULL, falls back to numerical gradient via numDeriv::grad.
#' @param hess_fn Optional Hessian function (second derivatives of log-likelihood).
#'                Signature: hess_fn(df, par, ...) returning a matrix.
#'                If NULL, falls back to numerical Hessian via numDeriv::hessian.
#' @return A `dfr_dist` object that inherits from `likelihood_model`.
#' @export
dfr_dist <- function(rate, par = NULL, eps = 0.01,
                     ob_col = "t", delta_col = "delta",
                     cum_haz_rate = NULL, score_fn = NULL,
                     hess_fn = NULL) {
    structure(
        list(rate = rate,
             par = par,
             eps = eps,
             ob_col = ob_col,
             delta_col = delta_col,
             cum_haz_rate = cum_haz_rate,
             score_fn = score_fn,
             hess_fn = hess_fn),
    class = c("dfr_dist", "likelihood_model", "univariate_dist", "dist"))
}

#' Function for determining whether an object is a `dfr_dist` object.
#' 
#' @param x The object to test.
#' @return A logical value indicating whether `x` is a `dfr_dist` 
#' object.
#' @export
is_dfr_dist <- function(x) {
    inherits(x, "dfr_dist")
}

#' Method for obtaining the hazard function of a `dfr_dist` object.
#'
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass into the `rate` function.
#' @return A function that computes the hazard function of the distribution.
#' It accepts `t`, the time at which to compute the hazard function, and
#' `par`, the parameters of the distribution. If `par` is `NULL`, then
#' the parameters of the `dfr_dist` object `x` are used. It also
#' accepts a `log` argument that determines whether to compute the log of
#' the hazard function. Finally, it passes any additional arguments to the
#' `rate` function of the `dfr_dist` object `x`.
#' @importFrom algebraic.dist hazard params
#' @export
hazard.dfr_dist <- function(x, ...) {
    function(t, par = NULL, ...) {
        par <- get_params(par, x$par)
        x$rate(t, par, ...)
    }
}

#' Method for obtaining the quantile (inverse cdf) of an object.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass into `cdf` constructor.
#' @return A function that computes the quantile of the distribution.
#' It accepts `p`, the probability at which to compute the quantile, 
#' `par`, the parameters of the distribution, and `...`, any additional
#' arguments to pass into the constructed cdf.
#' @importFrom stats uniroot
#' @importFrom algebraic.dist params cdf inv_cdf
#' @export
inv_cdf.dfr_dist <- function(x, ...) {

    cdf_fn <- cdf(x, ...)
    function(p, par = NULL, ...) {
        par <- get_params(par, x$par)
        uniroot_args <- list(
            f = function(t) {
                cdf_fn(t, par, ...) - p
            },
            interval = c(0, 1e3),
            extendInt = "upX")
        do.call(uniroot, uniroot_args)$root
    }
}

#' Method for obtaining the parameters of a `dfr_dist` object.
#'
#' @param x The object to obtain the parameters of.
#' @param ... Additional arguments (unused).
#' @return The parameters of the distribution.
#' @importFrom algebraic.dist params
#' @export
params.dfr_dist <- function(x, ...) {
    x$par
}

#' Sampling function for `dfr_dist` objects.
#' 
#' Since S(t,par) = exp(-cum_hz(t,par)), we can sample from the
#' distribution by letting t = 0 (or some other positive number if
#' we want to condition on T > t_min), sampling from an exponential
#' distribution with `lambda = rate(t, par)`, and then rejecting
#' the sample if `runif(1) > S(t, par)`. If accepted, add that
#' observation to the sample, otherwise reject it, let `t = t + eps`
#' where `eps` is some small number, and repeat. We continue this
#' process until we have `n` observations for the sample.
#'
#' @param x The object to obtain the sampler of.
#' @param ... Additional arguments to pass into the survival function
#' @return A function that samples from the distribution. It accepts
#' `n`, the number of samples to take, `t` is the time at which to start
#' sampling, `par` are the parameters of the distribution, and `eps` is
#' the  update for numerical integration. Finally, we pass additional
#' arguments `...` into the hazard function.
#' @importFrom algebraic.dist params surv sampler
#' @importFrom stats runif
#' @export
sampler.dfr_dist <- function(x, ...) {
    # Use inverse CDF sampling for reliability
    Q <- inv_cdf(x, ...)
    function(n, par = NULL, ...) {
        par <- get_params(par, x$par)
        u <- runif(n)
        sapply(u, function(p) Q(p, par, ...))
    }
}

#' Method for obtaining the cdf of a `dfr_dist` object.
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass into the `cum_haz`
#' constructor.
#' @return A function that computes the cdf of the distribution.
#' It accepts `t`, the time at which to compute the cdf, `par`,
#' the parameters of the distribution, `log.p` argument that
#' determines whether to compute the log of the cdf, `lower.limit`,
#' whether to compute the lower limit (F(t)) or upper limit
#' (S(t) = 1-F(t)). Finally, it passes any additional arguments `...`
#' to the `rate` function of the `dfr_dist` object `x`.
#' 
#' @importFrom algebraic.dist params cdf
#' @export
cdf.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    function(t, par = NULL, log.p = FALSE, lower.limit = TRUE, ...) {
        par <- get_params(par, x$par)
        haz <- H(t, par, ...)
        if (lower.limit) {
            p <- 1 - exp(-haz)
            ifelse(log.p, log(p), p)
        } else {
            ifelse(log.p, -haz, exp(-haz))
        }
    }
}

#' Method for obtaining the density (pdf) of a `dfr_dist` object.
#'
#' @param x The object to obtain the density of.
#' @param ... Additional arguments to pass.
#' @return A function that computes the density of the distribution.
#' It accepts `t`, the time at which to compute the density, `par` is
#' the parameters of the distribution, and `log` determines whether to
#' compute the log of the density. Finally, it passes any additional
#' arguments `...` to the `rate` function of the `dfr_dist` object `x`.
#' @importFrom stats density
#' @export
density.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    inner <- function(t, par, log, ...) {
        if (log) {
            log(x$rate(t, par, ...)) - H(t, par, ...)
        } else {
            x$rate(t, par, ...) * exp(-H(t, par, ...))
        }
    }
    function(t, par = NULL, log = FALSE, ...) {
        par <- get_params(par, x$par)
        sapply(t, function(ti) inner(ti, par, log, ...))
    }
}


#' Method for retrieving the support of an object `x`.
#'
#' @param x The object to obtain the support of.
#' @param ... Additional arguments to pass.
#' @return A support object for `x`, an interval (0,Inf).
#' @importFrom algebraic.dist interval sup
#' @export
sup.dfr_dist <- function(x, ...) {
    interval$new(0, Inf, FALSE, FALSE)
}

#' Method for obtaining the cumulative hazard function of a `dfr_dist` object.
#' @param x The object to obtain the cumulative hazard function of.
#' @param ... Additional arguments to pass into the `integrate` function
#'   (only used when no analytical cum_haz_rate is provided).
#' @return A function that computes the cumulative hazard H(t) of the distribution.
#' It accepts `t`, the time at which to compute the cumulative hazard, and
#' `par`, the parameters of the distribution. If `par` is `NULL`, then the
#' parameters of the `dfr_dist` object `x` are used. Finally, it passes any
#' additional arguments `...` to the rate function.
#' @details
#' If the `dfr_dist` object has an analytical `cum_haz_rate` function, that is
#' used directly for fast, exact computation. Otherwise, numerical integration
#' of the hazard function is performed.
#' @importFrom stats integrate
#' @importFrom utils modifyList
#' @export
cum_haz.dfr_dist <- function(x, ...) {
    if (!is.null(x$cum_haz_rate)) {
        return(function(t, par = NULL, ...) {
            par <- get_params(par, x$par)
            x$cum_haz_rate(t, par, ...)
        })
    }

    integrator_defaults <- list(
        lower = 0, subdivisions = 1000L, abs.tol = 1e-3)
    integrator <- modifyList(integrator_defaults, list(...))

    function(t, par = NULL, ...) {
        par <- get_params(par, x$par)
        res <- do.call(integrate,
            modifyList(integrator, list(
                upper = t,
                f = function(u) x$rate(u, par, ...))))
        if (res$message != "OK") {
            warning(res$message)
        }
        if (res$abs.error > integrator$abs.tol) {
            warning("Absolute error in cumulative hazard is greater than tolerance")
        }
        res$value
    }
}

#' Method for obtaining the survival function of a `dfr_dist` object.
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass into the `cum_haz`
#' constructor.
#' @return A function that computes the survival function of the
#' distribution.
#' It accepts `t`, the time at which to compute the survival, `par`,
#' the parameters of the distribution, `log.p` argument that
#' determines whether to compute the log of the survival, and
#' it passes any additional arguments into the `rate` function of
#' the `dfr_dist` object `x`.
#' @importFrom algebraic.dist params surv
#' @export
surv.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    function(t = 0, par = NULL, log.p = FALSE, ...) {
        par <- get_params(par, x$par)
        haz <- H(t, par, ...)
        ifelse(log.p, -haz, exp(-haz))
    }
}

#' Print method for `dfr_dist` objects.
#' @param x The `dfr_dist` object to print.
#' @param ... Additional arguments (not used)
#' @export
print.dfr_dist <- function(x, ...) {
    cat("Dynamic failure rate (DFR) distribution with failure rate:\n")
    print(x$rate)
    cat("It has a survival function given by:\n")
    cat("    S(t|rate) = exp(-H(t,...))\n")
    cat("where H(t,...) is the cumulative hazard function.\n")
    invisible(x)
}

#' Log-likelihood method for `dfr_dist` objects
#'
#' Returns a function that computes the log-likelihood of the data given
#' the distribution parameters. The log-likelihood for survival data is:
#'
#' For exact observations (uncensored): log(f(t)) = log(h(t)) - H(t)
#' For right-censored observations: log(S(t)) = -H(t)
#' For left-censored observations: log(F(t)) = log(1 - exp(-H(t)))
#'
#' where h(t) is the hazard function, H(t) is the cumulative hazard,
#' f(t) = h(t)*S(t) is the pdf, and S(t) = exp(-H(t)) is the survival function.
#'
#' @param model The `dfr_dist` object
#' @param ... Additional arguments to pass to the hazard and cumulative hazard
#' @return A function that computes the log-likelihood. It accepts:
#'         - `df`: A data frame with observation times and censoring indicators
#'           (delta: 1 = exact, 0 = right-censored, -1 = left-censored)
#'         - `par`: The parameters of the distribution
#'         - `...`: Additional arguments passed to internal functions
#' @importFrom likelihood.model loglik
#' @export
loglik.dfr_dist <- function(model, ...) {
    H <- cum_haz(model, ...)

    function(df, par = NULL, ...) {
        par <- get_params(par, model$par)

        t <- df[[model$ob_col]]

        if (model$delta_col %in% names(df)) {
            delta <- df[[model$delta_col]]
        } else {
            delta <- rep(1, length(t))
        }

        ll <- 0

        # Exact observations (delta = 1): log(h(t)) - H(t)
        exact_idx <- which(delta == 1)
        if (length(exact_idx) > 0) {
            t_exact <- t[exact_idx]
            h_exact <- sapply(t_exact, function(ti) model$rate(ti, par, ...))
            H_exact <- sapply(t_exact, function(ti) H(ti, par, ...))
            contrib <- sum(log(h_exact) - H_exact)
            if (is.nan(contrib)) return(-Inf)
            ll <- ll + contrib
        }

        # Right-censored observations (delta = 0): -H(t)
        right_idx <- which(delta == 0)
        if (length(right_idx) > 0) {
            t_right <- t[right_idx]
            H_right <- sapply(t_right, function(ti) H(ti, par, ...))
            ll <- ll - sum(H_right)
        }

        # Left-censored observations (delta = -1): log(1 - exp(-H(t)))
        left_idx <- which(delta == -1)
        if (length(left_idx) > 0) {
            t_left <- t[left_idx]
            H_left <- sapply(t_left, function(ti) H(ti, par, ...))
            ll <- ll + sum(log1p(-exp(-H_left)))
        }

        ll
    }
}

#' Score function (gradient of log-likelihood) for dfr_dist
#'
#' Returns a function that computes the score (gradient of log-likelihood)
#' with respect to parameters. Uses user-provided score function if available,
#' otherwise falls back to numerical differentiation via numDeriv::grad.
#'
#' @param model A dfr_dist object
#' @param ... Additional arguments passed to loglik
#' @return A function that computes the score vector
#' @importFrom likelihood.model score
#' @export
score.dfr_dist <- function(model, ...) {
    ll_fn <- loglik(model, ...)

    function(df, par = NULL, ...) {
        par <- get_params(par, model$par)
        if (!is.null(model$score_fn)) {
            return(model$score_fn(df, par, ...))
        }
        numDeriv::grad(function(p) ll_fn(df, par = p, ...), par)
    }
}

#' Hessian of log-likelihood for dfr_dist
#'
#' Returns a function that computes the Hessian matrix of the log-likelihood.
#' Uses user-provided Hessian function if available, otherwise falls back to
#' numerical differentiation via numDeriv::hessian.
#'
#' @param model A dfr_dist object
#' @param ... Additional arguments passed to loglik
#' @return A function that computes the Hessian matrix
#' @importFrom likelihood.model hess_loglik
#' @export
hess_loglik.dfr_dist <- function(model, ...) {
    ll_fn <- loglik(model, ...)

    function(df, par = NULL, ...) {
        par <- get_params(par, model$par)
        if (!is.null(model$hess_fn)) {
            return(model$hess_fn(df, par, ...))
        }
        numDeriv::hessian(function(p) ll_fn(df, par = p, ...), par)
    }
}

#' Retrieve the assumptions a DFR distribution makes about the data
#'
#' Returns a list of assumptions that the dynamic failure rate distribution
#' model makes about the underlying data and process.
#'
#' @param model A `dfr_dist` object
#' @param ... Additional arguments (ignored)
#' @return A character vector of model assumptions
#' @importFrom likelihood.model assumptions
#' @export
assumptions.dfr_dist <- function(model, ...) {
    c(
        "Non-negative hazard: h(t) >= 0 for all t > 0",
        "Cumulative hazard diverges: lim(t->Inf) H(t) = Inf",
        "Support is positive reals: t in (0, Inf)",
        "Observations are independent",
        "Censoring indicator convention: 1=exact, 0=right-censored, -1=left-censored",
        "Non-informative censoring: censoring mechanism independent of failure time"
    )
}

#' MLE solver for dfr_dist objects
#'
#' Returns a solver function that fits a dfr_dist model to survival data
#' using maximum likelihood estimation. The solver uses gradient-based
#' optimization (BFGS by default) with analytical or numerical gradients.
#'
#' @param object A `dfr_dist` object
#' @param ... Additional arguments passed to the log-likelihood, score, and
#'            Hessian constructors (e.g., integration parameters for cum_haz)
#' @return A solver function that accepts:
#'   - `df`: Data frame with observation times and censoring indicators
#'   - `par`: Initial parameter values (uses object's params if NULL)
#'   - `method`: Optimization method (default "BFGS")
#'   - `control`: Control parameters for optim()
#'   - `...`: Additional arguments passed to likelihood functions
#'
#' @details
#' The solver returns a `fisher_mle` object (from likelihood.model) containing:
#' - Parameter estimates
#' - Log-likelihood value at MLE
#' - Variance-covariance matrix (from Hessian)
#' - Convergence status
#'
#' Use methods like `coef()`, `vcov()`, `confint()`, `summary()` on the result.
#'
#' @examples
#' \donttest{
#' # Exponential distribution
#' exp_dist <- dfr_dist(
#'   rate = function(t, par, ...) rep(par[1], length(t)),
#'   par = c(lambda = 1)
#' )
#'
#' # Simulate data
#' set.seed(42)
#' df <- data.frame(t = rexp(100, rate = 2), delta = 1)
#'
#' # Fit model
#' solver <- fit(exp_dist)
#' result <- solver(df, par = c(1))
#' summary(result)
#' confint(result)
#' }
#'
#' @importFrom generics fit
#' @importFrom stats optim
#' @importFrom utils modifyList
#' @export
fit.dfr_dist <- function(object, ...) {
    ll <- loglik(object, ...)
    s <- score(object, ...)
    H <- hess_loglik(object, ...)

    function(df, par = NULL,
             method = c("BFGS", "Nelder-Mead", "L-BFGS-B", "CG", "SANN"),
             control = list(), ...) {
        if (is.null(par)) {
            par <- params(object)
            if (is.null(par)) {
                stop("Initial parameters required: specify 'par' argument ",
                     "or set parameters in dfr_dist()")
            }
        }

        control <- modifyList(list(fnscale = -1), control)
        method <- match.arg(method)

        sol <- optim(
            par = par,
            fn = function(p) ll(df, p, ...),
            gr = if (method == "SANN") NULL else function(p) s(df, p, ...),
            method = method,
            control = control
        )

        likelihood.model::fisher_mle(
            par = sol$par,
            loglik_val = sol$value,
            hessian = H(df, sol$par, ...),
            score_val = s(df, sol$par, ...),
            nobs = nrow(df),
            converged = (sol$convergence == 0),
            optim_result = sol
        )
    }
}
