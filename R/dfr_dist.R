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
#' @return A `dfr_dist` object that inherits from `likelihood_model`.
#' @export
dfr_dist <- function(rate, par = NULL, eps = 0.01,
                     ob_col = "t", delta_col = "delta") {
    structure(
        list(rate = rate,
             par = par,
             eps = eps,
             ob_col = ob_col,
             delta_col = delta_col),
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

    F <- cdf(x, ...)
    function(p, par = NULL, ...) {
        par <- get_params(par, x$par)
        # F(t) = p, so we want to find t such that F(t) = p.
        uniroot_args <- list(
            f = function(t) {
                F(t, par, ...) - p
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
    function (t, par = NULL, log.p = FALSE, lower.limit = TRUE, ...) {
        par <- get_params(par, x$par)
        haz <- H(t, par, ...)
        if (lower.limit) {
            p <- 1 - exp(-haz)
            return(ifelse(log.p, log(p), p))
        }
        else {            
            return(ifelse(log.p, -haz, exp(-haz)))
        }
    }
}

#' Method for obtaining the pdf of a `dfr_dist` object.
#' 
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @return A function that computes the pdf of the distribution.
#' It accepts `y`, the value at which to compute the pdf, `t`, the time
#' at which to compute the pdf, `par` is the parameters of the
#' distribution, and `log` determines whether to compute the log of
#' the pdf. Finally, it passes any additional arguments `...` to
#' the `rate` function of the `dfr_dist` object `x`.
#' @importFrom algebraic.dist pdf params
#' @export
pdf.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    # Inner function works on scalar t for integration compatibility
    inner <- function(t, par, log, ...) {
        if (log) {
            log(x$rate(t, par, ...)) - H(t, par, ...)
        } else {
            x$rate(t, par, ...) * exp(-H(t, par, ...))
        }
    }
    function(t, par = NULL, log = FALSE, ...) {
        par <- get_params(par, x$par)
        # Vectorize over t for integration compatibility
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

#' Method for obtaining the hazard function of a `dfr_dist` object.
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass into the `integrate` function.
#' @return A function that computes the hazard function of the distribution.
#' It accepts `t`, the time at which to compute the hazard function, and
#' `par`, the parameters of the distribution. If `par` is `NULL`, then the
#' parameters of the `dfr_dist` object `x` are used. Finally, it passes any
#' additional arguments `...` to the `rate` function of the `dfr_dist`
#' object `x`.
#' @importFrom stats integrate
#' @importFrom utils modifyList
#' @export
cum_haz.dfr_dist <- function(x, ...) {

    integrator_defaults <- list(
        lower = 0, subdivisions = 1000L, abs.tol = 1e-3)
    integrator <- modifyList(integrator_defaults, list(...))

    function(t, par = NULL, ...) {
        par <- get_params(par, x$par)
        # we call `integrate` with the options in `int_options` and
        # we integrate `rate` from 0 to `t` with respect to the parameters
        # `par` and any additional arguments `...`
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
    function (t = 0, par = NULL, log.p = FALSE, ...) {
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
}

#' Log-likelihood method for `dfr_dist` objects
#'
#' Returns a function that computes the log-likelihood of the data given
#' the distribution parameters. The log-likelihood for survival data is:
#'
#' For exact observations (uncensored): log(f(t)) = log(h(t)) - H(t)
#' For right-censored observations: log(S(t)) = -H(t)
#'
#' where h(t) is the hazard function, H(t) is the cumulative hazard,
#' f(t) = h(t)*S(t) is the pdf, and S(t) = exp(-H(t)) is the survival function.
#'
#' @param model The `dfr_dist` object
#' @param ... Additional arguments to pass to the hazard and cumulative hazard
#' @return A function that computes the log-likelihood. It accepts:
#'         - `df`: A data frame with observation times and censoring indicators
#'         - `par`: The parameters of the distribution
#'         - `...`: Additional arguments passed to internal functions
#' @importFrom likelihood.model loglik
#' @export
loglik.dfr_dist <- function(model, ...) {
    H <- cum_haz(model, ...)

    function(df, par = NULL, ...) {
        par <- get_params(par, model$par)

        ob_col <- model$ob_col
        delta_col <- model$delta_col

        # Extract observation times
        t <- df[[ob_col]]

        # Extract event indicators (1 = exact observation, 0 = right-censored)
        # If column doesn't exist, assume all exact observations
        if (delta_col %in% names(df)) {
            delta <- df[[delta_col]]
        } else {
            delta <- rep(1, length(t))
        }

        # Compute log-likelihood contributions
        ll <- 0

        # For exact observations (delta = 1): log(h(t)) - H(t)
        exact_idx <- which(delta == 1)
        if (length(exact_idx) > 0) {
            t_exact <- t[exact_idx]
            h_exact <- sapply(t_exact, function(ti) model$rate(ti, par, ...))
            H_exact <- sapply(t_exact, function(ti) H(ti, par, ...))
            ll <- ll + sum(log(h_exact) - H_exact)
        }

        # For censored observations (delta = 0): -H(t)
        censor_idx <- which(delta == 0)
        if (length(censor_idx) > 0) {
            t_censor <- t[censor_idx]
            H_censor <- sapply(t_censor, function(ti) H(ti, par, ...))
            ll <- ll - sum(H_censor)
        }

        ll
    }
}
