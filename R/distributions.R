#' Classic Survival Distribution Constructors
#'
#' Convenient constructors for commonly-used survival distributions.
#' Each provides the complete specification (rate, cum_haz_rate, score_fn,
#' and where practical, hess_fn) for optimal performance.
#'
#' @section Left-Censoring Note:
#' The analytical score and Hessian functions provided by these constructors
#' assume event indicators in \{0, 1\} (right-censored and exact observations).
#' For left-censored data (delta = -1), these functions are not applicable and
#' the package automatically falls back to numerical differentiation via
#' \code{numDeriv::grad} and \code{numDeriv::hessian} through the log-likelihood,
#' which handles all censoring types correctly.
#'
#' @name distributions
#' @family distributions
NULL

# =============================================================================
# Exponential Distribution
# =============================================================================

#' Exponential Distribution (Constant Hazard)
#'
#' Creates a DFR distribution with constant failure rate (exponential).
#' The exponential distribution is "memoryless" - the hazard does not depend

#' on time, making it appropriate for random failures unrelated to age.
#'
#' @param lambda Rate parameter (failure rate). If NULL, must be provided
#'   when calling methods or fitting. Must be positive.
#'
#' @details
#' The exponential distribution has:
#' \itemize{
#'   \item Hazard: \eqn{h(t) = \lambda}
#'   \item Cumulative hazard: \eqn{H(t) = \lambda t}
#'   \item Survival: \eqn{S(t) = e^{-\lambda t}}
#'   \item Mean time to failure: \eqn{1/\lambda}
#' }
#'
#' @section Reliability Interpretation:
#' Use exponential for:
#' \itemize{
#'   \item Electronic components during useful life (random failures)
#'   \item Systems with redundancy where failures are independent
#'   \item As a baseline model to test against more complex alternatives
#' }
#'
#' @return A `dfr_dist` object with analytical rate, cumulative hazard,
#'   and score function.
#'
#' @examples
#' # Component with MTBF of 1000 hours (lambda = 0.001)
#' comp <- dfr_exponential(lambda = 0.001)
#'
#' # Survival probability at 500 hours
#' S <- surv(comp)
#' S(500)  # ~60.6%
#'
#' # Fit to failure data
#' set.seed(42)
#' failures <- data.frame(t = rexp(50, rate = 0.001), delta = 1)
#' solver <- fit(comp)
#' result <- solver(failures, par = c(0.002))
#' coef(result)  # Should be close to 0.001
#'
#' @export
dfr_exponential <- function(lambda = NULL) {
    dfr_dist(
        rate = function(t, par, ...) {
            rep(par[[1]], length(t))
        },
        cum_haz_rate = function(t, par, ...) {
            par[[1]] * t
        },
        score_fn = function(df, par, ...) {
            delta <- get_delta(df)
            c(sum(delta == 1) / par[[1]] - sum(df$t))
        },
        hess_fn = function(df, par, ...) {
            delta <- get_delta(df)
            matrix(-sum(delta == 1) / par[[1]]^2, nrow = 1, ncol = 1)
        },
        par = lambda
    )
}

# =============================================================================
# Weibull Distribution
# =============================================================================

#' Weibull Distribution (Power-Law Hazard)
#'
#' Creates a DFR distribution with Weibull hazard function.
#' The Weibull is extremely versatile: it can model increasing (wear-out),
#' decreasing (infant mortality), or constant (exponential) failure rates.
#'
#' @param shape Shape parameter (k). Controls hazard behavior:
#'   k < 1: decreasing hazard (infant mortality)
#'   k = 1: constant hazard (exponential)
#'   k > 1: increasing hazard (wear-out)
#' @param scale Scale parameter (sigma). Controls time scale.
#'
#' @details
#' The Weibull distribution has:
#' \itemize{
#'   \item Hazard: \eqn{h(t) = (k/\sigma)(t/\sigma)^{k-1}}
#'   \item Cumulative hazard: \eqn{H(t) = (t/\sigma)^k}
#'   \item Survival: \eqn{S(t) = e^{-(t/\sigma)^k}}
#'   \item Characteristic life (63.2% failure): \eqn{\sigma}
#' }
#'
#' @section Reliability Interpretation:
#' \itemize{
#'   \item Shape < 1: Infant mortality (burn-in failures, defects)
#'   \item Shape = 1: Random failures (reduces to exponential)
#'   \item Shape = 2: Rayleigh distribution (linear hazard increase)
#'   \item Shape > 2: Accelerating wear-out (fatigue, corrosion)
#' }
#'
#' @section B-Life Calculation:
#' The B10 life (10% failure quantile) is commonly used in reliability:
#' \eqn{B10 = \sigma \cdot (-\log(0.9))^{1/k}}
#'
#' @return A `dfr_dist` object with analytical rate, cumulative hazard,
#'   and score function.
#'
#' @examples
#' # Bearing with wear-out failure (shape > 1)
#' bearing <- dfr_weibull(shape = 2.5, scale = 50000)
#'
#' # Hazard increases with time
#' h <- hazard(bearing)
#' h(10000)  # hazard at 10k hours
#' h(40000)  # much higher at 40k hours
#'
#' # B10 life calculation
#' Q <- inv_cdf(bearing)
#' B10 <- Q(0.10)  # 10% failure quantile
#'
#' # Fit to test data with right-censoring
#' set.seed(123)
#' test_data <- data.frame(
#'   t = pmin(rweibull(100, shape = 2.5, scale = 50000), 30000),
#'   delta = as.integer(rweibull(100, shape = 2.5, scale = 50000) <= 30000)
#' )
#' solver <- fit(dfr_weibull())
#' result <- solver(test_data, par = c(2, 40000))
#' coef(result)
#'
#' @export
dfr_weibull <- function(shape = NULL, scale = NULL) {
    par <- if (!is.null(shape) && !is.null(scale)) c(shape, scale) else NULL

    dfr_dist(
        rate = function(t, par, ...) {
            k <- par[[1]]
            sigma <- par[[2]]
            (k / sigma) * (t / sigma)^(k - 1)
        },
        cum_haz_rate = function(t, par, ...) {
            (t / par[[2]])^par[[1]]
        },
        score_fn = function(df, par, ...) {
            k <- par[[1]]
            sigma <- par[[2]]
            t <- df$t
            delta <- get_delta(df)

            n_events <- sum(delta == 1)
            t_ratio <- t / sigma
            log_t_ratio <- log(t_ratio)

            dk <- n_events / k + sum(delta * log_t_ratio) -
                sum(t_ratio^k * log_t_ratio)
            dsigma <- -n_events * k / sigma +
                (k / sigma) * sum(t_ratio^k)

            c(dk, dsigma)
        },
        hess_fn = function(df, par, ...) {
            k <- par[[1]]
            sigma <- par[[2]]
            t <- df$t
            delta <- get_delta(df)

            n_events <- sum(delta == 1)
            t_ratio <- t / sigma
            log_t_ratio <- log(t_ratio)
            t_ratio_k <- t_ratio^k

            d2k <- -n_events / k^2 - sum(t_ratio_k * log_t_ratio^2)
            d2k_sigma <- -n_events / sigma +
                (1 / sigma) * sum(t_ratio_k) +
                (k / sigma) * sum(t_ratio_k * log_t_ratio)
            d2sigma <- n_events * k / sigma^2 -
                k * (k + 1) / sigma^2 * sum(t_ratio_k)

            matrix(c(d2k, d2k_sigma, d2k_sigma, d2sigma), nrow = 2, ncol = 2)
        },
        par = par
    )
}

# =============================================================================
# Gompertz Distribution
# =============================================================================

#' Gompertz Distribution (Exponential Growth Hazard)
#'
#' Creates a DFR distribution with Gompertz hazard function.
#' The Gompertz models exponentially increasing failure rate, often used
#' for biological aging and wear-out processes that accelerate over time.
#'
#' @param a Initial hazard rate at t=0. Must be positive.
#' @param b Growth rate of the hazard. Must be positive.
#'
#' @details
#' The Gompertz distribution has:
#' \itemize{
#'   \item Hazard: \eqn{h(t) = a \cdot e^{bt}}
#'   \item Cumulative hazard: \eqn{H(t) = (a/b)(e^{bt} - 1)}
#'   \item Survival: \eqn{S(t) = \exp(-(a/b)(e^{bt} - 1))}
#' }
#'
#' @section Reliability Interpretation:
#' Use Gompertz for:
#' \itemize{
#'   \item Aging systems where failure rate grows exponentially
#'   \item Biological mortality (human lifespans)
#'   \item Corrosion/degradation with accelerating kinetics
#' }
#'
#' When b is small, Gompertz approximates exponential early in life.
#' As b increases, wear-out acceleration becomes more pronounced.
#'
#' @return A `dfr_dist` object with analytical rate, cumulative hazard,
#'   and score function.
#'
#' @examples
#' # Aging system: initial hazard 0.001, doubling every 1000 hours
#' # b = log(2)/1000 gives doubling time of 1000
#' system <- dfr_gompertz(a = 0.001, b = log(2)/1000)
#'
#' # Hazard at various ages
#' h <- hazard(system)
#' h(0)      # 0.001 (initial)
#' h(1000)   # 0.002 (doubled)
#' h(2000)   # 0.004 (quadrupled)
#'
#' # Survival probability
#' S <- surv(system)
#' S(5000)   # probability of surviving 5000 hours
#'
#' @export
dfr_gompertz <- function(a = NULL, b = NULL) {
    par <- if (!is.null(a) && !is.null(b)) c(a, b) else NULL

    dfr_dist(
        rate = function(t, par, ...) {
            par[[1]] * exp(par[[2]] * t)
        },
        cum_haz_rate = function(t, par, ...) {
            (par[[1]] / par[[2]]) * (exp(par[[2]] * t) - 1)
        },
        score_fn = function(df, par, ...) {
            a <- par[[1]]
            b <- par[[2]]
            t <- df$t
            delta <- get_delta(df)

            exp_bt <- exp(b * t)
            n_events <- sum(delta == 1)

            da <- n_events / a - (1 / b) * sum(exp_bt - 1)
            db <- sum(delta * t) + (a / b^2) * sum(exp_bt - 1) -
                (a / b) * sum(t * exp_bt)

            c(da, db)
        },
        par = par
    )
}

# =============================================================================
# Log-Logistic Distribution
# =============================================================================

#' Log-Logistic Distribution (Non-Monotonic Hazard)
#'
#' Creates a DFR distribution with log-logistic hazard function.
#' The log-logistic has a non-monotonic hazard that increases then decreases,
#' useful for modeling processes with an initial risk that diminishes.
#'
#' @param alpha Scale parameter. Median lifetime when beta > 1.
#' @param beta Shape parameter. Controls hazard shape:
#'   beta <= 1: monotonically decreasing hazard
#'   beta > 1: hazard increases to a peak then decreases
#'
#' @details
#' The log-logistic distribution has:
#' \itemize{
#'   \item Hazard: \eqn{h(t) = \frac{(\beta/\alpha)(t/\alpha)^{\beta-1}}{1 + (t/\alpha)^\beta}}
#'   \item Survival: \eqn{S(t) = \frac{1}{1 + (t/\alpha)^\beta}}
#'   \item Median: \eqn{\alpha} (when beta > 1)
#' }
#'
#' @section Reliability Interpretation:
#' The log-logistic is useful when:
#' \itemize{
#'   \item Initial failures decrease after screening period
#'   \item Risk peaks early then declines (therapy response)
#'   \item Hazard is not monotonic throughout lifetime
#' }
#'
#' Note: The cumulative hazard has no closed form and is computed numerically.
#' For efficiency with large datasets, consider providing `cum_haz_rate`
#' using numerical integration cached appropriately.
#'
#' @return A `dfr_dist` object with analytical rate function.
#'   Cumulative hazard uses numerical integration.
#'
#' @examples
#' # Component with peak hazard around t = alpha
#' comp <- dfr_loglogistic(alpha = 1000, beta = 2)
#'
#' # Non-monotonic hazard
#' h <- hazard(comp)
#' h(500)   # increasing phase
#' h(1000)  # near peak
#' h(2000)  # decreasing phase
#'
#' # Survival function
#' S <- surv(comp)
#' S(1000)  # 50% survival at median (alpha)
#'
#' @export
dfr_loglogistic <- function(alpha = NULL, beta = NULL) {
    par <- if (!is.null(alpha) && !is.null(beta)) c(alpha, beta) else NULL

    dfr_dist(
        rate = function(t, par, ...) {
            alpha <- par[[1]]
            beta <- par[[2]]
            t_ratio <- t / alpha
            (beta / alpha) * t_ratio^(beta - 1) / (1 + t_ratio^beta)
        },
        cum_haz_rate = function(t, par, ...) {
            log(1 + (t / par[[1]])^par[[2]])
        },
        score_fn = function(df, par, ...) {
            alpha <- par[[1]]
            beta <- par[[2]]
            t <- df$t
            delta <- get_delta(df)

            t_ratio <- t / alpha
            t_ratio_beta <- t_ratio^beta
            log_t_ratio <- log(t_ratio)

            term_alpha <- t_ratio_beta / (1 + t_ratio_beta)
            dalpha <- sum(-delta * beta / alpha) +
                sum((1 + delta) * (beta / alpha) * term_alpha)

            term_beta <- log_t_ratio * t_ratio_beta / (1 + t_ratio_beta)
            dbeta <- sum(delta / beta) + sum(delta * log_t_ratio) -
                sum((1 + delta) * term_beta)

            c(dalpha, dbeta)
        },
        par = par
    )
}
