#' Diagnostic Methods for DFR Distributions
#'
#' Methods for assessing model fit and visualizing survival distributions.
#'
#' @name diagnostics
#' @family dfr_dist
#' @importFrom graphics abline lines
#' @importFrom stats ppoints qexp qqplot residuals
NULL

# =============================================================================
# Residuals
# =============================================================================

#' Residuals for DFR Distribution Fits
#'
#' Computes residuals for assessing model fit of a DFR distribution to
#' survival data. Cox-Snell residuals should follow Exp(1) if the model
#' is correct. Martingale residuals identify observations poorly fit by
#' the model.
#'
#' @param object A `dfr_dist` object
#' @param data Data frame with survival data (must have time column and
#'   optionally delta column for censoring indicator)
#' @param par Parameter vector. If NULL, uses object's stored parameters.
#' @param type Type of residual:
#'   \describe{
#'     \item{"cox-snell"}{H(t_i) - should follow Exp(1) if model correct}
#'     \item{"martingale"}{delta_i - H(t_i) - useful for identifying outliers}
#'   }
#' @param ... Additional arguments passed to cum_haz
#'
#' @return Numeric vector of residuals, same length as data
#'
#' @details
#' **Cox-Snell residuals** are defined as r_i = H(t_i), the cumulative hazard
#' evaluated at the observation time. If the fitted model is correct, these
#' should follow an Exp(1) distribution (possibly censored).
#'
#' **Martingale residuals** are defined as M_i = delta_i - H(t_i), where
#' delta_i is the event indicator. They sum to zero and can identify
#' observations that are poorly fit. Large positive values indicate observations
#' that failed "too early" relative to the model; large negative values indicate
#' observations that survived "too long".
#'
#' @section Diagnostic Use:
#' \itemize{
#'   \item Q-Q plot of Cox-Snell residuals against Exp(1) to check overall fit
#'   \item Plot Martingale residuals vs. covariates to check functional form
#'   \item Plot Martingale residuals vs. fitted values to check homogeneity
#' }
#'
#' @examples
#' # Fit exponential to simulated data
#' set.seed(42)
#' df <- data.frame(t = rexp(100, rate = 0.5), delta = 1)
#' exp_dist <- dfr_exponential(lambda = 0.5)
#'
#' # Cox-Snell residuals
#' cs_resid <- residuals(exp_dist, df, type = "cox-snell")
#'
#' # Should follow Exp(1) - check with Q-Q plot
#' qqplot(qexp(ppoints(100)), sort(cs_resid),
#'        main = "Cox-Snell Residuals Q-Q Plot",
#'        xlab = "Theoretical Exp(1)", ylab = "Sample")
#' abline(0, 1, col = "red")
#'
#' # Martingale residuals
#' mart_resid <- residuals(exp_dist, df, type = "martingale")
#' summary(mart_resid)  # Should sum to approximately 0
#'
#' @export
residuals.dfr_dist <- function(object, data, par = NULL,
                                type = c("cox-snell", "martingale"), ...) {
    type <- match.arg(type)
    H <- cum_haz(object, ...)

    if (is.null(par)) {
        par <- object$par
        if (is.null(par)) {
            stop("Parameters required: provide via 'par' argument or in distribution object")
        }
    }

    if (!object$ob_col %in% names(data)) {
        stop(sprintf("Time column '%s' not found in data", object$ob_col))
    }
    t <- data[[object$ob_col]]

    H_vals <- sapply(t, function(ti) H(ti, par, ...))

    if (type == "cox-snell") {
        return(H_vals)
    }

    if (object$delta_col %in% names(data)) {
        delta <- data[[object$delta_col]]
    } else {
        delta <- rep(1, nrow(data))
    }

    delta - H_vals
}

# =============================================================================
# Plot Methods
# =============================================================================

#' Plot DFR Distribution Functions
#'
#' Visualizes the survival, hazard, or cumulative hazard function of a
#' DFR distribution. Optionally overlays empirical estimates from data.
#'
#' @param x A `dfr_dist` object
#' @param data Optional data frame with survival data for empirical overlay
#' @param par Parameter vector. If NULL, uses object's stored parameters.
#' @param what Which function to plot:
#'   \describe{
#'     \item{"survival"}{S(t) = exp(-H(t))}
#'     \item{"hazard"}{h(t) - instantaneous failure rate}
#'     \item{"cumhaz"}{H(t) - cumulative hazard}
#'   }
#' @param xlim x-axis limits. If NULL, determined from data or defaults to c(0, 10).
#' @param n Number of points for smooth curve (default 200)
#' @param add If TRUE, add to existing plot
#' @param col Line color for theoretical curve
#' @param lwd Line width for theoretical curve
#' @param empirical If TRUE and data provided, overlay Kaplan-Meier estimate
#' @param empirical_col Color for empirical curve
#' @param ... Additional arguments passed to plot()
#'
#' @return Invisibly returns the plotted values as a list with elements
#'   `t` (time points) and `y` (function values).
#'
#' @details
#' When `empirical = TRUE` and data is provided, overlays:
#' \itemize{
#'   \item For survival: Kaplan-Meier estimate (step function)
#'   \item For cumhaz: Nelson-Aalen estimate (step function)
#'   \item For hazard: Kernel-smoothed hazard estimate
#' }
#'
#' @examples
#' # Plot survival function for Weibull distribution
#' weib <- dfr_weibull(shape = 2, scale = 5)
#' plot(weib, what = "survival", xlim = c(0, 10))
#'
#' # Overlay hazard functions for different shapes
#' plot(weib, what = "hazard", xlim = c(0, 10), col = "blue")
#' weib_k1 <- dfr_weibull(shape = 1, scale = 5)  # Exponential
#' plot(weib_k1, what = "hazard", add = TRUE, col = "green")
#' weib_k3 <- dfr_weibull(shape = 3, scale = 5)  # Steeper wear-out
#' plot(weib_k3, what = "hazard", add = TRUE, col = "red")
#' legend("topleft", c("k=2", "k=1 (exp)", "k=3"),
#'        col = c("blue", "green", "red"), lwd = 2)
#'
#' # Compare fitted model to data
#' set.seed(123)
#' true_weib <- dfr_weibull(shape = 2.5, scale = 10)
#' sim_data <- data.frame(t = sampler(true_weib)(100), delta = 1)
#' solver <- fit(dfr_weibull())
#' result <- solver(sim_data, par = c(2, 8))
#' fitted_weib <- dfr_weibull(shape = coef(result)[1], scale = coef(result)[2])
#' plot(fitted_weib, data = sim_data, what = "survival",
#'      xlim = c(0, max(sim_data$t)), empirical = TRUE)
#'
#' @export
plot.dfr_dist <- function(x, data = NULL, par = NULL,
                           what = c("survival", "hazard", "cumhaz"),
                           xlim = NULL, n = 200, add = FALSE,
                           col = "black", lwd = 2, empirical = TRUE,
                           empirical_col = "steelblue", ...) {
    what <- match.arg(what)

    if (is.null(par)) {
        par <- x$par
        if (is.null(par)) {
            stop("Parameters required: provide via 'par' argument or in distribution object")
        }
    }

    if (is.null(xlim)) {
        xlim <- if (!is.null(data) && x$ob_col %in% names(data)) {
            c(0, max(data[[x$ob_col]]) * 1.1)
        } else {
            c(0, 10)
        }
    }

    t_seq <- seq(xlim[1] + 1e-6, xlim[2], length.out = n)

    fn <- switch(what,
                 survival = surv(x),
                 hazard = hazard(x),
                 cumhaz = cum_haz(x))

    y_vals <- sapply(t_seq, function(ti) fn(ti, par))

    ylab <- switch(what,
                   survival = "S(t)",
                   hazard = "h(t)",
                   cumhaz = "H(t)")

    if (!add) {
        plot(t_seq, y_vals, type = "l", col = col, lwd = lwd,
             xlim = xlim, xlab = "Time", ylab = ylab, ...)
    } else {
        lines(t_seq, y_vals, col = col, lwd = lwd)
    }

    if (!is.null(data) && empirical && !add && x$ob_col %in% names(data)) {
        t <- data[[x$ob_col]]
        delta <- if (x$delta_col %in% names(data)) {
            data[[x$delta_col]]
        } else {
            rep(1, length(t))
        }

        km <- kaplan_meier(t, delta)

        if (what == "survival") {
            lines(km$time, km$surv, type = "s", col = empirical_col, lwd = 1.5)
        } else if (what == "cumhaz") {
            lines(km$time, km$cumhaz, type = "s", col = empirical_col, lwd = 1.5)
        }
    }

    invisible(list(t = t_seq, y = y_vals))
}

# =============================================================================
# Kaplan-Meier / Nelson-Aalen (internal helper)
# =============================================================================

#' Compute Kaplan-Meier and Nelson-Aalen Estimates
#'
#' Internal function to compute non-parametric survival estimates.
#'
#' @param time Observation times
#' @param delta Event indicators (1 = event, 0 = censored)
#' @return List with time, surv, cumhaz
#' @keywords internal
kaplan_meier <- function(time, delta) {
    ord <- order(time)
    t_sorted <- time[ord]
    d_sorted <- delta[ord]

    unique_times <- unique(t_sorted[d_sorted == 1])

    if (length(unique_times) == 0) {
        return(list(time = numeric(0), surv = numeric(0), cumhaz = numeric(0)))
    }

    surv <- numeric(length(unique_times))
    cumhaz <- numeric(length(unique_times))
    S <- 1
    H <- 0

    for (i in seq_along(unique_times)) {
        ti <- unique_times[i]
        n_risk <- sum(t_sorted >= ti)
        d_i <- sum(t_sorted == ti & d_sorted == 1)

        S <- S * (1 - d_i / n_risk)
        H <- H + d_i / n_risk

        surv[i] <- S
        cumhaz[i] <- H
    }

    list(time = unique_times, surv = surv, cumhaz = cumhaz)
}

# =============================================================================
# Q-Q Plot for Residuals
# =============================================================================

#' Q-Q Plot for Cox-Snell Residuals
#'
#' Creates a Q-Q plot comparing Cox-Snell residuals to the theoretical
#' Exp(1) distribution. A good fit shows points along the diagonal.
#'
#' @param object A `dfr_dist` object
#' @param data Data frame with survival data
#' @param par Parameter vector. If NULL, uses object's stored parameters.
#' @param add_line If TRUE, adds reference line at y = x
#' @param ... Additional arguments passed to residuals and qqplot
#'
#' @return Invisibly returns Cox-Snell residuals
#'
#' @details
#' Cox-Snell residuals r_i = H(t_i) should follow an Exp(1) distribution
#' if the model is correctly specified. Departure from the diagonal line
#' in the Q-Q plot indicates model misspecification:
#'
#' \itemize{
#'   \item Points above the line: observations failed earlier than expected
#'   \item Points below the line: observations survived longer than expected
#'   \item Systematic curvature: wrong distributional form
#' }
#'
#' @examples
#' # Check fit of exponential model
#' set.seed(42)
#' df <- data.frame(t = rexp(100, rate = 0.5), delta = 1)
#' exp_dist <- dfr_exponential(lambda = 0.5)
#'
#' qqplot_residuals(exp_dist, df)
#'
#' # Check fit with wrong model (Weibull data, exponential fit)
#' df_weib <- data.frame(t = sampler(dfr_weibull(shape = 2, scale = 5))(100), delta = 1)
#' exp_fit <- dfr_exponential(lambda = 1 / mean(df_weib$t))  # Moment estimate
#' qqplot_residuals(exp_fit, df_weib)  # Should show systematic departure
#'
#' @export
qqplot_residuals <- function(object, data, par = NULL, add_line = TRUE, ...) {
    cs_resid <- residuals(object, data, par = par, type = "cox-snell", ...)

    if (object$delta_col %in% names(data)) {
        uncensored_resid <- cs_resid[data[[object$delta_col]] == 1]
    } else {
        uncensored_resid <- cs_resid
    }

    theoretical <- qexp(ppoints(length(uncensored_resid)))

    qqplot(theoretical, sort(uncensored_resid),
           main = "Cox-Snell Residuals Q-Q Plot",
           xlab = "Theoretical Exp(1) Quantiles",
           ylab = "Sample Quantiles",
           pch = 16, col = "steelblue")

    if (add_line) {
        abline(0, 1, col = "red", lwd = 2)
    }

    invisible(cs_resid)
}
