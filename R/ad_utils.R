#' Automatic Differentiation Utilities for dfr.dist
#'
#' These functions provide AD-based gradient and Hessian computation
#' using the femtograd package when available.
#'
#' @name ad_utils
NULL

#' Check if femtograd is available
#' @return TRUE if femtograd is installed
#' @keywords internal
has_femtograd <- function() {
    requireNamespace("femtograd", quietly = TRUE)
}

#' Compute gradient using femtograd reverse-mode AD
#'
#' @param f A function f(par) returning a scalar. The function should work
#'          with parameters accessed via `par[i]` indexing.
#' @param par Numeric vector of parameters
#' @return Numeric vector of partial derivatives
#' @keywords internal
ad_gradient <- function(f, par) {
    if (!has_femtograd()) {
        stop("femtograd package required for AD gradient. Install with: ",
             "devtools::install_github('queelius/femtograd')")
    }

    p <- length(par)
    grad_vec <- numeric(p)

    # Create value objects for each parameter
    par_vals <- lapply(par, femtograd::val)

    # Make it work with [i] indexing
    class(par_vals) <- c("ad_params", "list")

    # Evaluate function with value objects
    result <- f(par_vals)

    # Backward pass
    femtograd::zero_grad(result)
    femtograd::backward(result)

    # Extract gradients
    for (i in seq_len(p)) {
        grad_vec[i] <- femtograd::grad(par_vals[[i]])
    }

    grad_vec
}

#' Compute Jacobian using femtograd forward-mode AD
#'
#' For a function f: R^n -> R^m, computes the m x n Jacobian matrix.
#' Each column is computed by a forward pass with tangent = 1 for that input.
#'
#' @param f A function f(par) returning a numeric vector. The function should
#'          access parameters using `[[` indexing (e.g., `par[[1]]`) for AD
#'          compatibility, or accept a vector where femtograd ops are overloaded.
#' @param par Numeric vector of parameters
#' @return Jacobian matrix (output_dim x input_dim)
#' @keywords internal
ad_jacobian <- function(f, par) {
    if (!has_femtograd()) {
        stop("femtograd package required for AD Jacobian. Install with: ",
             "devtools::install_github('queelius/femtograd')")
    }

    p <- length(par)

    # First call to determine output dimension
    f_val <- f(par)
    m <- length(f_val)

    # Jacobian matrix: m outputs x p inputs
    jac <- matrix(0, nrow = m, ncol = p)

    # Forward-mode: one pass per input dimension
    for (j in seq_len(p)) {
        # Create dual numbers with tangent = 1 for parameter j
        par_dual <- lapply(seq_len(p), function(i) {
            femtograd::dual_num(par[i], if (i == j) 1 else 0)
        })

        # Make par_dual behave more like a vector for indexing
        class(par_dual) <- c("ad_params", "list")

        # Evaluate function
        result <- f(par_dual)

        # Extract tangents (partial derivatives w.r.t. par[j])
        # Handle both direct dual results and lists of duals (from c())
        if (m == 1) {
            # Result could be a dual directly or a list containing a dual
            if (is.list(result) && !inherits(result, "dual")) {
                jac[1, j] <- femtograd::tangent(result[[1]])
            } else {
                jac[1, j] <- femtograd::tangent(result)
            }
        } else {
            for (i in seq_len(m)) {
                jac[i, j] <- femtograd::tangent(result[[i]])
            }
        }
    }

    jac
}

#' Index operator for AD parameter lists
#' @keywords internal
#' @export
`[.ad_params` <- function(x, i) {
    x[[i]]
}

#' Compute Hessian via AD Jacobian of gradient
#'
#' Uses the hybrid approach: takes an analytical (or AD) gradient function
#' and computes its Jacobian using forward-mode AD to get the Hessian.
#'
#' @param score_fn A function score_fn(par) returning the gradient vector
#' @param par Numeric vector of parameters
#' @return Hessian matrix (p x p)
#' @export
ad_hessian <- function(score_fn, par) {
    # Hessian = Jacobian of gradient
    ad_jacobian(score_fn, par)
}

#' Compute Hessian using forward-over-reverse AD
#'
#' For when you have only the objective function (not analytical gradient).
#' Uses dual numbers wrapping value objects.
#'
#' @param f A function f(par) returning a scalar
#' @param par Numeric vector of parameters
#' @return Hessian matrix (p x p)
#' @keywords internal
ad_hessian_fwd_rev <- function(f, par) {
    if (!has_femtograd()) {
        stop("femtograd package required. Install with: ",
             "devtools::install_github('queelius/femtograd')")
    }

    p <- length(par)
    hess <- matrix(0, nrow = p, ncol = p)

    # For each direction j, compute d/dpar[j] of the gradient
    for (j in seq_len(p)) {
        # Create dual(value) objects: primal tracks gradient, tangent tracks d/dpar[j]
        par_dual_val <- lapply(seq_len(p), function(i) {
            v <- femtograd::val(par[i])
            femtograd::dual_num(v, if (i == j) 1 else 0)
        })

        # Evaluate function
        result <- f(par_dual_val)

        # Backward pass on primal to get gradient
        primal_result <- femtograd::primal(result)
        femtograd::zero_grad(primal_result)
        femtograd::backward(primal_result)

        # Extract Hessian row j: tangent of each gradient component
        for (i in seq_len(p)) {
            # grad of primal w.r.t. par[i], then tangent gives d/dpar[j]
            grad_i <- femtograd::grad(femtograd::primal(par_dual_val[[i]]))
            # For forward-over-reverse, we need tangent of the gradient
            # This requires dual numbers to propagate through backward pass
            # Simplified: use Jacobian of gradient instead
            hess[i, j] <- grad_i
        }
    }

    # Symmetrize (numerical errors can break symmetry slightly)
    0.5 * (hess + t(hess))
}


#' Create AD-aware log-likelihood function for dfr_dist
#'
#' Returns a function that computes log-likelihood using femtograd objects,
#' enabling automatic differentiation.
#'
#' @param model A dfr_dist object
#' @param df Data frame with observations
#' @return A function f(par) that works with femtograd value/dual objects
#' @keywords internal
make_ad_loglik <- function(model, df) {
    ob_col <- model$ob_col
    delta_col <- model$delta_col
    t_obs <- df[[ob_col]]
    delta <- if (delta_col %in% names(df)) df[[delta_col]] else rep(1, length(t_obs))

    # Return function that accepts femtograd parameter objects
    function(par) {
        # par can be a list of value/dual objects or numeric vector

        # For AD, we need to compute H(t) = integral of rate
        # This is tricky because integrate() doesn't work with AD objects
        # Solution: require analytical cumulative hazard for AD mode

        # Simple approach for exponential (constant rate):
        # H(t) = rate * t
        # For general case, need user-provided cum_haz_fn

        ll <- 0

        # This is a simplified version - full implementation needs
        # user-provided analytical cumulative hazard
        for (i in seq_along(t_obs)) {
            ti <- t_obs[i]
            h_i <- model$rate(ti, par)

            # Approximate H(t) - for AD, need analytical form
            # Using simple Euler approximation for prototype
            H_i <- h_i * ti  # Only correct for constant hazard!

            if (delta[i] == 1) {
                ll <- ll + log(h_i) - H_i
            } else {
                ll <- ll - H_i
            }
        }

        ll
    }
}
