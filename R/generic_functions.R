#' Method for obtaining the cumulative hazard function of an object.
#' @param x The object to obtain the cumulative hazard function of.
#' @param ... Additional arguments to pass.
#' @return A function that computes the cumulative hazard function of the
#' distribution.
#' @export
cum_haz <- function(x, ...) {
    UseMethod("cum_haz", x)
}
