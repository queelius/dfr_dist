get_params <- function(par, default = NULL) {
    if (is.null(par) || all(is.na(par))) {
        return(default)
    }
    if (!is.null(default) && length(par) != length(default)) {
        stop(sprintf(
            "Parameter length mismatch: got %d, expected %d",
            length(par), length(default)
        ))
    }
    if (!is.null(default)) {
        par[is.na(par)] <- default[is.na(par)]
    }
    par
}

require_params <- function(par, default) {
    par <- get_params(par, default)
    if (is.null(par)) {
        stop("Parameters required: provide via 'par' argument or in distribution object")
    }
    par
}

get_delta <- function(df, delta_col = "delta") {
    if (delta_col %in% names(df)) df[[delta_col]] else rep(1, nrow(df))
}
