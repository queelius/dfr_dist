get_params <- function(par, default = NULL) {
    if (is.null(par) || all(is.na(par))) {
        return(default)
    }
    if (!is.null(default)) {
        par[is.na(par)] <- default[is.na(par)]
    }
    par
}
