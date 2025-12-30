#' Message on attach
#' 
#' Checks if INLA is installed, and prints a warning if not.
#'
#' @param ... Not used
#' 
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
.onAttach <- function(...) {
  if (interactive()) {

    if (!requireNamespace("INLA", quietly = TRUE)) {
      warning(
        "The `INLA` package is required to the Bayesian methods in `BayesfMRI`.",
        "Please install it.", 
        call. = FALSE
      )
    }
  }

}
