#' Validate an individual session in a \code{"BfMRI.sess"} object.
#' 
#' Check if object is valid for a list entry in \code{"BfMRI.sess"}.
#'
#' A valid entry in a \code{"BfMRI.sess"} object is a list with the following 
#'  named fields:
#'  \itemize{
#'    \item{"BOLD"}{\eqn{T \times V} matrix of BOLD responses, rows are time points, columns are voxels}
#'    \item{"design"}{\eqn{T \times K} matrix containing the K task regressors}
#'    \item{"nuisance"}{Optional. \eqn{T \times L} matrix containing the L nuisance regressors}
#'  }
#'
#' @param x The putative entry in a \code{"BfMRI.sess"} object.
#'
#' @return Logical. Is \code{x} a valid entry in a \code{"BfMRI.sess"} object?
#'
#' @keywords internal
is.a_session <- function(x){

  # `x` is a list.
  if (!is.list(x)) { message("`x` must be a list."); return(FALSE) }

  # `x` has fields 'BOLD', 'design', and maybe 'nuisance'.
  fields <- c("BOLD", "design", "nuisance")
  if (length(x) == 2) {
    if (!all(names(x) == fields[seq(2)])) {
      message("`x` should have fields 'BOLD' and 'design'."); return(FALSE)
    }
  } else if (length(x) == 3) {
    if (!all(names(x) == fields[seq(3)])) {
      message("`x` should have fields 'BOLD', 'design', and 'nuisance'.")
      return(FALSE)
    }
  } else {
    message("`x` should be length 2 or 3."); return(FALSE)
  }
  has_nus <- length(x) == 3

  # The data types are ok.
  if (!(is.numeric(x$BOLD) && is.matrix(x$BOLD))) { 
    message("`x$BOLD` must be a numeric matrix."); return(FALSE)
  }
  des_pw <- (!has_nus) && (inherits(x$design, "dgCMatrix"))
  if (!des_pw) {
    if (!(is.numeric(x$design) && is.matrix(x$design))) { 
      message("`x$design` must be a numeric matrix.")
      if (!has_nus) { 
        message("(It may also be a 'dgCMatrix' if prewhitening has been done.)")
      }
      return(FALSE)
    }
  }
  if (has_nus) {
    if (!(is.matrix(x$nuisance))) { 
      message("`x$nuisance` must be a matrix."); return(FALSE)
    }
  }

  # The dimensions are ok.
  if (des_pw) {
    nvox <- sum(!is.na(x$BOLD[1,]))
    if (nrow(x$BOLD) != nrow(x$design)/nvox) {
      message(
        "If prewhitening has been performed, the number of rows in the design ",
        "matrix should be T*V, where T=nrow(BOLD) and V is the number of ",
        "columns of BOLD that are non-NA."
      )
      return(FALSE)
    }
  } else {
    if ((nrow(x$BOLD) != nrow(x$design))) { 
      message("'BOLD' and 'design' must have the same number of rows (time points).")
      return(FALSE)
    }
  }
  if (has_nus) {
    if (nrow(x$BOLD) != nrow(x$nuisance)) { 
      message("'BOLD' and 'nuisance' must have the same number of rows (time points).")
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Validate a \code{"BfMRI.sess"} object.
#' 
#' Check if object is valid for a \code{"BfMRI.sess"} object.
#'
#' A valid entry in a \code{"BfMRI.sess"} object is a list of lists, each with 
#'  the following named fields:
#'  \itemize{
#'    \item{"BOLD"}{\eqn{T \times V} matrix of BOLD responses, rows are time points, columns are voxels}
#'    \item{"design"}{\eqn{T \times K} matrix containing the K task regressors}
#'    \item{"nuisance"}{Optional. \eqn{T \times J} matrix containing the L nuisance regressors}
#'  }
#'  In addition, all sessions must have the same number of data locations and tasks.
#'
#' @param x The putative \code{"BfMRI.sess"} object.
#'
#' @return Logical. Is \code{x} a valid \code{"BfMRI.sess"} object?
#'
#' @export
is.BfMRI.sess <- function(x){
  if (!is.list(x)) { message("`x` must be a list."); return(FALSE) }

  nS <- length(x)

  # Check the first session.
  if (!is.a_session(x[[1]])) { 
    message("The first entry of `x` is not a valid session.")
    return(FALSE)
  }

  # We're done now, if there's only one session.
  if (nS < 2) { return(TRUE) }

  # Check the rest of the sessions.
  all_is_sess <- vapply(x, is.a_session, FALSE)
  if (!all_is_sess) {
    message(sum(!all_is_sess), " out of ", nS, " entries in `x` are not valid sessions.")
    return(FALSE)
  }

  # Check that all sessions have the same number of data locations and tasks.
  nV <- ncol(x[[1]]$BOLD)
  nK <- ncol(x[[1]]$design)
  for (ii in seq(2, nS)) {
    if (ncol(x[[ii]]$BOLD != nV)) {
      message(
        "The first session has ", nV, " locations, but session ", ii, 
        " (and maybe others) does not. All sessions must have the same number ",
        "of locations."
      )
      return(FALSE)
    }
    if (ncol(x[[ii]]$design != nK)) {
      message(
        "The first session has ", nK, " locations, but session ", ii, 
        " (and maybe others) does not. All sessions must have the same number ",
        "of locations."
      )
      return(FALSE)
    }
  }

  return(TRUE)
}