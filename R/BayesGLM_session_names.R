#' Get `session_names` for GLM
#' 
#' @keywords internal
BayesGLM_session_names <- function(nS, session_names, BOLD_names, design_names){
  # In order of priority:
  #   1) the session names argument
  #   2) the BOLD data 
  #   3) the design (quickly check if it's a list)
  if (!is.null(session_names)) {
    stopifnot(is.character(session_names))
    stopifnot(!any(duplicated(session_names)))
    if (length(session_names) != nS) { 
      stop("`session_names` length should match the number of session of BOLD data, ", nS, ".")
    }
    if (!is.null(BOLD_names) && !all(session_names == BOLD_names)) {
      warning("Using `session_names` rather than the names of `BOLD`.")
    }
    if (!is.null(design_names) && !all(session_names == design_names)) {
      warning("Using `session_names` rather than the names of `design`.")
    }
  } else if (!is.null(BOLD_names)) {
    stopifnot(!any(duplicated(BOLD_names)))
    #cat("Using names of `BOLD` as session names.\n")
    session_names <- BOLD_names
    if (length(BOLD_names) != nS) { 
      stop("`BOLD_names` length should match the number of session of BOLD data, ", nS, ".")
    }
    if (!is.null(design_names) && !all(session_names == design_names)) {
      warning("Using `session_names` rather than the names of `design`.")
    }
  } else if (!is.null(design_names)) {
    stopifnot(!any(duplicated(design_names)))
    #cat("Using names of `design` as session names.\n")
    session_names <- design_names
    if (length(design_names) != nS) { 
      stop("`design_names` length should match the number of session of BOLD data, ", nS, ".")
    }
  } else {
    session_names <- if (nS==1) { "single_sess" } else { paste0("sess_", seq(nS)) }
  }
  session_names
}