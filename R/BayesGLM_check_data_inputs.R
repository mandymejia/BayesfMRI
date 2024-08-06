#' Is a valid design?
#'
#' @param design The design matrix/array
#' @keywords internal
BayesGLM_is_valid_one_design <- function(design) {
  stopifnot(is.matrix.or.df(design) || length(dim(design))==3)
  stopifnot(all(apply(design, 2, function(q){all(is.na(q)) || is.numeric(q)})))
  stopifnot(all(dim(design) > 0))
  TRUE
}

#' Is a valid nuisance?
#'
#' @param nuisance The nuisance matrix
#' @keywords internal
BayesGLM_is_valid_one_nuisance <- function(nuisance) {
  stopifnot(is.matrix.or.df(nuisance))
  stopifnot(all(apply(nuisance, 2, is.numeric)))
  stopifnot(all(dim(nuisance) > 0))
  TRUE
}

#' Is a valid scrub?
#'
#' @param scrub The scrub matrix
#' @keywords internal
BayesGLM_is_valid_one_scrub <- function(scrub) {
  if (is.null(scrub)) { return(TRUE) }
  if (is.logical(scrub)) { return(TRUE) }

  stopifnot(is.matrix.or.df(scrub))
  stopifnot(all(apply(scrub, 2, is.numeric)))
  stopifnot(all(dim(scrub) > 0))
  stopifnot(all(colSums(scrub) == 1))
  stopifnot(all(colSums(scrub == 0) == nrow(scrub)-1))
  TRUE
}

#' Format design
#'
#' Format design for \code{BayesGLM}, \code{fit_bayesglm},
#'  \code{multiGLM}, and \code{multiGLM_fun}.
#' @param design The \code{design} argument input. Will be formatted to a
#'  \code{nS}-length list.
#' @param scale_design Scale the design matrix by dividing each column by its
#'  maximum and then subtracting the mean? Default: \code{TRUE}. If
#'  \code{FALSE}, the design matrix is centered but not scaled.
#' @param nS_expect The expected number of sessions, if known.
#' @param nT_expect The expected number of timepoints, if known. For
#'  multi-session data this is a session-length vector.
#' @param nD_expect The expected number of designs, if known. For per-location
#'  modeling this is equal to \code{nVd0}, the initial number of data locations.
#'  For multi-session data this is a session-length vector.
#' @param per_location_design \code{FALSE} if per-location modeling is not
#'  being performed (i.e. for multiGLM); \code{TRUE} if it is; or, \code{NULL}
#'  to infer based on the dimensions of \code{design} (\code{TRUE} if the
#'  design has three dimensions.)
#' @keywords internal
#' @return \code{design}
BayesGLM_format_design <- function(
  design, scale_design=TRUE,
  nS_expect=NULL, nT_expect=NULL, nD_expect=NULL,
  per_location_design=NULL
  ){

  # Make `design` a sessions-length list of matrices or arrays.
  if (inherits(design, "BfMRI_design")) { design <- design$design }
  if (!is.list(design) || is.data.frame(design)) { design <- list(design) }
  if (all(vapply(design, inherits, FALSE, "BfMRI_design"))) { design <- lapply(design, '[[', "design") }
  nS <- length(design)
  if (!is.null(nS_expect)) {
    stopifnot(fMRItools::is_1(nS_expect, "numeric") && nS_expect==round(nS_expect))
    if (length(design) != nS_expect) {
      stop("`design` is length ", nS, ", but ", nS_expect, " sessions are expected.")
    }
  }

  if (is.null(per_location_design)) {
    per_location_design <- !is.matrix(design[[1]]) && !is.data.frame(design[[1]])
  }

  for (ss in seq(nS)) {
    stopifnot(BayesGLM_is_valid_one_design(design[[ss]]))
    if (is.data.frame(design[[ss]])) { design[[ss]] <- as.matrix(design[[ss]]) }
  }

  must_match_msg <- "This must match across sessions. Please fix."

  # Number of dimensions.
  des_dims <- lapply(design, dim)
  des_nDims <- lapply(des_dims, length)
  if (length(unique(des_nDims)) != 1) {
    stop("Some designs are matrices, and others are arrays. ", must_match_msg)
  }
  des_nDims <- des_nDims[[1]]
  des_is_array <- des_nDims == 3

  # First dimension, nT
  des_nT <- vapply(des_dims, '[', 0, 1)
  # Note that designs are allowed to have different numbers of timepoints.
  #   because BOLD data may be of unequal lengths.
  if (!is.null(nT_expect)) {
    stopifnot(is.numeric(nT_expect))
    if (length(nT_expect)==1) { nT_expect <- rep(nT_expect, nS) }
    stopifnot(length(nT_expect) == nS)
    for (ss in seq(nS)) {
      if (des_nT[ss] != nT_expect[ss]) {
        stop("The design has ", des_nT[ss], " locations, but ", nT_expect[ss],
          " locations are expected.")
      }
    }
  }
  nT <- des_nT; rm(des_nT)

  # Second dimension, nK
  des_nK <- vapply(des_dims, '[', 0, 2)
  if (length(unique(des_nK)) != 1) {
    stop("The designs have different numbers of tasks. ", must_match_msg)
    # (To indicate missing tasks, please instead use a column of `NA` values.)
  }
  des_nK <- des_nK[1]
  nK <- des_nK; rm(des_nK)

  ### Field names
  field_names <- lapply(design, function(q){dimnames(q)[[2]]})
  if (length(unique(field_names)) > 1) {
    stop("Field names (second dim. names of `design`) ",
    "should not differ across sessions.")
  }
  field_names <- field_names[[1]]
  # Set `field_names` if not provided.
  if (is.null(field_names)) {
    cat("Setting field names to 'field_1'")
    if (nK>1) { cat(", 'field_2'") }
    if (nK>2) { cat(", and so on") }
    cat(".\n")
    field_names <- paste0("field_", seq(nK))
  }
  if (any(duplicated(field_names))) { stop("All field names should be unique.") }


  # Third dimension, nD
  if (des_is_array) {
    des_nD <- vapply(des_dims, '[', 0, 3)
    if(!all(des_nD == des_nD[1])) stop('Location-specific design matrix should have the same number of design matrices across all sessions.')
    lapply(design, function(q){dimnames(q)[[3]]})
    #`nD` can differ across sessions only for per-location modeling.
    if (per_location_design) {
      design_names <- NULL
      if (!is.null(nD_expect)) {
        stopifnot(is.numeric(nD_expect))
        if (length(nD_expect)==1) { nD_expect <- rep(nD_expect, nS) }
        stopifnot(length(nD_expect)==1)
        for (ss in seq(nS)) {
          if (des_nD[ss] != nD_expect[ss]) {
            stop("The design array's third dimension is size ", des_nD[ss],
                 ", but for per-location modeling it should match",
                 "the number of locations, ", nD_expect[ss], ".")
          }
        }
      }
    } else {
      # for multiGLM
      if (length(unique(design_names)) > 1) {
        stop("Design names (third dim. of the design) should not differ across sessions.")
      }
      design_names <- design_names[[1]]
      # Set `design_names` for multiGLM, if not provided.
      if (is.null(design_names)) {
        cat("Setting design names to 'design_1'")
        if (nK>1) { cat(", 'design_2'") }
        if (nK>2) { cat(", and so on") }
        cat(".\n")
        design_names <- paste0("design_", seq(des_nD))
      }
    }
  } else {
    if (per_location_design) {
      stop("For per-locations design, the designs must be arrays, not matrices.")
    }
    des_nD <- 1
    design_names <- NULL
  }
  nD <- des_nD[1]; rm(des_nD)

  # Replace the `dimnames` of all the designs.
  des_dimnames <- list(vol=NULL, field=field_names)
  if (des_is_array) {
    if (per_location_design) {
      des_dimnames <- c(des_dimnames, list(des_vol=NULL))
    } else {
      des_dimnames <- c(des_dimnames, list(design=design_names))
    }
  }
  for (ss in seq(nS)) { dimnames(design[[ss]]) <- des_dimnames }

  # Scale design matrix. -------------------------------------------------------
  stopifnot(is_1(scale_design, "logical"))
  if (scale_design) { stop() }
  # if (!des_is_array) {
  #   design <- if(scale_design) {
  #     lapply(design, scale_design_mat)
  #   } else {
  #     lapply(design, scale, scale = FALSE)
  #   }
  # } else {
  #   for (dd in seq(nD)) {
  #     for (ss in seq(nS)) {
  #       if (scale_design) {
  #         design[[ss]][,,dd] <- scale_design_mat(design[[ss]][,,dd])
  #       } else {
  #         design[[ss]][,,dd] <- scale(design[[ss]][,,dd], scale=FALSE)
  #       }
  #     }
  #   }
  # }

  # Identify any missing fields in `valid_cols` for bookkeeping. ---------------
  valid_cols <- if (!des_is_array) {
    array(NA, c(nS, nK),
      dimnames = list(session=NULL, field=field_names)
    )
  } else if (!per_location_design) {
    array(NA, c(nS, nK, nD),
      dimnames = list(session=NULL, field=field_names, design=design_names)
    )
  } else {
    array(NA, c(nS, nK, nD),
      dimnames = list(session=NULL, field=field_names, loc=NULL)
    )
  }
  for (ss in seq(nS)) {
    if (!des_is_array) {
      valid_cols[ss,] <- colSums(abs(design[[ss]])) > 0
    } else {
      valid_cols[ss,,] <- apply(abs(design[[ss]]), 2, sum) > 0
    }
  }
  if (per_location_design) { valid_cols <- apply(valid_cols, 3, all) }

  # Return ---------------------------------------------------------------------
  list(
    design=design,
    nT=nT,
    nK=nK,
    nD=nD,
    field_names=field_names,
    design_names=design_names,
    valid_cols=valid_cols,
    per_location_design=per_location_design
  )
}


#' Format nuisance
#'
#' Format nuisance for \code{BayesGLM}, \code{fit_bayesglm},
#'  \code{multiGLM}, and \code{multiGLM_fun}.
#' @param nuisance The \code{nuisance} argument input. Will be formatted to a
#'  \code{nS}-length list.
#' @param nS_expect The expected number of sessions, if known.
#' @param nT_expect The expected number of timepoints, if known. For
#'  multi-session data this is a session-length vector.
#' @keywords internal
#' @return \code{nuisance}
BayesGLM_format_nuisance <- function(
  nuisance, nS_expect=NULL, nT_expect=NULL
  ){

  # Make `nuisance` a sessions-length list of matrices or arrays.
  if (!is.list(nuisance) || is.data.frame(nuisance)) { nuisance <- list(nuisance) }
  nS <- length(nuisance)
  if (!is.null(nS_expect)) {
    stopifnot(fMRItools::is_1(nS_expect, "numeric") && nS_expect==round(nS_expect))
    if (length(nuisance) != nS_expect) {
      stop("`nuisance` is length ", nS, ", but ", nS_expect, " sessions are expected.")
    }
  }

  for (ss in seq(nS)) {
    stopifnot(BayesGLM_is_valid_one_nuisance(nuisance[[ss]]))
    if (is.data.frame(nuisance[[ss]])) { nuisance[[ss]] <- as.matrix(nuisance[[ss]]) }
  }

  # First dimension, nT
  nuis_nT <- vapply(nuisance, nrow, 0)
  # Note that nuisances are allowed to have different numbers of timepoints,
  #   because BOLD data may be of unequal lengths.
  if (!is.null(nT_expect)) {
    stopifnot(is.numeric(nT_expect))
    if (length(nT_expect)==1) { nT_expect <- rep(nT_expect, nS) }
    stopifnot(length(nT_expect) == nS)
    for (ss in seq(nS)) {
      if (nuis_nT[ss] != nT_expect[ss]) {
        stop("The `nuisance` for session ", ss, "has ", nuis_nT[ss],
          " locations, but ", nT_expect[ss], " locations are expected.")
      }
    }
  }

  # [TO DO] detect and warn user about columns that look like spike regressors

  nuisance
}


#' Format scrub
#'
#' Format scrub for \code{BayesGLM}, \code{fit_bayesglm},
#'  \code{multiGLM}, and \code{multiGLM_fun}.
#' @param scrub The \code{scrub} argument input. Will be formatted to a
#'  \code{nS}-length list.
#' @param nS_expect The expected number of sessions, if known.
#' @param nT_expect The expected number of timepoints, if known. For
#'  multi-session data this is a session-length vector.
#' @keywords internal
#' @return \code{scrub}
BayesGLM_format_scrub <- function(
  scrub, nS_expect=NULL, nT_expect=NULL
  ){

  # Make `scrub` a sessions-length list of matrices or arrays.
  if (!is.list(scrub) || is.data.frame(scrub)) { scrub <- list(scrub) }
  nS <- length(scrub)
  if (!is.null(nS_expect)) {
    stopifnot(fMRItools::is_1(nS_expect, "numeric") && nS_expect==round(nS_expect))
    if (length(scrub) != nS_expect) {
      stop("`scrub` is length ", nS, ", but ", nS_expect, " sessions are expected.")
    }
  }

  for (ss in seq(nS)) {
    stopifnot(BayesGLM_is_valid_one_scrub(scrub[[ss]]))
    if (is.logical(scrub[[ss]])) {
      stopifnot(!is.null(nT_expect))
      if (length(scrub[[ss]]) != nT_expect[ss]) {
        stop("Logical vectors in `scrub` should be an indicator vector of ",
          "which volumes to remove. But for session ", ss, ", the length of ",
          "this vector, ", length(scrub[[ss]]), " does not match the length ",
          "expected from the data ,", nT_expect[ss], ".")
      }
      # Convert from indicator vector to spikes
      spikes <- matrix(0, nrow=nT_expect[ss], ncol=sum(scrub[[ss]]))
      spikes[seq(0, (sum(scrub[[ss]])-1))*nT_expect[ss] + which(scrub[[ss]])] <- 1
      scrub[[ss]] <- spikes
    }
    if (is.data.frame(scrub[[ss]])) { scrub[[ss]] <- as.matrix(scrub[[ss]]) }
  }

  # First dimension, nT
  nuis_nT <- vapply(scrub, nrow, 0)
  # Note that scrubs are allowed to have different numbers of timepoints,
  #   because BOLD data may be of unequal lengths.
  if (!is.null(nT_expect)) {
    stopifnot(is.numeric(nT_expect))
    if (length(nT_expect)==1) { nT_expect <- rep(nT_expect, nS) }
    stopifnot(length(nT_expect) == nS)
    for (ss in seq(nS)) {
      if (nuis_nT[ss] != nT_expect[ss]) {
        stop("The `scrub` for session ", ss, "has ", nuis_nT[ss],
          " locations, but ", nT_expect[ss], " locations are expected.")
      }
    }
  }

  scrub
}
