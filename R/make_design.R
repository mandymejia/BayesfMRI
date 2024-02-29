#' Make design for BayesGLM
#'
#' Make the design matrix/array, or matrices/arrays, for
#'  \code{\link{BayesGLM_cifti}} or \code{\link{CompareGLM_cifti}},
#'  from information about the task design.
#'
#' Provide one of: (1) \code{design}; (2) \code{design_per_location}; (3)
#'  \code{design_compare}; or (4) all of \code{onsets}, \code{nTime}, \code{TR},
#'  and \code{dHRF}.
#'
#' @param design \eqn{T \times K} task design matrix. Each column
#'  represents the expected BOLD response for a field, a convolution of
#'  a hemodynamic response function (HRF) (or its derivative) with a task
#'  stimulus. Field names must be given by the column names.
#'
#'  If there are multiple sessions to model, this argument should be a list of
#'  such matrices. The list names should be the session names. Column names must
#'  match across sessions; if any field is missing from a given session, include
#'  a column of zeros in its place.
#'
#'  Tasks convolved with an HRF derivative may be included in \code{design} to
#'  model spatially; alternatively, provide these in \code{BayesGLM(_cifti)} as
#'  columns in \code{nuisance} to treat as nuisance signals. Note that INLA
#'  computation times increase if \eqn{K > 5}.
#' @param design_per_location Similar to \code{design}, but instead of a matrix,
#'  an \eqn{T \times K \times V} array, to use a unique design matrix for each
#'  of the \eqn{V} locations.
#' @param design_compare Similar to \code{design}, but instead of a matrix,
#'  an \eqn{T \times K \times D} array, to compare \eqn{D} design matrix choices.
#'
#   [TO DO]: for CompareGLM: Allow differing numbers of regressors across
#   \eqn{D} models, pad with NAs or 0.
#' @param onsets Task stimulus information, from which a design matrix (or
#'  matrices, for multi-session data) will be constructed. This is a list where
#'  each entry represents a task as a matrix of onsets (first column) and
#'  durations (second column) for each stimuli (each row) of the task, in
#'  seconds. List names should be the task names. \code{nTime} and \code{TR} are
#'  required with \code{onsets}.
#'
#'  For multi-session modeling, this argument should be a list of such lists.
#'  The list names should be the session names. Task names must match across
#'  sessions; if any task is missing from any session, include a list element
#'  equal to \code{NA} in its place.
#'
#'  An example of a properly-formatted single-session \code{onsets} is:
#'  \code{on_s1 <- list(taskA=cbind(on=c(1,9,17), dr=rep(1,3)), taskB=cbind(on=c(3,27), dr=rep(5,2)))}.
#'  In this session, there are two tasks: the first has three 1s-long stimuli,
#'  while the second has two 5s-long stimuli. A multi-session \code{onsets}
#'  would look like \code{on_multi <- list(s1=on_s1, s2=on_s2)},
#'  with \code{on_s2} formatted similarly to \code{on_s1}.
#' @param nTime,TR These are required if and only if \code{onsets} is provided.
#'
#'  \code{nTime} is the number of timepoints (volumes) in the task fMRI data.
#'  This can either be a single positive integer, or a vector of positive
#'  integers if there are multiple sessions each with different lengths.
#'
#'  \code{TR} is the temporal resolution of the data, in seconds.
#' @param dHRF Controls the extent of HRF derivatives modeling. Applicable if
#'  and only if \code{onsets} is provided.
#'
#'  Set to \code{0} to only model the main HRF regressor, and not include its
#'  derivatives; set to \code{1} (default) to model the temporal derivative too;
#'  or, set to \code{2} to model both the temporal and dispersion derivatives.
#'  If \code{dHRF==0}, there is one design column (field) per task. If
#'  \code{dHRF==1}, there are two fields per task. And if \code{dHRF==2}, there
#'  are three fields per task.
#'
#'  If there are several tasks and \code{dHRF>0}, the total number of design
#'  matrix columns may exceed five, which may require large computation times
#'  with INLA. The analysis can be adjusted by modeling the derivatives as
#'  nuisance signals rather than as fields. To do so, move the corresponding
#'  columns from the design matrix to the \code{nuisance} argument for
#'  \code{BayesGLM_cifti}.
#' @param scale_design Scale the design matrix by dividing each column by its
#'  maximum and then subtracting the mean? Default: \code{TRUE}. If
#'  \code{FALSE}, the design matrix is centered but not scaled.
#' @param session_names (Optional) Session names can be provided by the list
#'  names of \code{design}, \code{design_compare}, \code{design_per_location},
#'  or \code{onsets}; alternatively, session names can be provided by this
#'  argument.
#' @param ... Additional arguments to \code{\link{make_HRFs}}. Applicable if
#'  and only if \code{onsets} is provided.
#'
#' @return A \code{"BayesfMRI_design"} object: a list with elements
#'  \describe{
#'    \item{design}{The design matrix/array, or list of these for multi-session analysis.}
#'    \item{HRF_info}{Additional HRF modeling results, if \code{onsets} was provided.}
#'    \item{design_type}{The source of the design.}
#'    \item{dims}{The numbers of sessions, volumes, fields, and designs.}
#'    \item{field_names}{The names of each field.}
#'    \item{valid_cols}{Indicates sessions with tasks that had no stimuli.}
#' }
#' @export
make_design <- function(
  design=NULL, design_compare=NULL, design_per_location=NULL,
  onsets=NULL,
  nTime=NULL, TR=NULL, dHRF=c(1, 0, 2),
  scale_design = TRUE,
  session_names = NULL,
  ...
){

  # Identify how the task info has been provided. ------------------------------
  arg_provided <- vapply(
    list(design, design_compare, design_per_location, onsets),
    function(q){!is.null(q)}, FALSE
  )
  if (sum(arg_provided)!=1) { stop(
      "Provided exactly one of: `design`; `design_compare`; ",
      "`design_per_location`; `onsets` (with nTime and TR)."
  )}
  design_type <- c("design", "compare", "per_loc", "onsets")[arg_provided]

  # Merge the `design` args for code simplicity.
  if (design_type=="compare") {
    design <- design_compare; rm(design_compare)
  } else if (design_type=="per_loc") {
    design <- design_per_location; rm(design_per_location)
  }

  # Get `design` and `nS`. -----------------------------------------------------
  ### Format provided `design`... ----------------------------------------------
  if (!is.null(design)) {

    # Make `design` a named list of length `nS`.
    if (!is.list(design)) {
      design <- list(single_session=design)
    } else {
      if (is.null(names(design))) {
        names(design) <- paste0("session_", seq(length(design)))
      }
    }
    nS <- length(design)
    if (nS > 1 && design_type=="compare") { stop(
      "Multi-session modeling is not yet supported for `CompareGLM`."
    )}

    # `is.TRUE(design_is_array)` for CompareGLM or per-location modeling.
    nDD <- length(dim(design[[1]]))
    design_is_array <- nDD == 3
    if (nDD==1) {
      stop("The design should not be a 1D vector.")
    } else if (nDD==2) {
      if (design_type!="design") {
        stop("A 3D array is expected for `design_compare` or `design_per_location`.")
      }
    } else if (nDD==3) {
      if (design_type=="design") {
        stop(
          "`design` should be a matrix, not a 3-dimensional array. For design ",
          "comparison or per-location modeling, instead use `design_compare` ",
          "or `design_per_location`, respectively."
        )
      }
    } else {
      stop("The design has too many dimensions.")
    }

    # Checks that `design` is a numeric matrix or array.
    # For multi-session modeling, check that dimensions are consistent
    #   across sessions.
    stopifnot(all(vapply(design, is.numeric, FALSE)))
    stopifnot(length(unique(lapply(design, dim)))==1)
    if (!design_is_array) {
      stopifnot(all(vapply(design, is.matrix.or.df, FALSE)))
      stopifnot(all(vapply(design, ncol, 0)>0))
    } else {
      stopifnot(all(vapply(design, is.array, FALSE)))
      stopifnot(all(vapply(design, function(q){length(dim(q))==3}, FALSE)))
      stopifnot(all(vapply(design, function(q){dim(q)[2]}, 0)>0))
    }

    # Rename the dimnames of `design`. (Ignore whatever the user used.)
    dimnames_names <- c("vol", "field")
    if (design_type=="compare") { dimnames_names <- c(dimnames_names, "design") }
    if (design_type=="per_loc") { dimnames_names <- c(dimnames_names, "loc") }
    for (ss in seq(nS)) {
      if (is.null(dimnames(design[[ss]]))) {
        dimnames(design[[ss]]) <- vector("list", nDD)
      }
      names(dimnames(design[[ss]])) <- dimnames_names
    }

    # Second dim, `field`: provide if absent; check same across sess if present.
    field_names_all <- lapply(design, function(q){dimnames(q)[[2]]})
    if (length(unique(field_names_all)) > 1) {
      stop("Field names (second dim. of the design) should not differ across sessions.")
    }
    if (is.null(field_names_all[[1]])) {
      for (ss in seq(nS)) {
        dimnames(design[[ss]])[[2]] <- paste0("field_", seq(dim(design[[1]])[2]))
      }
    }
    rm(field_names_all)

    # Third dim, if applicable. Handled similarly to second dim.
    if (design_is_array) {
      design_names_all <- lapply(design, function(q){dimnames(q)[[3]]})
      if (length(unique(design_names_all)) > 1) {
        stop("Design names (third dim. of the design) should not differ across sessions.")
      }
      if (is.null(design_names_all[[1]])) {
        if (design_type=="compare") {
          for (ss in seq(nS)) {
            dimnames(design[[ss]])[[3]] <- paste0("design_", seq(dim(design[[1]])[3]))
          }
        } else if (design_type=="per_loc") {
          for (ss in seq(nS)) {
            dimnames(design[[ss]])[3] <- list(loc=NULL)
          }
        } else { stop() }
      }
      rm(design_names_all)
    }

    # Set `HRF_info` to `NULL` (only applicable for `onsets`)
    HRF_info <- NULL # only for `onsets`

  ### ... Or, format `onsets` and construct design from it. --------------------
  } else if (!is.null(onsets)) {
    design_is_array <- FALSE

    # Make `onsets` a named list of length `nS`.
    stopifnot(is.list(onsets))
    if (is.data.frame(onsets[[1]])) {
      onsets <- list(single_session=onsets)
    } else {
      if (is.null(names(onsets))) {
        names(onsets) <- paste0("session_", seq(length(onsets)))
      }
    }
    nS <- length(onsets)

    for (ss in seq(nS)) {
      stopifnot(all(vapply(onsets[[ss]], is_onsets, FALSE)))
    }

    # Construct design.
    x <- make_HRFs(onsets, nTime, TR, dHRF, ...)
    design <- lapply(x, '[[', "design")
    HRF_info <- lapply(x, '[', c("HRF", "stimulus")) # "FIR"
    rm(x)

  } else { stop }

  # Collect design matrix dims info. -------------------------------------------
  nTime <- dim(design[[1]])[1]

  if (!is.null(session_names)) {
    stopifnot(is.character(session_names))
    stopifnot(length(session_names) == nS)
    stopifnot(!any(duplicated(session_names)))
    names(design) <- session_names
    if (design_type == "onsets") { names(HRF_info) <- session_names }
  } else {
    session_names <- names(design)
  }

  # Note each `design` is TxK, or TxKxD (compare), or TxKxV (per-location).

  # Get the total number of fields and `field_names`.
  nK <- dim(design[[1]])[2]
  field_names <- dimnames(design[[1]])[[2]]

  # Get total number of designs and `design_names`.
  nD <- if (design_is_array) { dim(design[[1]])[3] } else { 1 }
  if (design_is_array) {
    design_names <- dimnames(design[[1]])[[3]]
  } else {
    design_names <- "single_design"
  }

  # Scale design matrix. -------------------------------------------------------
  stopifnot(is_1(scale_design, "logical"))
  if (!design_is_array) {
    design <- if(scale_design) {
      lapply(design, scale_design_mat)
    } else {
      lapply(design, scale, scale = FALSE)
    }
  } else {
    for (dd in seq(nD)) {
      for (ss in seq(nS)) {
        if (scale_design) {
          design[[ss]][,,dd] <- scale_design_mat(design[[ss]][,,dd])
        } else {
          design[[ss]][,,dd] <- scale(design[[ss]][,,dd], scale=FALSE)
        }
      }
    }
  }

  # Make `dims`: table of dimensions, for output. ------------------------------
  dims_nT <- if (length(unique(nTime))==1) {
    nTime
  } else {
    paste0("variable (", min(nTime), "-", max(nTime), ")")
  }
  dims = data.frame(
    row.names = c("sessions", "volumes", "fields", "design_matrices"),
    count = c(nS, dims_nT, nK, nD)
  )

  # Identify any missing fields in `valid_cols` for bookkeeping. ---------------
  valid_cols <- if (design_type %in% c("design", "onsets")) {
    array(NA, c(nS, nK),
      dimnames = list(session=session_names, field=field_names)
    )
  } else if (design_type == "compare") {
    array(NA, c(nS, nK, nD),
      dimnames = list(session=session_names, field=field_names, design=design_names)
    )
  } else if (design_type == "per_loc") {
    array(NA, c(nS, nK, nD),
      dimnames = list(session=session_names, field=field_names, loc=NULL)
    )
  } else { stop() }
  for (ss in seq(nS)) {
    if (!design_is_array) {
      valid_cols[ss,] <- colSums(abs(design[[ss]])) > 0
    } else {
      valid_cols[ss,,] <- apply(abs(design[[ss]]), 2, sum) > 0
    }
  }

  out <- list(
    design=design, HRF_info=HRF_info,
    design_type=design_type,
    dims=dims, field_names=field_names,
    valid_cols=valid_cols
  )
  class(out) <- "BfMRI_design"
  out
}

#' Is this a valid entry in `onsets`?
#'
#' Is this valid data for a single task's onsets? Expects a data.frame or
#'  numeric matrix with two numeric columns, onsets and durations, and at least
#'  one row.
#'
#' @param x The putative onsets matrix or data frame
#' @keywords internal
#' @return Length-one logical vector.
#'
is_onsets <- function(x){

  # First check if onsets is NA, which can be the case in multi-session analysis
  #   where not all fields are present in all sessions.
  if(length(x)==1 && is.na(x)) { return(TRUE) }

  is_nummat <- is.numeric(x) & is.matrix(x)
  is_df <- is.data.frame(x) && all(vapply(x, class, "") == "numeric")
  if (!(is_nummat || is_df)) { warning("The onsets are not a numeric matrix or data.frame."); return(FALSE) }

  if (nrow(x)<1) { warning("The onsets must have at least one row."); return(FALSE) }

  if (ncol(x) != 2) { warning("The onsets should have two columns, `onset` and `duration`."); return(FALSE) }
  if (!all(colnames(x) == c("onset", "duration"))) {
    warning("The onsets should have two columns, `onset` and `duration`."); return(FALSE)
  }

  TRUE
}
