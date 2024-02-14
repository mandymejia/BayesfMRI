#' Make design for BayesGLM
#'
#' Make the design matrix/array, or matrices/arrays, for \code{\link{BayesGLM}}
#'  or \code{\link{BayesGLM_cifti}}, from information about the tasks.
#'
#' Provide either (1) \code{design}, (2) \code{design_multiple}, or (3) all of
#'  \code{onsets}, \code{nTime}, \code{TR}, and \code{dHRF}.
#'
#' @param design \eqn{T \times K} task design matrix. Each column
#'  represents the expected BOLD response for a field, a convolution of
#'  a hemodynamic response function (HRF) (or its derivative) with a task stimulus.
#'  Field names must be given by the column names.
#'
#'  For multi-session modeling, this argument should be a list of such matrices,
#'  and column names must match across sessions. Session names must be given by
#'  the list names.
#'
#'  If any field is missing from any session, include a column of zeros in its
#'  place.
#'
#'  Tasks convolved with an HRF derivative may be included in
#'  \code{design} to model them spatially; alternatively, provide these in
#'  \code{BayesGLM(_cifti)} as columns in \code{nuisance} to treat as nuisance
#'  signals. Note that INLA computation times increase if the design matrix has
#'  more than five columns.
#' @param design_multiple Same as \code{design}, but rather than a matrix, use
#'  a \eqn{T \times K \times D} array of \eqn{D} different design matrices for
#'  model comparison.
# TO DO: Allow differing numbers of regressors across D models, pad with NAs or 0.
#' @param onsets A list where each element represents a
#'  task, and the value of each element is a matrix of onsets (first column) and
#'  durations (second column) for each stimuli (each row) of the corresponding
#'  task. The units of both columns is seconds. Task names must be given by
#'  the list names. \code{nTime} and \code{TR} are required if \code{onsets} is
#'  provided.
#'
#'  For multi-session modeling, this argument should be a list of such lists,
#'  and the task names must match across sessions. Session names must be given by
#'  the list names.
#'
#'  If any task is missing from any session, include a list element equal to
#'  \code{NA} in its place.
#' @param nTime,TR These are used if and only if \code{onsets} is provided.
#'
#'  \code{nTime} is the number of timepoints (volumes) in the task fMRI data.
#'  This can either be a single positive integer, or a vector of positive
#'  integers if there are multiple sessions each with different numbers of
#'  timepoints. There is no default value, so \code{nTime} is required if
#'  \code{onsets} is provided.
#'
#'  \code{TR} is the temporal resolution of the data, in seconds. There is no
#'  default value, so \code{TR} is required if \code{onsets} is provided.
#' @param scale_design Scale the design matrix by dividing each column by its
#'  maximum and then subtracting the mean? Default: \code{TRUE}. If
#'  \code{FALSE}, the design matrix is centered but not scaled.
#'
#' @return A list of two: the design matrix/arrays or list of them; and, extra
#'  results if \code{onsets} was provided.
#' @export
make_design <- function(
  design=NULL,
  design_multiple=NULL,
  onsets=NULL,
  nTime=NULL, TR=NULL,
  scale_design = TRUE
){

  # Identify how the tasks have been provided. ---------------------------------
  # Check: exactly one of `design`, `design_multiple`, or `onsets` (with `TR`).
  # Get: number of sessions, `nS`; `session_names`.
  # Get either:
  #   number of tasks, `nJ`; `task_names`.      (N/A for `design`, `design_multiple`)
  #   number of fields, `nK`; `field_names`.    (determined later for `onsets`)
  ok_design_args <- TRUE
  is.matrix.or.df <- function(q){is.matrix(q)||is.data.frame(q)}
  #   `design`
  if (!is.null(design)) {
    if (!is.null(design_multiple)) { ok_design_args <- FALSE }
    if (!is.null(onsets)) { ok_design_args <- FALSE }
    if (is.list(design)) {
      nS <- length(design)
      stopifnot(all(vapply(design, is.numeric, FALSE)))
      stopifnot(all(vapply(design, is.matrix.or.df, FALSE)))
      stopifnot(all(vapply(design, ncol, 0)>0))
      stopifnot(length(unique(lapply(design, dim)))==1)
    } else {
      nS <- 1
      stopifnot(is.numeric(design) && is.matrix.or.df(design))
      design <- list(single_session=design)
    }
    field_names <- colnames(design[[1]])
    nK <- length(field_names)
    session_names <- names(design)
  #   `design_multiple`
  } else if (!is.null(design_multiple)) {
    if (!is.null(design)) { ok_design_args <- FALSE }
    if (!is.null(onsets)) { ok_design_args <- FALSE }
    if (is.list(design_multiple)) {
      nS <- length(design_multiple)
      stopifnot(all(vapply(design_multiple, is.numeric, FALSE)))
      stopifnot(all(vapply(design_multiple, is.array, FALSE)))
      stopifnot(all(vapply(design_multiple, function(q){dim(q)[2]}, 0)>0))
      stopifnot(all(vapply(design_multiple, function(q){length(dim(q))==3}, FALSE)))
      stopifnot(length(unique(lapply(design_multiple, dim)))==1)
    } else {
      nS <- 1
      stopifnot(is.numeric(design_multiple) && is.array(design_multiple))
      design_multiple <- list(single_session=design_multiple)
    }
    field_names <- dimnames(design_multiple[[1]])[[2]]
    nK <- length(field_names)
    session_names <- names(design_multiple)
  #   `onsets`
  } else if (!is.null(onsets)) {
    if (!is.null(design)) { ok_design_args <- FALSE }
    if (!is.null(design_multiple)) { ok_design_args <- FALSE }
    if (is.list(onsets[[1]])) {
      nS <- length(onsets)
      stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(r){all(apply(r, 2, is.numeric))}, FALSE))}, FALSE)))
      stopifnot(all(vapply(onsets, function(q){all(vapply(q, is.matrix.or.df, FALSE))}, FALSE)))
      stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(v){ncol(v)==2}, FALSE))}, FALSE)))
      stopifnot(length(unique(lapply(onsets, length)))==1)
      stopifnot(length(unique(lapply(onsets, names)))==1)
    } else {
      nS <- 1
      stopifnot(is.list(onsets))
      stopifnot(all(vapply(onsets, function(v){ncol(v)==2}, FALSE)))
      onsets <- list(single_session=onsets)
    }
    task_names <- names(onsets[[1]])
    # nJ <- length(task_names) # never used.
    # nK will depend on `dHRF`, `dHRF_as` in `BayesGLM_cifti`.
    session_names <- names(onsets)
    if (is.null(nTime)) { stop("`nTime` is required if `onsets` is provided.") }
    stopifnot(is.numeric(nTime))
    stopifnot(all(nTime > 0))
    stopifnot(all(nTime == round(nTime)))
    if (length(nTime)==1) {
      nTime <- rep(nTime, nS)
    } else {
      stopifnot(length(nTime)==nS)
    }
    if (is.null(TR)) { stop("`TR` is required if `onsets` is provided.") }
    stopifnot(is.numeric(TR))
    stopifnot(TR > 0)
  #   Error if none of the three are provided
  } else {
    ok_design_args <- FALSE
  }
  if (!ok_design_args) {
    stop("Exactly one must be provided: `design`, `design_multiple`, or `onsets`. Please consult the help documentation for `make_design`.")
  }

  # Check the other arguments. -------------------------------------------------
  stopifnot(is_1(scale_design, "logical"))

  # Construct design matrix for `onsets`. --------------------------------------
  if (!is.null(onsets)) {
    x <- make_design_from_onsets(onsets, nTime, TR)
    design <- x$design
    onsets_misc <- x$misc
  } else { onsets_misc <- NULL }

  ## Scale design matrix. ------------------------------------------------------
  if (!is.null(design)) {
    design <- if(scale_design) {
      sapply(design, scale_design_mat, simplify = FALSE)
    } else {
      sapply(design, scale, scale = FALSE, simplify = FALSE)
    }
  }

  if (!is.null(design_multiple)) {
    for (dd in seq(dim(design_multiple[[1]])[3])) {
      for (ss in seq(length(design))) {
        if (scale_design) {
          design[[ss]][,,dd] <- scale_design_mat(design[[ss]][,,dd])
        } else {
          design[[ss]][,,dd] <- scale(design[[ss]][,,dd], scale=FALSE)
        }
      }
    }
  }

  out <- list(design=design, onsets_misc=onsets_misc)
  class(out) <- "BfMRI_design"
  out
}

#' Make design from onsets
#'
#' Helper function to \code{\link{make_design}} for constructing a design matrix
#'  from \code{onsets}, \code{nTime}, and \code{TR}.
#' @param onsets,nTime,TR See \code{\link{make_design}}.
#' @return A list of two: the list of design matrices; and, extra results.
#' @keywords internal
make_design_from_onsets <- function(onsets, nTime, TR){

  # Prep. ----------------------------------------------------------------------
  session_names <- names(onsets)
  nS <- length(onsets)
  task_names <- names(onsets[[1]])
  nJ <- length(task_names)
  # [NOTE] Argument checks are in `make_design`.

  # Compute over sessions. -----------------------------------------------------
  # Get: `stimulus`, `HRF`, `HRFs` & derivatives, `FIR`, `design`, `design_FIR`.

  # Do `make_HRFs` for each session.
  x <- lapply(seq(nS), function(ss){
    make_HRFs(onsets[[ss]], nT=nTime[ss], TR=TR, dHRF=2)
  })

  # Aggregate results.
  design <- lapply(x, '[[', "design") # TxJ
  stimulus <- lapply(x, '[[', "stimulus") #TxJ
  HRF <- x[[1]]$HRF # (upsampled T)x4
  HRFs <- lapply(x, '[[', "HRF_convolved") #TxKx3 array

  HRFs_d <- if (is.null(HRFs[[1]])) { NULL } else {
    lapply(HRFs, function(q){
      q <- as.matrix(q[,,2], ncol=nJ)
      colnames(q) <- paste0(task_names, "_dHRF")
      q
    })
  }
  HRFs_dd <- if (is.null(HRFs[[1]])) { NULL } else {
    lapply(HRFs, function(q){
      q <- as.matrix(q[,,3], ncol=nJ)
      colnames(q) <- paste0(task_names, "_ddHRF")
      q
    })
  }
  HRFs <- if (is.null(HRFs[[1]])) { NULL } else {
    lapply(HRFs, function(q){
      q <- as.matrix(q[,,1], ncol=nJ)
      colnames(q) <- paste0(task_names)
      q
    })
  }

  FIR <- NULL # lapply(x, '[[', "FIR") #T x K x nFIR -- just return this for now, don't use in modeling
  design_FIR <- if (is.null(FIR)) { NULL } else {
    lapply(FIR, function(q){
      # [NOTE] Untested
      q <- as.matrix(q, ncol=nJ)
      colnames(q) <- c(outer(paste0(task_names, "_FIR"), seq(nJ), paste0))
      q
    })
  }

  list(
    design = design,
    misc = list(
      stimulus=stimulus,
      HRF=HRF,
      HRFs=HRFs,
      HRFs_d=HRFs_d,
      HRFs_dd=HRFs_dd,
      FIR=FIR,
      design_FIR=design_FIR
    )
  )
}
