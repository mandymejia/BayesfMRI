#' Remove from design
#'
#' Remove field from result of \code{\link{make_design}}.
#'
#' @param design A \code{"BayesfMRI_design"} object from \code{\link{make_design}}.
#' @param rm_name,rm_pattern The fields to remove. Either provide the names to
#'  \code{rm_name} or the patterns to match with \code{\link[base]{grepl}} to
#'  \code{rm_pattern}.
#'
#'  By default, \code{rm_pattern} specifies the removal of fields ending in
#'  "_dHRF" or "_ddHRF". If \code{design} was constructed from \code{onsets},
#'  then these are the HRF derivatives.
#' @return A list of two: the modified \code{"BayesfMRI_design"} object, and
#'  the removed columns.
remove_from_design <- function(
  design, rm_name=NULL, rm_pattern=".*_dHRF|.*_ddHRF"){

  stopifnot(inherits(design, "BfMRI_design"))
  stopifnot(xor(is.null(rm_name), is.null(rm_pattern)))

  to_rm <- if (!is.null(rm_name)) {
    dimnames(design$design[[1]])[[2]] == rm
  } else if (!is.null(rm_pattern)) {
    grepl(rm_pattern, dimnames(design$design[[1]])[[2]])
  } else { stop() }

  removed <- lapply(design$design, function(q){ q[,to_rm,drop=FALSE] })

  nDD <- length(dim(design$design[[1]]))
  for (ss in seq(length(design$design))) {
    design$design[[ss]] <- if (nDD==2) {
      design$design[[ss]][,!to_rm,drop=FALSE]
    } else {
      design$design[[ss]][,!to_rm,,drop=FALSE]
    }
  }

  design$dims["fields","count"] <- sum(!to_rm)
  design$valid_cols <- if (nDD==2) {
    design$valid_cols[,!to_rm,drop=FALSE]
  } else {
    design$valid_cols[,!to_rm,,drop=FALSE]
  }

  list(design=design, removed=removed)
}
