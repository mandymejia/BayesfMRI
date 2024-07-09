#' Intersection mask for BayesGLM or activations result
#'
#' @param x The list of \code{"fit_bglm"}, \code{"BGLM"}, or \code{"act_BGLM"} objects.
#' @return The intersections masks
#'
#' @keywords internal
intersect_mask <- function(x) {
  what <- if (inherits(x[[1]], "fit_bglm")) {
    "fit_bglm"
  } else if (inherits(x[[1]]$BGLMs[[1]], "fit_bglm")) {
    "BGLM"
  } else if (inherits(x[[1]], "act_BGLM")) {
    "act_BGLM"
  } else { stop("`x` must be a list of 'BGLM' or 'act_BGLM' objects.") }
  if (!all(vapply(x, inherits, FALSE, what))) {
    if (!all(vapply(x, function(q){inherits(q$BGLMs[[1]], "fit_bglm")}, FALSE))) {
      stop("`x` must be a list of 'BGLM' or 'act_BGLM' objects.")
    }
  }

  brainstructures <- switch(what,
    fit_bglm = "unknown",
    BGLM = names(x[[1]]$BGLMs),
    act_BGLM = names(x[[1]]$spatial)
  )
  nB <- length(brainstructures)
  # [TO DO] check that brainstructures match across sessions.

  # Get intersection mask for each brainstructure.
  Masks <- setNames(vector("list", nB), brainstructures)
  for (bb in seq(nB)) {
    bs <- brainstructures[bb]
    spatial_type_bb <- switch(bs,
      cortexL = "vertex",
      cortexR = "vertex",
      subcort = "voxel",
      unknown = x[[bb]]$spatial$spatial_type # BGLM case
    )

    masks <- if (what == "fit_bglm") {
      ## nB == 1
      do.call(rbind, lapply(x, function(q){
        as.logical(q$spatial$maskMdat)
      }))
    } else if (what == "BGLM") {
      do.call(rbind, lapply(x, function(q){
        as.logical(q$BGLMs[[bb]]$spatial$maskMdat)
      }))
    } else if (what == "act_BGLM") {
      do.call(rbind, lapply(x, function(q){
        as.logical(q$spatial[[bb]]$maskMdat)
      }))
    } else { stop() }
    browser()
    Masks[[bb]] <- apply(masks, 2, all)
  }

  Masks
}
