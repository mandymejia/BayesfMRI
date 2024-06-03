#' Format BayesGLM_fun results into \code{"xifti"} objects
#'
#' Format beta estimates of the BayesGLM_fun results into list of \code{"xifti"}
#'  objects for each session.
#'
#' @param BGLMs The list of \code{"BGLM"} objects, from \code{BayesGLM_fun}.
#' @param do,spatial,submeta,session_names,field_names See \code{BayesGLM}.
#' @param method \code{"classical"} or \code{"Bayesian"}
#' @return The list of \code{"xifti"} objects
#' @keywords internal
#'
BayesGLM_format_cifti <- function(
  BGLMs, do,
  spatial, submeta,
  session_names, field_names,
  method=c("classical", "Bayesian")){

  nS <- length(session_names)
  method <- match.arg(method, c("classical", "Bayesian"))

  result_xii <- setNames(vector("list", nS), session_names)
  datL <- datR <- datSub <- NULL
  for (ss in seq(nS)) {
    if (do$left) {
      datL <- switch(method,
        classical = BGLMs$cortexL$result_classical[[ss]]$estimates,
        Bayesian = BGLMs$cortexL$field_estimates[[ss]]
      )
      # Update mwall b/c `mask2` in `BayesGLM_fun` can change the medial wall.
      mwallL <- BGLMs$cortexL$spatial$mask
      colnames(datL) <- NULL
    }
    if (do$right) {
      datR <- switch(method,
        classical = BGLMs$cortexR$result_classical[[ss]]$estimates,
        Bayesian = BGLMs$cortexR$field_estimates[[ss]]
      )
      mwallR <- BGLMs$cortexR$spatial$mask
      colnames(datR) <- NULL
    }
    if (do$sub) {
      datSub <- switch(method,
        classical = BGLMs$subcort$result_classical[[ss]]$estimates,
        Bayesian = BGLMs$subcort$field_estimates[[ss]]
      )
      colnames(datSub) <- NULL
    }
    result_xii[[ss]] <- as.xifti(
      cortexL = datL,
      cortexL_mwall = if (do$left) { mwallL } else { NULL },
      cortexR = datR,
      cortexR_mwall = if (do$right) { mwallR } else { NULL },
      subcortVol = datSub,
      subcortLabs = submeta$labels,
      subcortMask = submeta$mask
    )
    result_xii[[ss]]$meta$subcort$trans_mat <- submeta$trans_mat
    result_xii[[ss]]$meta$subcort$trans_units <- submeta$trans_units
    result_xii[[ss]]$meta$cifti$names <- field_names
  }

  result_xii
}
