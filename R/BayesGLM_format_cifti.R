#' Format fit_bayesglm results into \code{"xifti"} objects
#'
#' Format beta estimates or RSS of the \code{fit_bayesglm} results into list of
#'  \code{"xifti"} objects for each session.
#'
#' @param BGLMs The list of \code{"BGLM"} objects, from \code{fit_bayesglm}.
#' @param do,spatial,submeta,session_names,field_names See \code{BayesGLM}.
#' @param method \code{"classical"} or \code{"Bayesian"}
#' @return The list of \code{"xifti"} objects
#' @importFrom ciftiTools as.xifti
#' @keywords internal
#'
BayesGLM_format_cifti <- function(
  BGLMs,
  do, spatial, submeta, session_names, field_names,
  method=c("classical", "Bayesian")){

  nS <- length(session_names)
  method <- match.arg(method, c("classical", "Bayesian"))

  out <- list(estimates=NULL, RSS=NULL)
  for (what in c("estimates", "RSS")) {
    dat_name <- paste0(method, "_", what)
    result_xii <- setNames(vector("list", nS), session_names)
    datL <- datR <- datSub <- NULL
    for (ss in seq(nS)) {
      if (do$left) {
        datL <- switch(dat_name,
          classical_estimates = BGLMs$cortexL$result_classical[[ss]]$estimates,
          Bayesian_estimates = BGLMs$cortexL$field_estimates[[ss]],
          classical_RSS = BGLMs$cortexL$result_classical[[ss]]$RSS,
          Bayesian_RSS = BGLMs$cortexL$RSS[[ss]],
        )
        mwallL <- BGLMs$cortexL$spatial$maskIn
        colnames(datL) <- NULL
      }
      if (do$right) {
        datR <- switch(dat_name,
          classical_estimates = BGLMs$cortexR$result_classical[[ss]]$estimates,
          Bayesian_estimates = BGLMs$cortexR$field_estimates[[ss]],
          classical_RSS = BGLMs$cortexR$result_classical[[ss]]$RSS,
          Bayesian_RSS = BGLMs$cortexR$RSS[[ss]]
        )
        mwallR <- BGLMs$cortexR$spatial$maskIn
        colnames(datR) <- NULL
      }
      if (do$sub) {
        datSub <- switch(dat_name,
          classical_estimates = BGLMs$subcort$result_classical[[ss]]$estimates,
          Bayesian_estimates = BGLMs$subcort$field_estimates[[ss]],
          classical_RSS = BGLMs$subcort$result_classical[[ss]]$RSS,
          Bayesian_RSS = BGLMs$subcort$RSS[[ss]]
        )
        subMask <- BGLMs$subcort$spatial$maskIn
        subLabs <- BGLMs$subcort$spatial$labels
        colnames(datSub) <- NULL
      }
      result_xii[[ss]] <- ciftiTools::as.xifti(
        cortexL = datL,
        cortexL_mwall = if (do$left) { mwallL } else { NULL },
        cortexR = datR,
        cortexR_mwall = if (do$right) { mwallR } else { NULL },
        subcortVol = if (do$sub) { as.matrix(datSub) } else { NULL },
        subcortLabs = if (do$sub) { subLabs } else { NULL },
        subcortMask = if (do$sub) { subMask } else { NULL },
      )
      result_xii[[ss]]$meta$subcort$trans_mat <- submeta$trans_mat
      result_xii[[ss]]$meta$subcort$trans_units <- submeta$trans_units
      result_xii[[ss]]$meta$cifti$names <- field_names
    }

    out[[what]] <- result_xii
  }

  out
}
