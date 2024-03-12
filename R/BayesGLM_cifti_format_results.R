#' Format BayesGLM results into \code{"xifti"} objects
#'
#' Format beta estimates of the BayesGLM results into list of \code{"xifti"}
#'  objects for each session.
#'
#' @param BayesGLM_results The list of BayesGLM results
#' @param do,spatial,submeta,session_names,field_names See \code{BayesGLM_cifti}.
#' @param method \code{"classical"} or \code{"Bayesian"}
#' @return The list of \code{"xifti"} objects
#' @keywords internal
#'
BayesGLM_cifti_format_results <- function(
  BayesGLM_results, do,
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
        classical = BayesGLM_results$cortexL$result_classical[[ss]]$estimates,
        Bayesian = BayesGLM_results$cortexL$field_estimates[[ss]]
      )
      # Update mwall b/c `mask2` in `BayesGLM` can change the medial wall.
      mwallL <- BayesGLM_results$cortexL$spatial$mask
      colnames(datL) <- NULL
    }
    if (do$right) {
      datR <- switch(method,
        classical = BayesGLM_results$cortexR$result_classical[[ss]]$estimates,
        Bayesian = BayesGLM_results$cortexR$field_estimates[[ss]]
      )
      mwallR <- BayesGLM_results$cortexR$spatial$mask
      colnames(datR) <- NULL
    }
    if (do$sub) {
      datSub <- switch(method,
        classical = BayesGLM_results$subcort$result_classical[[ss]]$estimates,
        Bayesian = BayesGLM_results$subcortical$field_estimates[[ss]]
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

# #' Format BayesGLM results into \code{"xifti"} objects for the multi-model case
# #'
# #' Format beta estimates of the BayesGLM results for the multi-model case into
# #'  list of \code{"xifti"} objects for each session.
# #'
# #' @param BayesGLM_results The list of BayesGLM results
# #' @param session_names The session names
# #' @param method \code{"classical"} or \code{"Bayesian"}
# #' @return The list of \code{"xifti"} objects
# #' @keywords internal
# #'
# BayesGLM_cifti_format_results_multi <- function(
#   BayesGLM_results, session_names, method=c("classical", "Bayesian")){

#     datL <- datR <- datSub <- NULL #index of best-fitting model
#     betaL <- betaR <- betaSub <- NULL #beta estimates for best-fitting model
#     sigma2L <- sigma2R <- sigma2Sub <- NULL #residual var of models

#     for (ss in seq(nS)) {
#       # INDEX OF BEST MODEL
#       if (do$left) {
#         datL <- BayesGLM_results$cortex_left$result_multiple[[ss]]$bestmodel #index of best model
#         betaL <- BayesGLM_results$cortex_left$result_multiple[[ss]]$beta_estimates #V x K x P (P = number of models tested)
#         sigma2L <- BayesGLM_results$cortex_left$result_multiple[[ss]]$sigma2 #V x P (P = number of models tested)
#         #only save beta estimates for the best fitting model
#         betaL <- apply(betaL, 1, as.matrix, simplify=FALSE) #form into a list of length V, each a K x P matrix
#         betaL <- t(mapply(function(matrix, index) matrix[, index, drop = FALSE], betaL, datL, SIMPLIFY = TRUE)) #beta estimates (VxK) for the best model
#         # Update mwall b/c `mask2` in `BayesGLM` can change the medial wall.
#         mwallL <- !is.na(datL)
#         datL <- datL[mwallL]
#         betaL <- betaL[mwallL,]
#         sigma2L <- sigma2L[mwallL,]
#       }
#       if (do$right) {
#         datR <- BayesGLM_results$cortex_right$result_multiple[[ss]]$bestmodel
#         betaR <- BayesGLM_results$cortex_right$result_multiple[[ss]]$beta_estimates #V x K x P (P = number of models tested)
#         sigma2R <- BayesGLM_results$cortex_right$result_multiple[[ss]]$sigma2 #V x P (P = number of models tested)
#         #only save beta estimates for the best fitting model
#         betaR <- apply(betaR, 1, as.matrix, simplify=FALSE) #form into a list of length V, each a K x P matrix
#         betaR <- t(mapply(function(matrix, index) matrix[, index, drop = FALSE], betaR, datR, SIMPLIFY = TRUE)) #beta estimates for the best model
#         mwallR <- !is.na(datR)
#         datR <- datR[mwallR]
#         betaR <- betaR[mwallR,]
#         sigma2R <- sigma2R[mwallR,]
#       }
#       if (do$sub) {
#         #[TO DO]: do this for subcortex, as for L and R above
#         #datSub <- BayesGLM_results$subcortical$result_multiple[[ss]]$bestmodel
#         #colnames(datSub) <- NULL
#       }

#       bestmodel_xii[[ss]] <- as.xifti(
#         cortexL = datL,
#         cortexL_mwall = if (do$left) { mwallL } else { NULL },
#         cortexR = datR,
#         cortexR_mwall = if (do$right) { mwallR } else { NULL },
#         subcortVol = datSub,
#         subcortLabs = submeta$labels,
#         subcortMask = submeta$mask
#       )
#       bestmodel_xii[[ss]]$meta$subcort$trans_mat <- submeta$trans_mat
#       bestmodel_xii[[ss]]$meta$subcort$trans_units <- submeta$trans_units
#       bestmodel_xii[[ss]]$meta$cifti$names <- field_names

#       results_xii$classical[[ss]] <- as.xifti(
#         cortexL = betaL,
#         cortexL_mwall = if (do$left) { mwallL } else { NULL },
#         cortexR = betaR,
#         cortexR_mwall = if (do$right) { mwallR } else { NULL },
#         subcortVol = betaSub,
#         subcortLabs = submeta$labels,
#         subcortMask = submeta$mask
#       )
#       results_xii$classical[[ss]]$meta$subcort$trans_mat <- submeta$trans_mat
#       results_xii$classical[[ss]]$meta$subcort$trans_units <- submeta$trans_units
#       results_xii$classical[[ss]]$meta$cifti$names <- field_names

#       sigma2_xii[[ss]] <- as.xifti(
#         cortexL = sigma2L,
#         cortexL_mwall = if (do$left) { mwallL } else { NULL },
#         cortexR = sigma2R,
#         cortexR_mwall = if (do$right) { mwallR } else { NULL },
#         subcortVol = sigma2Sub,
#         subcortLabs = submeta$labels,
#         subcortMask = submeta$mask
#       )
#       sigma2_xii[[ss]]$meta$subcort$trans_mat <- submeta$trans_mat
#       sigma2_xii[[ss]]$meta$subcort$trans_units <- submeta$trans_units
#       sigma2_xii[[ss]]$meta$cifti$names <- field_names
#     }

#   list(
#     bestmodel_xii = bestmodel_xii,
#     result_xii = result_xii,
#     sigma2_xii = sigma2_xii
#   )
# }
