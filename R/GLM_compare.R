#' Classical GLM for multiple models
#'
#' Classical GLM for multiple models
#'
# @param data,spatial,spatial_type,session_names,field_names,design_type See \code{fit_bayesglm}.
# @param valid_cols,nT,nD,var_resid,sqrtInv_all See \code{fit_bayesglm}.
# @return A list of results
#' @keywords internal
GLM_compare <- function(...){NULL}
#   data, spatial, spatial_type,
#   session_names, field_names, design_type,
#   valid_cols, nT, nD,
#   var_resid, sqrtInv_all
# ){

#   nS <- length(session_names)
#   nK <- length(field_names)
#   nV <- get_nV(spatial)

#   out <- setNames(vector('list', length=nS), session_names)
#   for (ss in seq(nS)) {

#     beta_hat_ss <- array(
#       NA,
#       dim=c(nV$T, nK, nD),
#       dimnames = list(loc = 1:nV$T, field = field_names, model = seq(nD))
#     )
#     # Keep track of residual SD (proxy for R^2 or AIC)
#     sigma2_ss <- matrix(NA, nrow=nV$T, ncol=nD)

#     for (dd in seq(nD)) {
#       cat(paste0('\tFitting model ',dd,'\n'))

#       x <- GLM_classical(
#         data, spatial, spatial_type,
#         session_names, field_names, design_type,
#         valid_cols[ss,], nT[ss],
#         do_pw=FALSE, compute_SE=FALSE
#       )
#       beta_hat_ss[,,dd] <- x$estimates
#       sigma2_ss[spatial$maskD,dd] <- sqrt(
#         colSums(x$resids^2)/(nT[ss] - sum(spatial$maskD)) # [TO DO] check!!!
#       )
#     }

#     #determine best model (minimum residual error)
#     bestmodel_ss <- apply(sigma2_ss, 1, function(x){
#       wm <- which.min(x)
#       varx <- var(x, na.rm=TRUE)
#       if(is.na(varx)) varx <- 0
#       if(varx==0) wm <- NA
#       wm
#     })

#     out[[ss]] <- list(
#       beta_estimates = beta_hat_ss,
#       bestmodel = bestmodel_ss,
#       sigma2 = sigma2_ss
#     )
#   }

#   result <- list(
#     field_estimates = lapply(out, '[[', "beta_estimates"),
#     result_multiple = out#,
#     # mesh = mesh,
#     # spde = spde,
#     # mesh_orig = mesh_orig,
#     # mask = mask,
#     # mask_orig = mask_orig, #[TO DO] return the params passed into the function instead?
#     # mask_qc = mask_qc,
#     # design = lapply(data, function(ss){ss$design}),
#     # field_names = field_names,
#     # session_names = session_names,
#     # n_sess_orig = nS_orig,
#     # call = match.call()
#   )
#   class(result) <- "CompareGLM"

#   result
#}
