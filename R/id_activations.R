#' Identify task activations
#'
#' Identify areas of activation for each task from the result of \code{BayesGLM}
#'  or \code{BayesGLM_cifti}.
#'
#' @param model_obj Result of \code{BayesGLM} or \code{BayesGLM_cifti} model
#'  call, of class \code{"BayesGLM"} or \code{"BayesGLM_cifti"}.
# @param method The method to be used for identifying activations, either 'posterior' (Default) or '2means'
#' @param tasks The task(s) to identify activations for. Give either the name(s)
#'  as a character vector, or the numerical indices. If \code{NULL} (default),
#'  analyze all tasks.
#' @param sessions The session(s) to identify activations for. Give either the
#'  name(s) as a character vector, or the numerical indices. If \code{NULL}
#' (default), analyze the first session.
#'
#'  Currently, if multiple sessions are provided, activations are identified
#'  separately for each session. (Information is not combined between the
#'  different sessions.)
#' @param method \code{"Bayesian"} (default) or \code{"classical"}. If
#'  \code{model_obj} does not have Bayesian results because \code{Bayes} was set
#'  to \code{FALSE}, only the \code{"classical"} method can be used.
#' @param alpha Significance level. Default: \code{0.05}.
#' @param threshold Activation threshold, for example \code{1} for 1\% signal
#'  change if \code{scale_BOLD=="mean"} during model estimation. Setting a
#'  \code{threshold} is required for the Bayesian method; \code{NULL} (default)
#'  will use a \code{threshold} of zero for the classical method.
# @param excur_method For method = 'Bayesian' only: Either \code{EB} (empirical Bayes) or \code{QC} (Quantile Correction), depending on the method that should be used to find the
#   excursions set. Note that if any contrasts (including averages across sessions) are used in the modeling, the method chosen must be \code{EB}.
#   The difference in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
# @param area.limit For method = 'Bayesian' only: Below this value, clusters of activations will be considered spurious.  If NULL (default), no limit.
#' @param correction For the classical method only: Type of multiple comparisons
#'  correction: \code{"FWER"} (Bonferroni correction, the default), \code{"FDR"}
#'  (Benjamini Hochberg), or \code{"none"}.
# @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
#  depending on the method that should be used to find the excursions set. Note that to ID
#  activations for averages across sessions, the method chosen must be \code{EB}. The difference
#  in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#' @inheritParams verbose_Param
# @param type For method='2means' only: The type of 2-means clustering to perform ('point' or 'sequential')
# @param n_sample The number of samples to generate if the sequential 2-means type is chosen. By default, this takes a value of 1000.
#'
# @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @importFrom ciftiTools convert_xifti
#' @importFrom fMRItools is_posNum is_1
#'
#' @return An \code{"act_BayesGLM"} or \code{"act_BayesGLM_cifti"} object, a 
#'  list which indicates the activated locations along with related information.
#' @export
#'
id_activations <- function(
  model_obj,
  tasks=NULL,
  sessions=NULL,
  method=c("Bayesian", "classical"),
  alpha=0.05,
  threshold=NULL,
  #area.limit=NULL,
  correction = c("FWER", "FDR", "none"),
  #excur_method = c("EB", "QC"),
  verbose = 1){

  # If 'BayesGLM_cifti', we will loop over the brain structures.
  is_cifti <- inherits(model_obj, "BayesGLM_cifti")
  if (is_cifti) {
    cifti_obj <- model_obj
    model_obj <- cifti_obj$BayesGLM_results
  } else {
    if (!inherits(model_obj, "BayesGLM")) {
      stop("`model_obj` is not a `'BayesGLM'` or 'BayesGLM_cifti' object.")
    }
    model_obj <- list(obj=model_obj)
  }
  idx1 <- min(which(!vapply(model_obj, is.null, FALSE)))

  # Argument checks.
  stopifnot(is.null(tasks) || is.character(tasks) || is.numeric(tasks))
  if (is.character(tasks)) { stopifnot(all(tasks %in% model_obj[[idx1]]$task_names)) }
  if (is.numeric(tasks)) {
    stopifnot(all.equal(tasks, round(tasks))==TRUE)
    stopifnot(min(tasks) >= 1)
    stopifnot(max(tasks) <= length(model_obj[[idx1]]$task_names))
    stopifnot(length(tasks) == length(unique(tasks)))
    tasks <- model_obj[[idx1]]$task_names[tasks]
  }
  if (is.null(tasks)) { tasks <- model_obj[[idx1]]$task_names }
  stopifnot(is.null(sessions) || is.character(sessions) || is.numeric(sessions))
  if (is.character(sessions)) { stopifnot(all(sessions %in% model_obj[[idx1]]$session_names)) }
  if (is.numeric(sessions)) {
    stopifnot(all.equal(sessions, round(sessions))==TRUE)
    stopifnot(min(sessions) >= 1)
    stopifnot(max(sessions) <= length(model_obj[[idx1]]$session_names))
    stopifnot(length(sessions) == length(unique(sessions)))
    sessions <- model_obj[[idx1]]$session_names[sessions]
  }
  if (is.null(sessions)) { sessions <- model_obj[[idx1]]$session_names[1] }
  method <- match.arg(method, c("Bayesian", "classical"))
  stopifnot(is_posNum(alpha) && alpha < 1)
  stopifnot(is.null(threshold) || is_posNum(threshold, zero_ok=TRUE))
  correction <- match.arg(correction, c("FWER", "FDR", "none"))
  stopifnot(is_posNum(verbose, zero_ok=TRUE))

  # Check that Bayesian results are available, if requested.
  if (method=="Bayesian" && is.null(model_obj[[idx1]]$INLA_model_obj)) {
    warning("`method=='Bayesian'` but only classical model results are available. Setting `method` to `'classical'`.")
    method <- "classical"
  }

  # Check threshold was set, if computing Bayesian activations.
  if (is.null(threshold)) {
    if (method=='Bayesian') stop("Must specify an activation `threshold` when `method=='Bayesian'`.")
    if (method=='classical') threshold <- 0
  }

  # Initialize list of activations.
  nModels <- length(model_obj)
  activations <- vector('list', length=nModels)
  names(activations) <- names(model_obj)

  actFUN <- switch(method,
    Bayesian = id_activations.posterior,
    classical = id_activations.classical
  )
  actArgs <- list(tasks=tasks, alpha=alpha, threshold=threshold)
  if (method == "classical") {
    actArgs <- c(actArgs, list(correction=correction))
  } else {
    correction <- "not applicable"
  }

  # Loop over model objects (brain structures, if `is_cifti`).
  for (mm in seq(nModels)) {
    if (is.null(model_obj[[mm]])) { next }
    if (method=="Bayesian" && identical(attr(model_obj[[mm]]$INLA_model_obj, "format"), "minimal")) {
      stop("Bayesian activations are not available because `return_INLA` was set to `'minimal'` in the `BayesGLM` call. Request the classical activations, or re-run `BayesGLM`.")
    }
    name_obj_mm <- names(model_obj)[mm]
    if (method=="Bayesian" && name_obj_mm!="obj") {
      if (verbose>0) cat(paste0("Identifying Bayesian GLM activations in ",name_obj_mm,'\n'))
    }
    # Loop over sessions
    activations[[mm]] <- vector("list", length(sessions))
    names(activations[[mm]]) <- sessions
    for (session in sessions) {
      actArgs_ms <- c(actArgs, list(model_obj=model_obj[[mm]], session=session))
      activations[[mm]][[session]] <- do.call(actFUN, actArgs_ms)
    }
  }

  result <- list(
    activations = activations,
    method=method,
    alpha=alpha,
    threshold=threshold,
    correction=correction,
    task_names = model_obj[[idx1]]$task_names,
    session_names = model_obj[[idx1]]$session_names
    #excur_method = c("EB", "QC")
  )

  # If BayesGLM, return.
  if (!is_cifti) {
    result$activations <- result$activations[[1]]
    class(result) <- "act_BayesGLM"
    return(result)
  }

  # If BayesGLM_cifti, create 'xifti' with activations.
  act_xii <- vector("list", length(sessions))
  names(act_xii) <- sessions
  for (session in sessions) {
    the_xii <- cifti_obj$task_estimates_xii$classical[[session]]
    act_xii_ss <- 0*select_xifti(the_xii, match(tasks, the_xii$meta$cifti$names))
    for (bs in names(the_xii$data)) {
      if (bs=="subcort") { next }
      if (!is.null(the_xii$data[[bs]])) {
        dat <- 1*activations[[bs]][[session]]$active
        colnames(dat) <- NULL
        if (method=="classical") { dat <- dat[!is.na(dat[,1]),,drop=FALSE] }
        act_xii_ss$data[[bs]] <- dat
      }
    }
    # if (!is.null(model_obj$subcortical)) {
    #   if(method == "EM") {
    #     for(m in 1:length(activations$subcortical)){
    #       datS <- 1*activations$subcortical[[m]]$active
    #       act_xii_ss$data$subcort[!is.na(model_obj$betas_EM[[1]]$data$subcort[,1]),] <-
    #         datS
    #     }
    #   }
    #   if(method == "classical") {
    #     datS <- 1*activations$subcortical$active
    #     act_xii_ss$data$subcort[!is.na(model_obj$betas_classical[[1]]$data$subcort[,1]),] <-
    #       datS
    #   }
    #   act_xii_ss$data$subcort <- matrix(datS, ncol=length(tasks))
    # }
    act_xii_ss <- convert_xifti(
      act_xii_ss, "dlabel", 
      values=c(Inactive=0, Active=1), colors='red'
    )
    act_xii[[session]] <- act_xii_ss
  }

  result <- c(list(activations_xii=act_xii), result)
  class(result) <- "act_BayesGLM_cifti"
  result
}

#' Identify activations using joint posterior probabilities
#'
#' Identifies areas of activation given an activation threshold and significance
#'  level using joint posterior probabilities
#'
#' For a given latent field, identifies locations that exceed a certain activation
#'  threshold (e.g. 1 percent signal change) at a given significance level, based on the joint
#'  posterior distribution of the latent field.
#'
#' @param model_obj Result of \code{BayesGLM}, of class \code{"BayesGLM"}.
#' @param tasks,session,alpha,threshold See \code{\link{id_activations}}.
# @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
# depending on the method that should be used to find the excursions set. Note that to ID
# activations for averages across sessions, the method chosen must be \code{EB}. The difference
# in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#'
#' @return A list with two elements: \code{active}, which gives a matrix of zeros
#' and ones of the same dimension as \code{model_obj$task_estimates${session}},
#' and \code{excur_result}, an object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for
#'  more information).
#'
#' @importFrom excursions excursions.inla
#'
#' @keywords internal
id_activations.posterior <- function(
  model_obj,
  tasks, session,
  alpha=0.05, threshold){
  #area.limit=NULL,
  #excur_method = c("EB","QC")

  stopifnot(inherits(model_obj, "BayesGLM"))
  #excur_method <- match.arg(excur_method, c("EB","QC"))

  sess_ind <- which(model_obj$session_names == session)
  mesh <- model_obj$mesh
  n_vox <- mesh$n
  #indices of beta vector corresponding to specified session
  inds <- (1:n_vox) + (sess_ind-1)*n_vox

	#loop over latent fields
	excur <- vector('list', length=length(tasks))
	act <- matrix(NA, nrow=n_vox, ncol=length(tasks))
	colnames(act) <- tasks
	for(f in tasks){

		#if(is.null(area.limit)){
			res.exc <- excursions.inla(
        model_obj$INLA_model_obj,
        name=f, ind=inds, u=threshold, type='>', alpha=alpha, method="EB"
      )
		#} else {
		#	res.exc <- excursions.inla.no.spurious(model_obj$INLA_model_obj, mesh=mesh, name=f, ind=inds, u=threshold, type='>', method=excur_method, alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
		#}
	  which_f <- which(tasks==f)
		act[,which_f] <- res.exc$E[inds] == 1
    excur[[which_f]] <- res.exc
	}
	result <- list(active=act, excursions_result=excur)

  #compute size of activations
  areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
  areas_act <- apply(act, 2, function(x) sum(areas_all[x==1]))

  result$areas_all <- areas_all
  result$areas_act <- areas_act

	result
}

#' Identification of areas of activation in a General Linear Model using classical methods
#'
#' @param model_obj A \code{BayesGLM} object
#' @param tasks,session,alpha,threshold See \code{\link{id_activations}}.
#' @param correction (character) Either 'FWER' or 'FDR'. 'FWER' corresponds to the
#'   family-wise error rate with Bonferroni correction, and 'FDR' refers to the
#'   false discovery rate using Benjamini-Hochberg.
#' @param mesh (Optional) An \code{"inla.mesh"} object (see \code{\link{make_mesh}} for
#'  surface data). Only necessary for computing surface areas of identified activations.
#'
#' @return A matrix corresponding to the
#'   0-1 activation status for the model coefficients.
#'
#' @importFrom stats sd pt p.adjust
#' @importFrom matrixStats colVars
#'
#' @keywords internal
id_activations.classical <- function(model_obj,
                                     tasks,
                                     session,
                                     alpha = 0.05,
                                     threshold = 0,
                                     correction = c("FWER", "FDR", "none"),
                                     mesh = NULL) {

  # Argument checks ------------------------------------------------------------
  if (!inherits(model_obj, "BayesGLM")) {
    stop(paste0(
      "The model object is of class ",
      paste0(class(model_obj), collapse=", "),
      " but should be of class 'BayesGLM'."
    ))
  }
  # tasks, session, alpha, threshold checked in `id_activations`
  correction <- match.arg(correction, c("FWER","FDR","none"))

  beta_est <- model_obj$result_classical[[session]]$estimates
  se_beta <- model_obj$result_classical[[session]]$SE_estimates
  DOF <- model_obj$result_classical[[session]]$DOF

  nvox <- nrow(beta_est)
  #if(any(!(tasks %in% 1:K))) stop(paste0('tasks must be between 1 and the number of tasks, ',K))
  beta_est <- matrix(beta_est[,tasks], nrow=nvox) #need matrix() in case beta_est[,tasks] is a vector
  se_beta <- matrix(se_beta[,tasks], nrow=nvox) #need matrix() in case beta_est[,tasks] is a vector
  K <- length(tasks)

  #Compute t-statistics and p-values
  t_star <- (beta_est - threshold) / se_beta
  if(!is.matrix(t_star)) t_star <- matrix(t_star, nrow=nvox)
  #perform multiple comparisons correction
  p_values <- p_values_adj <- active <- matrix(NA, nvox, K)

  for (kk in 1:K) {
    p_values_k <- sapply(t_star[,kk], pt, df = DOF, lower.tail = F)
    p_vals_adj_k <- switch(correction,
      FWER = p.adjust(p_values_k, method='bonferroni'),
      FDR = p.adjust(p_values_k, method='BH'),
      none = p_values_k
    )
    p_values[,kk] <- p_values_k
    p_values_adj[,kk] <- p_vals_adj_k
    active[,kk] <- (p_vals_adj_k < alpha)
  }

  #compute size of activations
  if(!is.null(mesh)){
    mask <- model_obj$result_classical[[session]]$mask
    if(sum(mask) != nrow(mesh$loc)) stop('Supplied mesh is not consistent with mask in model_obj.')
    areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
    areas_act <- apply(active[mask==1,,drop=FALSE], 2, function(x) sum(areas_all[x==1]))
  } else {
    areas_all <- areas_act <- NULL
  }

  na_pvalues <- which(is.na(p_values[,1]))

  result <- list(
    p_values = p_values[-na_pvalues,, drop = F],
    p_values_adj = p_values_adj[-na_pvalues,, drop = F],
    active = active[-na_pvalues,, drop = F],
    correction = correction,
    alpha = alpha,
    threshold = threshold,
    areas_all = areas_all,
    areas_act = areas_act
  )
  result
}

# #' Identify activations using joint posterior probabilities with EM results
# #'
# #' Identifies areas of activation given an activation threshold and significance
# #'  level using joint posterior probabilities
# #'
# #' For a given latent field, identifies locations that exceed a certain activation
# #'  threshold (e.g. 1 percent signal change) at a given significance level, based on the joint
# #'  posterior distribution of the latent field.
# #'
# #' @param model_obj An object of class \code{"BayesGLM"}, a result of a call
# #'  to \code{BayesGLMEM}.
# #' @param tasks Name of latent field or vector of names on which to identify activations.  By default, analyze all tasks.
# #' @param sessions (character) The name of the session that should be examined.
# #' If \code{NULL} (default), the first session is used.
# #' @param alpha Significance level (e.g. 0.05)
# #' @param threshold Activation threshold (e.g. 1 for 1 percent signal change if scale=TRUE in model estimation)
# #' @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
# #'
# #'
# #' @return A list with two elements: \code{active}, which gives a matrix of zeros
# #' and ones of the same dimension as \code{model_obj$task_estimates${sessions}},
# #' and \code{excur_result}, an object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for
# #'  more information).
# #'
# #' @importFrom excursions excursions
# #' @importFrom stats na.omit
# #'
# #' @keywords internal
# id_activations.em <-
#   function(model_obj,
#            tasks = NULL,
#            sessions = NULL,
#            alpha = 0.05,
#            threshold,
#            area.limit = NULL) {
#     if (!inherits(model_obj, "BayesGLM"))
#       stop(paste0(
#         "The model object is of class ",
#         class(model_obj),
#         " but should be of class 'BayesGLM'."
#       ))
#     #check sessions argument

#     #if only one session, analyze that one
#     all_sessions <- model_obj$sessions
#     n_sess <- length(all_sessions)
#     if (n_sess == 1) {
#       sessions <- all_sessions
#     }

#     #check sessions not NULL, check that a valid session name
#     if(!is.null(sessions)){
#       if(!(sessions %in% all_sessions)) stop(paste0('sessions does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
#       sess_ind <- which(all_sessions == sessions)
#     }

#     #if averages not available and sessions=NULL, pick first session and return a warning
#     # has_avg <- is.matrix(model_obj$avg_task_estimates)
#     has_avg <- !is.null(model_obj$avg_task_estimates)
#     if(is.null(sessions) & !has_avg){
#       sessions <- all_sessions[1]
#       warning(paste0("Your model object does not have averaged beta estimates. Using the first session instead. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
#     }

#     if(is.null(sessions) & has_avg){
#       sessions <- 'avg'
#     }

#     # Make an indicator for subcortical data
#     is_subcort <- "BayesGLMEM_vol3D" %in% as.character(model_obj$call)

#     #if sessions is still NULL, use average and check excur_method argument
#     # if(is.null(sessions)){
#     #   if(has_avg & excur_method != "EB") {
#     #     excur_method <- 'EB'
#     #     warning("To id activations for averaged beta estimates, only the excur_method='EB' is supported. Setting excur_method to 'EB'.")
#     #   }
#     # }

#     #check tasks argument
#     if(is.null(tasks)) tasks <- model_obj$tasks
#     if(!any(tasks %in% model_obj$tasks)) stop(paste0("Please specify only field names that corresponds to one of the latent fields: ",paste(model_obj$tasks, collapse=', ')))

#     #check alpha argument
#     if(alpha > 1 | alpha < 0) stop('alpha value must be between 0 and 1, and it is not')

#     mesh <- model_obj$mesh
#     if(is_subcort) {
#       n_subcort_models <- length(model_obj$spde_obj)
#       result <- vector("list", length = n_subcort_models)
#       for(m in 1:n_subcort_models){
#         # mesh <- make_mesh(model_obj$spde_obj[[m]]$vertices[[1]],
#         #                   model_obj$spde_obj[[m]]$faces[[1]])
#         # n_vox <- mesh$n
#         # excur <- vector('list', length=length(tasks))
#         # act <- matrix(NA, nrow=n_vox, ncol=length(tasks))
#         # colnames(act) <- tasks

#         if(sessions == "avg") {
#           beta_est <- c(model_obj$EM_result_all[[m]]$posterior_mu)
#         }
#         if(sessions != "avg") {
#           beta_mesh_est <- c(model_obj$task_estimates[[which(all_sessions == sessions)]])
#           beta_est <- beta_mesh_est[!is.na(beta_mesh_est)]
#         }
#         Phi_k <- model_obj$EM_result_all[[m]]$mesh$Amat
#         # Phi <-
#         #   Matrix::bdiag(rep(
#         #     list(Phi_k),
#         #     length(tasks)
#         #   ))
#         # w <- beta_est %*% Phi
#         Sig_inv <- model_obj$EM_result_all[[m]]$posterior_Sig_inv
#         # Q_est <- Matrix::tcrossprod(Phi %*% Sig_inv, Phi)

#         ### use the ind argument for excursions to get the marginal posterior
#         ### excursions sets (much faster)

#         n_vox <- length(beta_est) / length(tasks)
#         act <- matrix(NA, nrow=n_vox, ncol=length(tasks))

#         V <- Matrix::diag(INLA::inla.qinv(Sig_inv))

#         for(f in tasks) {
#           which_f <- which(tasks==f)
#           f_inds <- (1:n_vox) + (which_f - 1)*n_vox
#           res.exc <-
#             excursions(
#               alpha = alpha,
#               u = threshold,
#               mu = beta_est,
#               Q = Sig_inv,
#               vars = V,
#               ind = f_inds,
#               type = ">", method = "EB"
#             )

#           act[,which_f] <- res.exc$E[f_inds]
#         }
#         act <- as.matrix(Phi_k %*% act)
#         result[[m]] <- list(active = act, excursions_result=res.exc)
#       }
#     }

#     if(!is_subcort) {
#       n_vox <- mesh$n

#       #for a specific session
#       if(!is.null(sessions)){

#         # inds <- (1:n_vox) + (sess_ind-1)*n_vox #indices of beta vector corresponding to session v

#         #loop over latent fields
#         excur <- vector('list', length=length(tasks))
#         act <- matrix(NA, nrow=n_vox, ncol=length(tasks))
#         colnames(act) <- tasks
#         res.exc <-
#           excursions(
#             alpha = alpha,
#             u = threshold,
#             mu = stats::na.omit(
#               c(model_obj$task_estimates[[sessions]])
#             ),
#             Q = model_obj$posterior_Sig_inv,
#             type = ">", method = "EB"
#           )
#         for(f in tasks) {
#           which_f <- which(tasks==f)
#           f_inds <- (1:n_vox) + (which_f - 1)*n_vox
#           act[,which_f] <- res.exc$E[f_inds]
#         }
#       }
#       result <- list(active=act, excursions_result=res.exc)
#     }

#     return(result)
#   }
