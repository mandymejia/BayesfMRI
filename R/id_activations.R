#' Identify Task Activations
#'
#' Identify areas of activation for each task
#'
#' @param model_obj Result of \code{BayesGLM_cifti} model call, of class
#'  \code{"BayesGLM_cifti"}
# @param method The method to be used for identifying activations, either 'posterior' (Default) or '2means'
#' @param field_names The name(s) of the latent field(s) on which to 
#'  identify activations. If \code{NULL} (default), analyze all tasks.
#' @param session_name The name of the session to examine. 
#' If \code{NULL} (default), the average across all sessions is used.
#' @param alpha Significance level. Default: \code{0.05}.
#' @param method \code{"Bayesian"} (default) or \code{"classical"}.
#' @param threshold Activation threshold, for example \code{1} for 1\% signal
#'  change if \code{scale_BOLD=="mean"} during model estimation. Setting a 
#'  threshold is required for the Bayesian method; \code{NULL} (default) will 
#'  use a threshold of zero for the classical method.
# @param excur_method For method = 'Bayesian' only: Either \code{EB} (empirical Bayes) or \code{QC} (Quantile Correction), depending on the method that should be used to find the
#   excursions set. Note that if any contrasts (including averages across sessions) are used in the modeling, the method chosen must be \code{EB}.
#   The difference in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
# @param area.limit For method = 'Bayesian' only: Below this value, clusters of activations will be considered spurious.  If NULL (default), no limit.
#' @param correction For the classical method only: Type of multiple comparisons
#'  correction: \code{"FWER"} (Bonferroni correction, the default), \code{"FDR"}
#'  (Benjamini Hochberg), or \code{"none"}.
# @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
#' depending on the method that should be used to find the excursions set. Note that to ID
#' activations for averages across sessions, the method chosen must be \code{EB}. The difference
#' in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#' @param verbose Print progress updates? Default: \code{TRUE}.
# @param type For method='2means' only: The type of 2-means clustering to perform ('point' or 'sequential')
# @param n_sample The number of samples to generate if the sequential 2-means type is chosen. By default, this takes a value of 1000.
#'
# @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @importFrom ciftiTools convert_xifti
#' @export
#'
#'
id_activations_cifti <- function(
  model_obj,
  field_names=NULL,
  session_name=NULL,
  alpha=0.05,
  method=c("Bayesian", "classical"),
  threshold=NULL,
  #area.limit=NULL,
  correction = c("FWER", "FDR", "none"),
  #excur_method = c("EB", "QC"),
  verbose = TRUE){

  # Simple argument checks


  if (!inherits(model_obj, "BayesGLM_cifti")) {
    stop("`model_obj` must be of class 'BayesGLM_cifti' (output of `BayesGLM_cifti` function).")
  }

  ## HERE -- figure out how to trigger NULL if Bayes=FALSE  (also below when adding the mesh)
  ## If no Bayesian results but method=Bayesian, should return error (or switch to classical with a warning)

  if(is.null(model_obj$task_estimates[[1]])) method <- 'classical'
  method <- match.arg(method, c('Bayesian','classical'))
  if(!(method %in% c('Bayesian','classical'))) stop("The method argument should be 'Bayesian', 'EM', or 'classical'.")

  # [TO DO]: check that the requested method(s) actually exist in model_obj?
  # what's the best way to check this? all entries in e.g. `$betas_Bayesian` are `NULL`?

  if(is.null(threshold)){
    if(method=='classical') threshold <- 0
    if(method=='Bayesian') stop("Must specify an activation threshold when method='Bayesian'.")
  }

  GLM_list <- model_obj$BayesGLM_results
  num_models <- length(GLM_list)
  activations <- vector('list', length=num_models)
  names(activations) <- names(GLM_list)
  do_left <- !is.null(GLM_list$cortexL)
  do_right <- !is.null(GLM_list$cortexR)
  do_sub <- !is.null(GLM_list$subcortical)

  if (is.null(field_names))
    field_names <- model_obj$task_names
  if (any(!(field_names %in% model_obj$task_names)))
    stop(paste0(
      'All elements of field_names must appear in model_obj$task_names: ',
      paste(model_obj$task_names, collapse = ',')
    ))
  field_inds <- which(model_obj$task_names %in% field_names)

  if(method == "Bayesian" && !("INLA_model_obj" %in% do.call(c, lapply(GLM_list, names)))) {
    method <- "EM"
  }

  if(method=='Bayesian'){
    #loop over hemispheres/structures
    for(mm in 1:num_models){
      if(is.null(GLM_list[[mm]])) next
      model_m <- GLM_list[[mm]]
      if(verbose) cat(paste0('Identifying Bayesian GLM activations in ',names(GLM_list)[mm],'\n'))
      act_m <- id_activations.posterior(model_obj=model_m,
                                        field_names=field_names,
                                        session_name=session_name,
                                        threshold=threshold,
                                        alpha=alpha)
                                        #area.limit=area.limit,
                                        #excur_method=excur_method)

      activations[[mm]] <- act_m
    }
  }

  if(method=='classical'){
    #if(verbose) cat('Identifying classical GLM activations')
    for(mm in 1:num_models){
      if(is.null(GLM_list[[mm]])) next
      model_m <- GLM_list[[mm]]
      mesh_m <- model_obj$GLMs_Bayesian[[mm]]$mesh #will be NULL if GLMs_Bayesian is NULL (no Bayesian model was fit)
      act_m <- id_activations.classical(model_obj=model_m,
                                        field_inds=field_inds,
                                        session_name=session_name,
                                        threshold=threshold,
                                        alpha=alpha,
                                        correction=correction,
                                        mesh=mesh_m)

      activations[[mm]] <- act_m
    }
  }

  # if(method=='EM'){
  #   for(mm in 1:num_models){
  #     if(is.null(GLM_list[[mm]])) next
  #     model_m <- GLM_list[[mm]]
  #     if(verbose) cat(paste0('Identifying EM GLM activations in ',names(GLM_list)[mm],'\n'))
  #     act_m <- id_activations.em(model_obj=model_m,
  #                                field_names=field_names,
  #                                session_name=session_name,
  #                                threshold=threshold,
  #                                alpha=alpha,
  #                                area.limit=NULL)

  #     activations[[mm]] <- act_m
  #   }
  # }

  #map results to xifti objects
  activations_xifti <- 0*model_obj[[paste0("betas_", method[1])]][[1]]
  if(length(field_names) != length(model_obj$task_names)) activations_xifti$meta$cifti$names <- field_names
  if(do_left) {
    datL <- 1*activations$cortexL$active
    if(method=='classical') datL <- datL[!is.na(datL[,1]),] #remove medial wall locations
    #datL[datL==0] <- NA
    activations_xifti$data$cortex_left <- matrix(datL, ncol=length(field_names))
  }
  if(do_right) {
    datR <- 1*activations$cortexR$active
    if(method=='classical') datR <- datR[!is.na(datR[,1]),] #remove medial wall locations
    #datR[datR==0] <- NA
    activations_xifti$data$cortex_right <- matrix(datR, ncol=length(field_names))
  }
  if(do_sub) {
    if(method == "EM") {
      for(m in 1:length(activations$subcortical)){
        datS <- 1*activations$subcortical[[m]]$active
        activations_xifti$data$subcort[!is.na(model_obj$betas_EM[[1]]$data$subcort[,1]),] <-
          datS
      }

    }
    if(method == "classical") {
      datS <- 1*activations$subcortical$active
      activations_xifti$data$subcort[!is.na(model_obj$betas_classical[[1]]$data$subcort[,1]),] <-
        datS
    }
    activations_xifti$data$subcort <- matrix(datS, ncol=length(field_names))
  }

  activations_xifti <- convert_xifti(activations_xifti, "dlabel", colors='red')

  result <- list(activations = activations,
                 activations_xifti = activations_xifti)

  return(result)
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
#' @param model_obj An object of class ‘"BayesGLM"’, a result of a call to BayesGLM
#' @param field_names Name of latent field or vector of names on which to identify activations.  By default, analyze all tasks.
#' @param session_name (character) The name of the session that should be examined.
#' If \code{NULL} (default), the average across all sessions is used.
#' @param alpha Significance level (e.g. 0.05)
#' @param threshold Activation threshold (e.g. 1 for 1 percent signal change if scale=TRUE in model estimation)
# @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
# @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
# depending on the method that should be used to find the excursions set. Note that to ID
# activations for averages across sessions, the method chosen must be \code{EB}. The difference
# in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#'
#'
#' @return A list with two elements: \code{active}, which gives a matrix of zeros
#' and ones of the same dimension as \code{model_obj$task_estimates${session_name}},
#' and \code{excur_result}, an object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for
#'  more information).
#'
#' @importFrom excursions excursions.inla
#'
#' @export
id_activations.posterior <- function(model_obj,
                                     field_names=NULL,
                                     session_name=NULL,
                                     alpha=0.05,
                                     threshold){
                                     #area.limit=NULL,
                                     #excur_method = c("EB","QC")

  if (!inherits(model_obj, "BayesGLM")) {
    stop(paste0(
      "The model object is of class ",
      paste0(class(model_obj), collapse=", "),
      " but should be of class 'BayesGLM'."
    ))
  }

  #excur_method <- match.arg(excur_method, c("EB","QC"))

  ### Check session_name argument

  #if only one session, analyze that one
  all_sessions <- model_obj$session_names
  n_sess <- length(all_sessions)
  #if only one session available, use that one
  if(n_sess == 1){ session_name <- all_sessions }
  #if session_name=NULL and multiple sessions available, pick first session and return a warning
  if(is.null(session_name)){
    session_name <- all_sessions[1]
    warning(paste0("Using the first session. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
  }
  #check that session_name is valid
  if(!(session_name %in% all_sessions)) stop(paste0('session_name does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
  sess_ind <- which(all_sessions == session_name)

  #check field_names argument
  if(is.null(field_names)) field_names <- model_obj$task_names
  if(any(!(field_names %in% model_obj$task_names))) stop(paste0("Please specify only field names that corresponds to one of the latent fields: ",paste(model_obj$task_names, collapse=', ')))

  #check alpha argument
	if(alpha > 1 | alpha < 0) stop('alpha value must be between 0 and 1, and it is not')

  mesh <- model_obj$mesh
  n_vox <- mesh$n
  #indices of beta vector corresponding to specified session
  inds <- (1:n_vox) + (sess_ind-1)*n_vox

	#loop over latent fields
	excur <- vector('list', length=length(field_names))
	act <- matrix(NA, nrow=n_vox, ncol=length(field_names))
	colnames(act) <- field_names
	for(f in field_names){

		#if(is.null(area.limit)){
			res.exc <- excursions.inla(model_obj$INLA_model_obj, name=f, ind=inds, u=threshold, type='>', alpha=alpha, method="EB")
		#} else {
		#	res.exc <- excursions.inla.no.spurious(model_obj$INLA_model_obj, mesh=mesh, name=f, ind=inds, u=threshold, type='>', method=excur_method, alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
		#}
	  which_f <- which(field_names==f)
		act[,which_f] <- res.exc$E[inds]
    excur[[which_f]] <- res.exc
	}
	result <- list(active=act, excursions_result=excur)


  #compute size of activations
  areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
  areas_act <- apply(act, 2, function(x) sum(areas_all[x==1]))

  result$areas_all <- areas_all
  result$areas_act <- areas_act

	return(result)
}

#' Identification of areas of activation in a General Linear Model using classical methods
#'
#' @param model_obj A \code{BayesGLM} object
#' @param field_inds Indices of tasks for which to identify activations
#' @param session_name (character) The name of the session that should be
#'   examined. If \code{NULL} (default), the average across all sessions is used.
#' @param alpha A significance level to be used in both the family-wise error
#'   rate (FWER) or false discovery rate (FDR) multiple correction methods of
#'   determining significance.
#' @param threshold Activation threshold (e.g. 0.01 for 1 percent signal change)
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
#' @export
id_activations.classical <- function(model_obj,
                                     field_inds = NULL,
                                     session_name = NULL,
                                     alpha = 0.05,
                                     threshold = 0,
                                     correction = c("FWER","FDR","none"),
                                     mesh = NULL) {

  if (!inherits(model_obj, "classicalGLM")) {
    stop(paste0(
      "The model object is of class ",
      paste0(class(model_obj), collapse=", "),
      " but should be of class 'classicalGLM'."
    ))
  }

  correction <- match.arg(correction, c("FWER","FDR","none"))

  if(is.null(threshold)) threshold <- 0

  ### Check session_name argument

  #if only one session, analyze that one
  all_sessions <- names(model_obj)
  n_sess <- length(all_sessions)
  #if only one session available, use that one
  if(n_sess == 1){ session_name <- all_sessions }
  #if session_name=NULL and multiple sessions available, pick first session and return a warning
  if(is.null(session_name)){
    session_name <- all_sessions[1]
    warning(paste0("Using the first session instead. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
  }
  #check that session_name is valid
  if(!(session_name %in% all_sessions)) stop(paste0('session_name does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
  sess_ind <- which(all_sessions == session_name)

  # #check field_inds argument
  num_fields <- ncol(model_obj[[sess_ind]]$estimates)
  if(is.null(field_inds)){
    field_inds <- seq(num_fields)
    if(num_fields > 1) message(paste0('Since field_inds=NULL, I will analyze all ', num_fields, ' tasks.'))
  } else {
    if(any(!(field_inds %in% seq(num_fields)))) stop(paste0("Please specify only field inds between 1 and ",num_fields))
  }

  beta_est <- model_obj[[sess_ind]]$estimates
  se_beta <- model_obj[[sess_ind]]$SE_estimates
  DOF <- model_obj[[sess_ind]]$DOF

  nvox <- nrow(beta_est)
  #if(any(!(field_inds %in% 1:K))) stop(paste0('field_inds must be between 1 and the number of tasks, ',K))
  beta_est <- matrix(beta_est[,field_inds], nrow=nvox) #need matrix() in case beta_est[,field_inds] is a vector
  se_beta <- matrix(se_beta[,field_inds], nrow=nvox) #need matrix() in case beta_est[,field_inds] is a vector
  K <- length(field_inds)

  #Compute t-statistics and p-values
  t_star <- (beta_est - threshold) / se_beta
  if(!is.matrix(t_star)) t_star <- matrix(t_star, nrow=nvox)
  #perform multiple comparisons correction
  p_values <- p_values_adj <- active <- matrix(NA, nvox, K)

    for(k in 1:K){
    p_values_k <- sapply(t_star[,k], pt, df = DOF, lower.tail = F)
    if(correction == "FWER") p_vals_adj_k <- p.adjust(p_values_k, method='bonferroni')
    if(correction == "FDR") p_vals_adj_k <- p.adjust(p_values_k, method='BH')
    if(correction == "none") p_vals_adj_k <- p_values_k

    p_values[,k] <- p_values_k
    p_values_adj[,k] <- p_vals_adj_k
    active[,k] <- (p_vals_adj_k < alpha)
  }

  #compute size of activations
  if(!is.null(mesh)){
    mask <- model_obj[[sess_ind]]$mask
    if(sum(mask) != nrow(mesh$loc)) stop('Supplied mesh is not consistent with mask in model_obj.')
    areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
    areas_act <- apply(active[mask==1,,drop=FALSE], 2, function(x) sum(areas_all[x==1]))
  } else {
    areas_all <- areas_act <- NULL
  }

  na_pvalues <- which(is.na(p_values[,1]))
  p_values <- p_values[-na_pvalues,, drop = F]
  p_values_adj <- p_values_adj[-na_pvalues,, drop = F]
  active <- active[-na_pvalues,, drop = F]

  result <- list(p_values = p_values,
                 p_values_adj = p_values_adj,
                 active = active,
                 correction = correction,
                 alpha = alpha,
                 threshold = threshold,
                 areas_all = areas_all,
                 areas_act = areas_act)
  return(result)
}


#' Identify activations using joint posterior probabilities with EM results
#'
#' Identifies areas of activation given an activation threshold and significance
#'  level using joint posterior probabilities
#'
#' For a given latent field, identifies locations that exceed a certain activation
#'  threshold (e.g. 1 percent signal change) at a given significance level, based on the joint
#'  posterior distribution of the latent field.
#'
#' @param model_obj An object of class ‘"BayesGLM"’, a result of a call to BayesGLMEM
#' @param field_names Name of latent field or vector of names on which to identify activations.  By default, analyze all tasks.
#' @param session_name (character) The name of the session that should be examined.
#' If \code{NULL} (default), the average across all sessions is used.
#' @param alpha Significance level (e.g. 0.05)
#' @param threshold Activation threshold (e.g. 1 for 1 percent signal change if scale=TRUE in model estimation)
#' @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
#'
#'
#' @return A list with two elements: \code{active}, which gives a matrix of zeros
#' and ones of the same dimension as \code{model_obj$task_estimates${session_name}},
#' and \code{excur_result}, an object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for
#'  more information).
#'
#' @importFrom excursions excursions
#' @importFrom stats na.omit
#'
#' @export
id_activations.em <-
  function(model_obj,
           field_names = NULL,
           session_name = NULL,
           alpha = 0.05,
           threshold,
           area.limit = NULL) {
    if (!inherits(model_obj, "BayesGLM"))
      stop(paste0(
        "The model object is of class ",
        class(model_obj),
        " but should be of class 'BayesGLM'."
      ))
    #check session_name argument

    #if only one session, analyze that one
    all_sessions <- model_obj$session_names
    n_sess <- length(all_sessions)
    if (n_sess == 1) {
      session_name <- all_sessions
    }

    #check session_name not NULL, check that a valid session name
    if(!is.null(session_name)){
      if(!(session_name %in% all_sessions)) stop(paste0('session_name does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
      sess_ind <- which(all_sessions == session_name)
    }

    #if averages not available and session_name=NULL, pick first session and return a warning
    # has_avg <- is.matrix(model_obj$avg_task_estimates)
    has_avg <- !is.null(model_obj$avg_task_estimates)
    if(is.null(session_name) & !has_avg){
      session_name <- all_sessions[1]
      warning(paste0("Your model object does not have averaged beta estimates. Using the first session instead. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
    }

    if(is.null(session_name) & has_avg){
      session_name <- 'avg'
    }

    # Make an indicator for subcortical data
    is_subcort <- "BayesGLMEM_vol3D" %in% as.character(model_obj$call)

    #if session_name is still NULL, use average and check excur_method argument
    # if(is.null(session_name)){
    #   if(has_avg & excur_method != "EB") {
    #     excur_method <- 'EB'
    #     warning("To id activations for averaged beta estimates, only the excur_method='EB' is supported. Setting excur_method to 'EB'.")
    #   }
    # }

    #check field_names argument
    if(is.null(field_names)) field_names <- model_obj$task_names
    if(!any(field_names %in% model_obj$task_names)) stop(paste0("Please specify only field names that corresponds to one of the latent fields: ",paste(model_obj$task_names, collapse=', ')))

    #check alpha argument
    if(alpha > 1 | alpha < 0) stop('alpha value must be between 0 and 1, and it is not')

    mesh <- model_obj$mesh
    if(is_subcort) {
      n_subcort_models <- length(model_obj$spde_obj)
      result <- vector("list", length = n_subcort_models)
      for(m in 1:n_subcort_models){
        # mesh <- make_mesh(model_obj$spde_obj[[m]]$vertices[[1]],
        #                   model_obj$spde_obj[[m]]$faces[[1]])
        # n_vox <- mesh$n
        # excur <- vector('list', length=length(field_names))
        # act <- matrix(NA, nrow=n_vox, ncol=length(field_names))
        # colnames(act) <- field_names

        if(session_name == "avg") {
          beta_est <- c(model_obj$EM_result_all[[m]]$posterior_mu)
        }
        if(session_name != "avg") {
          beta_mesh_est <- c(model_obj$task_estimates[[which(all_sessions == session_name)]])
          beta_est <- beta_mesh_est[!is.na(beta_mesh_est)]
        }
        Phi_k <- model_obj$EM_result_all[[m]]$mesh$Amat
        # Phi <-
        #   Matrix::bdiag(rep(
        #     list(Phi_k),
        #     length(field_names)
        #   ))
        # w <- beta_est %*% Phi
        Sig_inv <- model_obj$EM_result_all[[m]]$posterior_Sig_inv
        # Q_est <- Matrix::tcrossprod(Phi %*% Sig_inv, Phi)

        ### use the ind argument for excursions to get the marginal posterior
        ### excursions sets (much faster)

        n_vox <- length(beta_est) / length(field_names)
        act <- matrix(NA, nrow=n_vox, ncol=length(field_names))

        V <- Matrix::diag(INLA::inla.qinv(Sig_inv))

        for(f in field_names) {
          which_f <- which(field_names==f)
          f_inds <- (1:n_vox) + (which_f - 1)*n_vox
          res.exc <-
            excursions(
              alpha = alpha,
              u = threshold,
              mu = beta_est,
              Q = Sig_inv,
              vars = V,
              ind = f_inds,
              type = ">", method = "EB"
            )

          act[,which_f] <- res.exc$E[f_inds]
        }
        act <- as.matrix(Phi_k %*% act)
        result[[m]] <- list(active = act, excursions_result=res.exc)
      }
    }

    if(!is_subcort) {
      n_vox <- mesh$n

      #for a specific session
      if(!is.null(session_name)){

        # inds <- (1:n_vox) + (sess_ind-1)*n_vox #indices of beta vector corresponding to session v

        #loop over latent fields
        excur <- vector('list', length=length(field_names))
        act <- matrix(NA, nrow=n_vox, ncol=length(field_names))
        colnames(act) <- field_names
        res.exc <-
          excursions(
            alpha = alpha,
            u = threshold,
            mu = stats::na.omit(
              c(model_obj$task_estimates[[session_name]])
            ),
            Q = model_obj$posterior_Sig_inv,
            type = ">", method = "EB"
          )
        for(f in field_names) {
          which_f <- which(field_names==f)
          f_inds <- (1:n_vox) + (which_f - 1)*n_vox
          act[,which_f] <- res.exc$E[f_inds]
        }
      }
      result <- list(active=act, excursions_result=res.exc)
    }

    return(result)
  }
