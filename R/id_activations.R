#' Identify Task Activations
#'
#' Identify areas of activation for each task
#'
#' @param model_obj Result of BayesGLM_cifti model call (of class BayesGLM_cifti)
# @param method The method to be used for identifying activations, either 'posterior' (Default) or '2means'
#' @param field_names Name of latent field, or vector of names, on which to identify activations. By default, analyze all tasks.
#' @param session_name (character) The name of the session that should be
#' examined. If \code{NULL} (default), the average across all sessions is used.
#' @param alpha Significance level (e.g. 0.05)
#' @param method Either 'Bayesian' or 'classical'
#' @param threshold Activation threshold (e.g. 1 for 1 percent signal change if scale=TRUE during model estimation)
#' @param excur_method For method = 'Bayesian' only: Either \code{EB} (empirical Bayes) or \code{QC} (Quantile
#'   Correction), depending on the method that should be used to find the
#'   excursions set. Note that if any contrasts (including averages across
#'   sessions) are used in the modeling, the method chosen must be \code{EB}.
#'   The difference in the methods is that the \code{EB} method assumes Gaussian
#'    posterior distributions for the parameters.
#' @param area.limit For method = 'Bayesian' only: Below this value, clusters of activations will be considered spurious.  If NULL (default), no limit.
#' @param correction For method = 'classical' only: Type of multiple comparisons correction, 'FDR' (Benjamini Hochberg) or 'FWER' (Bonferroni correction)
#' @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
#' depending on the method that should be used to find the excursions set. Note that to ID
#' activations for averages across sessions, the method chosen must be \code{EB}. The difference
#' in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#' @param verbose (Logical) If TRUE, print progress updates.
# @param type For method='2means' only: The type of 2-means clustering to perform ('point' or 'sequential')
# @param n_sample The number of samples to generate if the sequential 2-means type is chosen. By default, this takes a value of 1000.
#'
# @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @importFrom ciftiTools convert_xifti
#' @export
#'
#'
id_activations_cifti <- function(model_obj,
                                 field_names=NULL,
                                 session_name=NULL,
                                 alpha=0.05,
                                 method=c('Bayesian','classical'),
                                 threshold=NULL,
                                 area.limit=NULL,
                                 correction = c("FWER","FDR","permutation","none"),
                                 excur_method = c("EB","QC"),
                                 verbose = TRUE){

  if(class(model_obj) != 'BayesGLM_cifti') stop('The model_obj argument must be of class BayesGLM_cifti (output of BayesGLM_cifti function), but it is not.')

  method <- match.arg(method, c('Bayesian','classical'))
  if(!(method %in% c('Bayesian','classical'))) stop("The method argument should only be 'Bayesian' or 'classical'.")

  if(is.null(threshold)){
    if(method=='classical') threshold <- 0
    if(method=='Bayesian') stop("Must specify an activation threshold when method='Bayesian'.")
  }

  if(method=='Bayesian') GLM_list <- model_obj$GLMs_Bayesian
  if(method=='classical') GLM_list <- model_obj$GLMs_classical
  num_models <- length(GLM_list)
  activations <- vector('list', length=num_models)
  names(activations) <- names(GLM_list)
  do_left <- !is.null(GLM_list$cortexL)
  do_right <- !is.null(GLM_list$cortexR)

  if(is.null(field_names)) field_names <- model_obj$beta_names
  if(any(!(field_names %in% model_obj$beta_names))) stop(paste0('All elements of field_names must appear in model_obj$beta_names: ', paste(model_obj$beta_names, collapse=',')))
  field_inds <- which(model_obj$beta_names %in% field_names)

  if(method=='Bayesian'){
    for(mm in 1:num_models){
      if(is.null(GLM_list[[mm]])) next
      model_m <- GLM_list[[mm]]
      if(verbose) cat(paste0('Identifying Bayesian GLM activations in ',names(GLM_list)[mm],'\n'))
      act_m <- id_activations.posterior(model_obj=model_m,
                                                    field_names=field_names,
                                                    session_name=session_name,
                                                    threshold=threshold,
                                                    alpha=alpha,
                                                    area.limit=area.limit,
                                                    excur_method=excur_method)

      activations[[mm]] <- act_m
    }
  }

  if(method=='classical'){
    #if(verbose) cat('Identifying classical GLM activations')
    for(mm in 1:num_models){
      if(is.null(GLM_list[[mm]])) next
      model_m <- GLM_list[[mm]]
      act_m <- id_activations.classical(model_obj=model_m,
                                        field_inds=field_inds,
                                        session_name=session_name,
                                        threshold=threshold,
                                        alpha=alpha,
                                        correction=correction)

      activations[[mm]] <- act_m
    }
  }

  #map results to xifti objects
  activations_xifti <- 0*model_obj[[paste0("betas_", method[1])]][[1]]
  if(length(field_names) != length(model_obj$beta_names)) activations_xifti$meta$cifti$names <- field_names
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

  activations_xifti <- convert_xifti(activations_xifti, to="dlabel", colors='red')

  result <- list(activations = activations,
                 activations_xifti = activations_xifti)

  return(result)
}







  # model_result <- result$model_result
  # active_xifti <- clear_data(result$subjICmean_xifti)
  # nleft <- nrow(result$subjICmean_xifti$data$cortex_left)
  # nright <- nrow(result$subjICmean_xifti$data$cortex_right)
  #
  # if(class(model_result)=='stICA' & is.null(spatial_model)) spatial_model <- TRUE
  # if(class(model_result)=='tICA' & is.null(spatial_model)) spatial_model <- FALSE
  # if((spatial_model==TRUE) & (class(model_result) == 'tICA')) {
  #   warning('spatial_model set to TRUE but class of model result is tICA. Setting spatial_model = FALSE, performing inference using standard template ICA.')
  #   spatial_model <- FALSE
  # }
  #
  # #if spatial model available but spatial_model set to FALSE, grab standard template ICA model result
  # if(class(model_result)=='stICA' & spatial_model==FALSE){
  #   model_result <- model_result$result_tICA
  # }
  #
  # #run activations function
  # activations_result <- activations(model_result, u=u, alpha=alpha, type=type, method_p=method_p, verbose=verbose, which.ICs=which.ICs, deviation=deviation)
  #
  # #construct xifti object for activation maps
  # act_viz <- activations_result$active*1
  # act_viz[act_viz==0] <- NA
  # active_xifti$data$cortex_left <- act_viz[1:nleft,]
  # active_xifti$data$cortex_right <- act_viz[nleft+(1:nright),]
  #
  # #include convert_xifti function
  #
  # return(list(activations = activations_result,
  #             active_xifti = active_xifti))





## HERE: Update this function to use field_inds instead of field_names

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
#' @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
#' @param excur_method Either \code{"EB"} (empirical Bayes) or \code{"QC"} (Quantile Correction),
#' depending on the method that should be used to find the excursions set. Note that to ID
#' activations for averages across sessions, the method chosen must be \code{EB}. The difference
#' in the methods is that the \code{EB} method assumes Gaussian posterior distributions for the parameters.
#'
#'
#' @return A list with two elements: \code{active}, which gives a matrix of zeros
#' and ones of the same dimension as \code{model_obj$beta_estimates${session_name}},
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
                                     threshold,
                                     area.limit=NULL,
                                     excur_method = c("EB","QC")){

  if(class(model_obj) != "BayesGLM") stop(paste0("The model object is of class ",class(model_obj)," but should be of class 'BayesGLM'."))

  excur_method <- match.arg(excur_method, c("EB","QC"))

  #check session_name argument

  #if only one session, analyze that one
  all_sessions <- model_obj$session_names
  n_sess <- length(all_sessions)
  if(n_sess == 1){ session_name <- all_sessions }

  #if averages not available and session_name=NULL, pick first session and return a warning
  has_avg <- is.matrix(model_obj$avg_beta_estimates)
  if(is.null(session_name) & !has_avg){
    session_name <- all_sessions[1]
    warning(paste0("Your model object does not have averaged beta estimates. Using the first session instead. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
  }

  #if session_name is still NULL, use average and check excur_method argument
  if(is.null(session_name)){
    if(has_avg & excur_method != "EB") {
      excur_method <- 'EB'
      warning("To id activations for averaged beta estimates, only the excur_method='EB' is supported. Setting excur_method to 'EB'.")
    }
  }

  #check session_name not NULL, check that a valid session name
  if(!is.null(session_name)){
    if(!(session_name %in% all_sessions)) stop(paste0('session_name does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
    sess_ind <- which(all_sessions == session_name)
  }

  #check field_names argument
  if(is.null(field_names)) field_names <- model_obj$beta_names
  if(!any(field_names %in% model_obj$beta_names)) stop(paste0("Please specify only field names that corresponds to one of the latent fields: ",paste(model_obj$beta_names, collapse=', ')))

  #check alpha argument
	if(alpha > 1 | alpha < 0) stop('alpha value must be between 0 and 1, and it is not')

  mesh <- model_obj$mesh
  n_vox <- mesh$n

	#for a specific session
	if(!is.null(session_name)){

	  inds <- (1:n_vox) + (sess_ind-1)*n_vox #indices of beta vector corresponding to session v

		#loop over latent fields
		excur <- vector('list', length=length(field_names))
		act <- matrix(NA, nrow=n_vox, ncol=length(field_names))
		colnames(act) <- field_names
		for(f in field_names){

  		if(is.null(area.limit)){
  			res.exc <- excursions.inla(model_obj$INLA_result, name=f, ind=inds, u=threshold, type='>', method=excur_method, alpha=alpha)
  		} else {
  			res.exc <- excursions.inla.no.spurious(model_obj$INLA_result, mesh=mesh, name=f, ind=inds, u=threshold, type='>', method=excur_method, alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
  		}
		  which_f <- which(field_names==f)
  		act[,which_f] <- res.exc$E[inds]
      excur[[which_f]] <- res.exc
		}
		result <- list(active=act, excursions_result=excur)
	}

	#for the average over sessions
	if(is.null(session_name)){

	  avg_inds <- model_obj$INLA_result$misc$configs$config[[1]]$pred_idx
	  mu_avg <- model_obj$INLA_result$misc$configs$config[[1]]$mean[avg_inds]
	  Q_avg <- model_obj$INLA_result$misc$configs$config[[1]]$Q[avg_inds,avg_inds]
	  avg_exc <- excursions(alpha = alpha,
	                        u = threshold,
	                        mu = mu_avg,
	                        Q = Q_avg,
	                        type = ">",
	                        method = "EB")
	  act <- matrix(NA, nrow=n_vox, ncol=length(field_names))
	  for(f in field_names){
	    which_f <- which(field_names==f)
	    f_inds <- (1:n_vox) + (which_f - 1)*n_vox
	    act[,which_f] <- avg_exc$E[f_inds]
	  }

	  result <- list(
	    active = act,
	    excursions_result = avg_exc
	  )
	}

  #compute size of activations
  areas_all <- diag(inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
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
                                     correction = c("FWER","FDR","permutation","none"),
                                     mesh = NULL) {
  # Bring in data for debugging
  # model_obj <- readRDS("~/Desktop/id_activations_classical_input.rds")
  # field_inds = NULL
  # session_name = NULL
  # alpha = 0.05
  # threshold = 0
  # correction = 'permutation'
  # mesh = NULL
  # Argument check
  if(class(model_obj) != "classicalGLM") stop(paste0("The model object is of class ",class(model_obj)," but should be of class 'classicalGLM'."))

  correction <- match.arg(correction, c("FWER","FDR","permutation","none"))

  if(is.null(threshold)) threshold <- 0

  #check session_name argument

  #if only one session, analyze that one
  all_sessions <- names(model_obj)
  n_sess <- length(all_sessions)
  if(n_sess == 1){ session_name <- all_sessions }

  # #if averages not available and session_name=NULL, pick first session and return a warning
  has_avg <- ('avg' %in% names(model_obj))
  if(is.null(session_name) & !has_avg){
    session_name <- all_sessions[1]
    warning(paste0("Your model object does not have averaged beta estimates. Using the first session instead. For a different session, specify a session name from among: ", paste(all_sessions, collapse = ', ')))
  }

  # If averages are available and no session name is defined, use the averages
  if(is.null(session_name) & has_avg){
    session_name <- 'avg'
    message("No session name defined. Using the average of sessions.")
    sess_ind <- which(all_sessions == session_name)
  }

  #check session_name not NULL, check that a valid session name
  if(!is.null(session_name)){
    if(!(session_name %in% all_sessions)) stop(paste0('session_name does not appear in the list of sessions: ', paste(all_sessions, collapse=', ')))
    sess_ind <- which(all_sessions == session_name)
  }

  # #check field_inds argument
  if(is.null(field_inds)){
    num_fields <- ncol(model_obj[[sess_ind]]$estimates)
    field_inds <- seq(num_fields)
    if(num_fields > 1) message(paste0('Since field_inds=NULL, I will analyze all ', num_fields, ' tasks.'))
  }

  beta_est <- model_obj[[sess_ind]]$estimates
  se_beta <- model_obj[[sess_ind]]$SE_estimates
  DOF <- model_obj[[sess_ind]]$DOF

  nvox <- nrow(beta_est)
  K <- ncol(beta_est)
  if(any(!(field_inds %in% 1:K))) stop(paste0('field_inds must be between 1 and the number of tasks, ',K))
  beta_est <- matrix(beta_est[,field_inds], nrow=nvox) #need matrix() in case beta_est[,field_inds] is a vector
  se_beta <- matrix(se_beta[,field_inds], nrow=nvox) #need matrix() in case beta_est[,field_inds] is a vector
  K <- length(field_inds)

  #Compute t-statistics and p-values
  t_star <- (beta_est - threshold) / se_beta
  if(!is.matrix(t_star)) t_star <- matrix(t_star, nrow=nvox)
  #perform multiple comparisons correction
  if(correction != "permutation") {
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
  na_pvalues <- which(is.na(p_values[,1]))
  p_values <- p_values[-na_pvalues,, drop = F]
  p_values_adj <- p_values_adj[-na_pvalues,, drop = F]
  active <- active[-na_pvalues,, drop = F]
  }

  if(correction == "permutation") {
    if(sum(c("null_estimates","null_SE_estimates") %in% names(model_obj[[sess_ind]])) != 2) stop("This model_obj is from an older version and is not compatible with permutation testing.")
    if(is.null(model_obj[[sess_ind]]$null_estimates) | is.null(model_obj[[sess_ind]]$null_SE_estimates)) stop("Permutations were not completed in the analysis that created the model_obj. Please choose a different correction method.")
    p_values <- p_values_adj <- NULL
    t_stats <- (model_obj[[sess_ind]]$null_estimates - threshold) / model_obj[[sess_ind]]$null_SE_estimates
    max_tstats <- apply(t_stats,2:3,max, na.rm = TRUE)
    null_thresholds <- apply(max_tstats, 1, quantile, probs = (1 - alpha))
    active <- mapply(function(ts,null_ts) {
      return(ts > null_ts)
    }, ts = split(t_star, col(t_star)), null_ts = null_thresholds)
    if(!"matrix" %in% class(active)) active <- as.matrix(active)
    active <- active[!is.na(active[,1]),]
  }

  #compute size of activations
  if(!is.null(mesh)){
    mask <- model_obj[[1]]$mask
    if(sum(mask) != nrow(mesh$loc)) stop('Supplied mesh is not consistent with mask in model_obj.')
    areas_all <- diag(inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
    areas_act <- apply(active[mask==1,], 2, function(x) sum(areas_all[x==1]))
  } else {
    areas_all <- areas_act <- NULL
  }



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
