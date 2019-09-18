#' Identifies areas of activation
#'
#' @description For each latent field, identifies locations that are activated
#'
#' @param model_obj An object of class ‘"BayesGLM"’, a result of a call to BayesGLM
#' @param method The method to be used for identifying activations, either 'posterior' (Default) or '2means'
#' @param field_name Name of latent field, or vector of names, on which to identify activations
#' @param session_names Names of sessions included in INLA model that resulted in object
#' @param threshold For method='posterior' only: activation threshold (e.g. 0.01 for 1 percent signal change)
#' @param alpha For method='posterior' only: Significance level (e.g. 0.05)
#' @param area.limit For method='posterior' only: Below this value, activations will be considered spurious.  If NULL (default), no limit.
#' @param type For method='2means' only: The type of 2-means clustering to perform ('point' or 'sequential')
#'
#' @details Put additional details here.
#'
#' @return A nested list, where the first layer separates by session, and the second layer is another list of two elements: `active`, which gives a matrix of zeros and ones of the same dimention as `model_obj$beta_estimates${session_name}`, and `excur_result`, an object of class excurobj if method='posterior' (see `help(excursions.inla)` for more information) and is NULL if method='2means'.
#'
#' @examples \dontrun{}
id_activations <- function(model_obj, method=c('posterior', '2means'), field_name=NULL, threshold=NULL, alpha=NULL, area.limit=NULL, type=c('point','sequential')){

  if(class(model_obj) != 'BayesGLM') stop('The model_obj argument must be of class BayesGLM, but it is not.')

  if(length(method)>1) method <- method[1]
  if(length(type)>1) type <- type[1]

  if(method=='posterior'){

    message('Identifying activations using joint posterior probabilities')
    if(is.null(threshold)) stop('For identifying activations based on posterior probability method, must specify threshold.')
    if(is.null(alpha)) stop('For identifying activations based on posterior probability method, must specify alpha.')
    if(!is.null(type)) stop('For identifying activations based on posterior probability method, do not specify type.')
    out <- id_activations.posterior(model_obj=model_obj, field_name=field_name, threshold=threshold, alpha=alpha, area.limit=area.limit)

  } else if(method=='2means'){

    message('Identifying activations using 2-means clustering')
    if(is.null(type)) stop('For identifying activations based on 2means, you must specify type.')
    if(!is.null(threshold)) stop('For identifying activations based on 2means, do not specify threshold.')
    if(!is.null(alpha)) stop('For identifying activations based on posterior probability method, do not specify alpha.')
    out <- id_activations.2means(model_obj=model_obj, field_name=field_name, type=type)

  } else {
    stop('Must specify method = "posterior" or "2means"')
  }

  return(out)
}

#' Identifies areas of activation given an activation threshold and significance level using joint posterior probabilities
#'
#' @description For a given latent field, identifies locations that exceed a certain activation
#' threshold (e.g. 1 percent signal change) at a given significance level, based on the joint
#' posterior distribution of the latent field.
#'
#' @param model_obj An object of class ‘"BayesGLM"’, a result of a call to BayesGLM
#' @param field_name Name of latent field or vector of names on which to identify activations
#' @param threshold Activation threshold (e.g. 0.01 for 1 percent signal change)
#' @param alpha Significance level (e.g. 0.05)
#' @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
#'
#' @details Put additional details here.
#'
#' @return A nested list, where the first layer separates by session, and the second layer is another list of two elements: `active`, which gives a matrix of zeros and ones of the same dimention as `model_obj$beta_estimates${session_name}`, and `excur_result`, an object of class excurobj (see `help(excursions.inla)` for more information).
#'
#' @examples \dontrun{}
id_activations.posterior <- function(model_obj, field_name=NULL, threshold, alpha=0.05, area.limit=NULL){

  session_names <- model_obj$session_names
	n_sess <- length(session_names)
	n_vox <- model_obj$mesh$n

	if(is.null(field_name)) field_name <- model_obj$beta_names
	if(!any(field_name %in% model_obj$beta_names)) stop("Please specify only field names that corresponds to one of the latent fields (i.e. 'bbeta1').")

	if(alpha > 1 | alpha < 0) stop('alpha value must be between 0 and 1, and it is not')

	result <- vector('list', n_sess)
	names(result) <- session_names
	for(v in 1:n_sess){
		sess_name <- session_names[v]
		inds_v <- (1:n_vox) + (v-1)*n_vox #indices of beta vector corresponding to session v

		#loop over latent fields
		excur_v <- vector('list', length=length(field_name))
		act_v <- matrix(NA, nrow=n_vox, ncol=length(field_name))
		colnames(act_v) <- field_name
		for(f in field_name){

  		if(is.null(area.limit)){
  			res.exc <- excursions.inla(model_obj$INLA_result, name=f, ind=inds_v, u=threshold, type='>', method='QC', alpha=alpha, F.limit=0.2)
  		} else {
  			res.exc <- excursions.inla.no.spurious(model_obj$INLA_result, mesh=mesh, name=f, ind=inds_v, u=threshold, type='>', method='QC', alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
  		}
		  which_f <- which(field_name==f)
  		act_v[,which_f] <- res.exc$E[inds_v]
  		#act_v[is.na(act_v),which_f] <- 0
      excur_v[[which_f]] <- res.exc
		}
		result[[v]] <- list(active=act_v, excursions_result=excur_v)
	}

	out <- result
	return(out)
}



#' Identify activations using 2-means clustering methods
#'
#' @param model_obj An object of class `BayesGLM`
#' @param field_name Name of latent field or vector of names on which to identify activations
#' @param type A string that should be either "point" or "sequential". The "point" type does a simple 2-means clustering to determine areas of activation. The "sequential" type uses the sequential 2-means variable selection method, as described in Li and Pati (2017). The "sequential" method takes significantly longer, but should do a better job of accounting for posterior variance.
#' @param n_sample The number of samples to generate if the sequential 2-means type is chosen. By default, this takes a value of 1000.
#'
#' @return A nested list, where the first layer separates by session, and the second layer is another list of two elements: `active`, which gives a matrix of zeros and ones of the same dimention as `model_obj$beta_estimates${session_name}`, and `excur_result`, which is NULL for the 2means method.
#' @export
#'
#' @examples \dontrun{}
#' @md
id_activation.2means <- function(model_obj, field_name = NULL, type = "point", n_sample = NULL) {
  if(!type %in% c("point","sequential")) stop("The type needs to be either 'point' or 'sequential'.")
  if(type == "point") {
    out <- sapply(model_obj$beta_estimates, function(est_type) {
      if(is.null(field_name)) field_name <- colnames(est_type)
      if(!any(field_name %in% colnames(est_type))) stop("Please specify a field name that corresponds to one of the output latent field estimate names (i.e. bbeta1).")
      est_type <- est_type[,which(colnames(est_type) %in% field_name)]
      out2 <- sapply(split(est_type,col(est_type)), function(beta_est) {
        vector_beta <- c(beta_est)
        if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
        km_beta <- kmeans(abs(vector_beta),2)
        which_nonzero <- which.max(km_beta$centers[,1])
        keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
        return(keep_nonzero)
      }, simplify = TRUE)
      colnames(out2) <- colnames(est_type)
      return(out2)
    }, simplify = FALSE)
    final_out <- list(active = out,
                      excur_result = NULL)
    return(final_out)
  }
  if(type == "sequential") {
    if(is.null(n_sample)) n_sample <- 1000
    n_remain <- n_sample
    if(is.null(field_names)) field_name <- model_obj$beta_names
    select_vars <- sapply(field_name, function(each_var) 0,
                          simplify = FALSE)
    complete_sample <- vector("list",length(model_obj$beta_names))
    names(complete_sample) <- model_obj$beta_names
    while(n_remain > 0) {
      cat("Sampling,", n_remain, "samples remain of", n_sample, "\n")
      next_sample_size <- min(n_remain,100)
      test_sample <- inla.posterior.sample(n = next_sample_size,
                                           result = model_obj$INLA_result,
                                           selection = select_vars)
      next_sample <- sapply(model_obj$beta_names, function(each_var) {
        sapply(test_sample, function(each_sample) {
          var_name <- gsub(":[0-9]*","",rownames(each_sample$latent))
          var_sample <- each_sample$latent[var_name == each_var,]
          return(var_sample)
        }, simplify = "array")
      },simplify = FALSE)
      complete_sample <- mapply(function(cs,ns) {cbind(cs,ns)},
                                cs = complete_sample, ns = next_sample, SIMPLIFY = FALSE)
      n_remain <- n_remain - next_sample_size
    }
    b_estimate <- sapply(complete_sample, function(betas){
      median(apply(betas,1,sd))
    })
    final_nums <- mapply(s2m_B,B = complete_sample,
                         sigma = b_estimate, SIMPLIFY = TRUE)
    final_nums[final_nums != 0] <- 1
    final_out <- list(active = final_nums,
                      excur_result = NULL)
    return(final_out)
  }
}
