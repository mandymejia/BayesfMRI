#' Identify activations using 2-means clustering methods
#'
#' @param model_obj An object of class \code{BayesGLM}
#' @param field_name Name of latent field or vector of names on which to identify
#'  activations
#' @param type_2means A string that should be either "point" or "sequential". The
#'  "point" type_2means does a simple 2-means clustering to determine areas of activation.
#'  The "sequential" type_2means uses the sequential 2-means variable selection method,
#'  as described in Li and Pati (2017). The "sequential" method takes significantly
#'  longer, but should do a better job of accounting for posterior variance.
#' @param n_sample The number of samples to generate if the sequential 2-means
#'  type_2means is chosen. By default, this takes a value of 1000.
#'
#' @return A nested list, where the first layer separates by session, and the
#'  second layer is another list of two elements: \code{active}, which gives a
#'  matrix of zeros and ones of the same dimension as
#'  \code{model_obj$beta_estimates${session_name}}, and \code{excur_result},
#'  which is \code{NULL} for the "2means" method.
#'
#' @importFrom stats kmeans dist
#' @importFrom INLA inla.posterior.sample
#'
#' @export
#' @md
id_activations.2means <- function(
  model_obj,
  field_name = NULL,
  type_2means = "point",
  n_sample = NULL) {

  if (!type_2means %in% c("point", "sequential"))
    stop("The type_2means needs to be either 'point' or 'sequential'.")
  if (type_2means == "point") {
    out <- sapply(model_obj$beta_estimates, function(est_type) {
      if (is.null(field_name))
        field_name <- colnames(est_type)
      if (!any(field_name %in% colnames(est_type)))
        stop(
          "Please specify a field name that corresponds to one of the output latent field estimate names (i.e. bbeta1)."
        )
      est_type <-
        est_type[, which(colnames(est_type) %in% field_name)]
      out2 <-
        sapply(split(est_type, col(est_type)), function(beta_est) {
          vector_beta <- c(beta_est)
          if (any(is.na(vector_beta)))
            vector_beta <- vector_beta[!is.na(vector_beta)]
          km_beta <- kmeans(abs(vector_beta), 2)
          which_nonzero <- which.max(km_beta$centers[, 1])
          keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
          return(keep_nonzero)
        }, simplify = TRUE)
      colnames(out2) <- colnames(est_type)
      return(list(active = out2,
                  excur_result = NULL))
    }, simplify = FALSE)
    return(out)
  }
  if (type == "sequential") {
    if (is.null(n_sample)){
      n_sample <- 1000
    }
    n_remain <- n_sample
    if (is.null(field_name))
      field_name <- model_obj$beta_names
    select_vars <- sapply(field_name, function(each_var)
      0,
      simplify = FALSE)
    b_dists <- sapply(model_obj$beta_estimates, function(est_type) {
      if (is.null(field_name))
        field_name <- colnames(est_type)
      if (!any(field_name %in% colnames(est_type)))
        stop(
          "Please specify a field name that corresponds to one of the output latent field estimate names (i.e. bbeta1)."
        )
      est_type <-
        est_type[, which(colnames(est_type) %in% field_name)]
      out2 <-
        sapply(split(est_type, col(est_type)), function(beta_est) {
          vector_beta <- c(beta_est)
          if (any(is.na(vector_beta)))
            vector_beta <- vector_beta[!is.na(vector_beta)]
          km_beta <- kmeans(abs(vector_beta), 2)
          which_nonzero <- which.max(km_beta$centers[, 1])
          nz_beta <- vector_beta[km_beta$cluster == which_nonzero]
          # keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
          # return(keep_nonzero)
          return(dist(range(nz_beta)))
        }, simplify = TRUE)
      names(out2) <- colnames(est_type)
      return(out2)
    }, simplify = FALSE)
    complete_sample <- sapply(model_obj$session_names, function(ses) {
      out <- vector("list", length = length(field_name))
      names(out) <- field_name
      return(out)
    }, simplify = FALSE)
    # complete_sample <- vector("list", length(model_obj$beta_names))
    # names(complete_sample) <- model_obj$beta_names
    while (n_remain > 0) {
      cat("Sampling,", n_remain, "samples remain of", n_sample, "\n")
      next_sample_size <- min(n_remain, 100)
      test_sample <- inla.posterior.sample(
        n = next_sample_size,
        result = model_obj$INLA_result,
        selection = select_vars
      )
      next_sample <-
        sapply(model_obj$session_names, function(each_ses) {
          sapply(model_obj$beta_names, function(each_var) {
            sapply(test_sample, function(each_sample) {
              var_name <- gsub(":[0-9]*", "", rownames(each_sample$latent))
              var_sample <- each_sample$latent[var_name == each_var, ]
              return(var_sample)
            }, simplify = "array")
          }, simplify = FALSE)
        }, simplify = FALSE)

      complete_sample <- mapply(function(complete, next_s) {
        mapply(function(cs, ns) {
          cbind(cs, ns)
        },
        cs = complete,
        ns = next_s,
        SIMPLIFY = FALSE)
      }, complete = complete_sample, next_s = next_sample, SIMPLIFY = FALSE)
      n_remain <- n_remain - next_sample_size
    }
    # b_estimate <- sapply(complete_sample, function(betas) {
    #   median(apply(betas, 1, sd))
    # })
    # b_estimate <- sapply(complete_sample, function(each_ses) {
    #   sapply(each_ses, function(betas) {
    #     # median(apply(betas, 1, sd))
    #     min(apply(betas, 1, sd))
    #   })
    # }, simplify = FALSE)
    b_estimate <- b_dists
    # final_nums <- mapply(s2m_B,
    #                      B = complete_sample,
    #                      sigma = b_estimate,
    #                      SIMPLIFY = TRUE)
    final_nums <- mapply(function(each_ses_beta,each_ses_b) {
      out_id <- mapply(s2m_B,B = each_ses_beta,
                       sigma = each_ses_b, SIMPLIFY = TRUE)
      out_id[out_id != 0] <- 1
      return(list(active = out_id,
                  excur_result = NULL))
    }, each_ses_beta = complete_sample, each_ses_b = b_estimate, SIMPLIFY = FALSE)
    # final_nums[final_nums != 0] <- 1
    # final_out <- list(active = final_nums,
    #                   excur_result = NULL)
    # return(final_out)
    return(final_nums)
  }
}

#' Identify areas of activation
#'
#' For each latent field, identifies locations that are activated
#'
#' @param model_obj An object of class ‘"BayesGLM"’, a result of a call to BayesGLM
# @param method The method to be used for identifying activations, either 'posterior' (Default) or '2means'
#' @param field_names Name of latent field, or vector of names, on which to identify activations
#' @param threshold For method='posterior' only: activation threshold (e.g. 0.01 for 1 percent signal change)
#' @param alpha For method='posterior' only: Significance level (e.g. 0.05)
#' @param area.limit For method='posterior' only: Below this value, activations will be considered spurious.  If NULL (default), no limit.
# @param type_2means For method='2means' only: The type of 2-means clustering to perform ('point' or 'sequential')
# @param nsamp_2means The number of samples to generate if the sequential 2-means type is chosen. By default, 1000.
#'
#' @return A nested list, where the first layer separates by session, and the
#'  second layer is another list of two elements: \code{active}, which gives a
#'  matrix of zeros and ones of the same dimension as
#'  \code{model_obj$beta_estimates${session_name}}, and \code{excur_result}, an
#'  object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for more information)
#  if \code{method='posterior'} and is \code{NULL} if \code{method='2means'}.
#'
#' @export
id_activations <- function(model_obj,
                           #method=c('posterior', '2means'),
                           field_names=NULL,
                           threshold=NULL,
                           alpha=NULL,
                           area.limit=NULL){
  #type_2means=c('point','sequential'),
  #nsamp_2means = NULL){

  if(class(model_obj) != 'BayesGLM') stop('The model_obj argument must be of class BayesGLM, but it is not.')

  #if(length(method)>1) method <- method[1]
  #if(length(type_2means)>1) type_2means <- type_2means[1]

  # if(method=='posterior'){

  message('Identifying activations using joint posterior probabilities')
  if(is.null(threshold)) stop('For identifying activations based on posterior probability method, must specify threshold.')
  if(is.null(alpha)) stop('For identifying activations based on posterior probability method, must specify alpha.')
  #if(!is.null(type_2means)) stop('For identifying activations based on posterior probability method, do not specify type_2means.')
  out <- id_activations.posterior(model_obj=model_obj, field_names=field_names, threshold=threshold, alpha=alpha, area.limit=area.limit)

  # } else if(method=='2means'){
  #
  #   message('Identifying activations using 2-means clustering')
  #   if(is.null(type_2means)) stop('For identifying activations based on 2means, you must specify type_2means.')
  #   if(!is.null(threshold)) stop('For identifying activations based on 2means, do not specify threshold.')
  #   if(!is.null(alpha)) stop('For identifying activations based on posterior probability method, do not specify alpha.')
  #   out <- id_activations.2means(model_obj=model_obj, field_names=field_names, type_2means=type_2means, n_sample=n_sample)
  #
  # } else {
  #   stop('Must specify method = "posterior" or "2means"')
  # }

  return(out)
}

