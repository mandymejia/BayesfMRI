#' Find nonzero element in a matrix using 2-means clustering
#'
#' @param beta_est A vector or matrix of values from which values close to zero should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values close to zero are assigned the value of zero. The closeness of a value to zero is found by performing two-means clustering on the absolute values of beta_est, and
#' @export
#'
#' @examples
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}

#' Sequential 2-means variable selection
#'
#' @param x A vector consisting of all variables of interest for a single draw from a posterior distribution
#' @param b A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return The number of nonzero values detected within x
#' @export
#'
#' @examples
s2m <- function(x,b){
  two_means <- kmeans(abs(x),2)
  zero_idx <- which(two_means$cluster == which.min(two_means$centers))
  A <- x[zero_idx]
  two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
  iterations <- 1
  while(abs(two_centers$centers[1, 1] - two_centers$centers[2, 1]) > b) {
    zero_idx <- which(two_centers$cluster == which.min(two_centers$centers))
    A <- A[zero_idx]
    two_centers <- kmeans(abs(A),2,algorithm=c("Lloyd"))
    iterations <- iterations + 1
  }
  num_nonzero <- length(x) - length(A)
  return(num_nonzero)
}

#' Sequential 2-means on array B
#'
#' @param B An array of posterior samples (typically a matrix), in which the last margin corresponds to a single posterior sample
#' @param sigma  A scale parameter used to determine at what distance cluster centers are considered to be the same.
#'
#' @return An array of dimension `head(dim(B),-1)` with a point estimate of B based on the sequential 2-means method
#' @export
#'
#' @examples
#' @md
s2m_B <- function(B,sigma){
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}

#' Identify activations using 2-means clustering methods
#'
#' @param model_obj An object of class `BayesGLM`
#' @param field_name Name of latent field or vector of names on which to identify activations
#' @param type A string that should be either "point" or "sequential". The "point" type does a simple 2-means clustering to determine areas of activation. The "sequential" type uses the sequential 2-means variable selection method, as described in Li and Pati (2017). The "sequential" method takes significantly longer, but should do a better job of accounting for posterior variance.
#' @param n_sample The number of samples to generate if the sequential 2-means type is chosen. By default, this takes a value of 1000.
#'
#' @return A list, where the first element gives the areas of activation in a nested list. The first list layer separates by session, and the second list layer is a list of two elements: `active`, which gives a matrix of zeros and ones of the same dimention as `model_obj$beta_estimates${session_name}`, and `excur_result`, which will take the value `NULL`.
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
