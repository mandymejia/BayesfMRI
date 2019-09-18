#' Find which values in a matrix appear to be nonzero
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

s2m_B <- function(B,sigma){
  nonzero_nums <- sapply(asplit(B,length(dim(B))),function(B_s) s2m(c(B_s),sigma))
  num_nonzero <- ceiling(median(nonzero_nums))
  median_B <- apply(B,seq(length(dim(B)) - 1),median)
  cutoff <- quantile(c(abs(median_B)),1 - num_nonzero/length(median_B))
  out <- median_B
  out[which(out < cutoff)] <- 0
  return(out)
}

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
