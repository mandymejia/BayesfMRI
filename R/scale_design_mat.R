#' Scale the design matrix
#'
#' @param design_mat The original (unscaled) design matrix that is T x K, where
#'     T is the number of time points, and k is the number of task covariates
#'
#' @return A scaled design matrix
#' 
#' @keywords internal
#' 
scale_design_mat <- function(design_mat) {
  stopifnot(is.matrix(design_mat))
  apply(design_mat,2,function(task) {
    (task - mean(task)) / max(abs(task))
  })
}
