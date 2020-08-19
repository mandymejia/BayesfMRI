#' Scale the design matrix
#'
#' @param design_mat The original (unscaled) design matrix that is T x K, where
#'     T is the number of time points, and k is the number of task covariates
#'
#' @return A scaled design matrix
#' @export
#'
#' @examples
#' library(neuRosim)
#' # Task 1
#' t1 <-
#'   specifydesign(
#'     onsets = seq(0, 200, by = 40),
#'     durations = 1,
#'     totaltime = 200,
#'     TR = 1,
#'     effectsize = 1.3,
#'     conv = "double-gamma",
#'     param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
#'   )
#' # Task 2
#' t2 <-
#'   specifydesign(
#'     onsets = seq(20, 200, by = 40),
#'     durations = 1,
#'     totaltime = 200,
#'     TR = 1,
#'     effectsize = 1.3,
#'     conv = "double-gamma",
#'     param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
#'   )
#' A <- cbind(t1,t2) # This is the design matrix
#' B <- scale_design_mat(A)
scale_design_mat <- function(design_mat) {
  if(!"matrix" %in% class(design_mat)) stop("The design matrix must be a matrix class object.")
  output_mat <- apply(design_mat,2,function(task) {
    returned_col <- task / max(task)
    returned_col <- returned_col - mean(returned_col)
    return(returned_col)
  })
  return(output_mat)
}
