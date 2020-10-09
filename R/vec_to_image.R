#' Transform vector data to an image
#'
#' This fills in parts of a template with values from \code{vec_data}.
#'
#' @param vec_data A V by p matrix, where V is the number of voxels within a
#'   mask and p is the number of vectors to transform into matrix images
#' @param template_image A binary matrix in which V entries are 1 and the rest
#'   of the entries are zero
#'
#' @return A list of masked values from \code{vec_data}
#' @export
vec2image <- function(vec_data, template_image) {
  each_col <- sapply(split(vec_data, col(vec_data)), function(vd) {
    out <- template_image
    out[out == 1] <- vd
    out[out == 0] <- NA
    return(out)
  }, simplify = F)
  return(each_col)
}
