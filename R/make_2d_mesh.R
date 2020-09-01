#' Make a 2-D Mesh
#'
#' @param mask A binary (0/1) matrix that should be used to make a mesh for two-dimensional data
#'
#' @return  A \code{inla.mesh} object that can be used to run the \code{BayesGLM_surface} function with two-dimensional data
#' @export
make_2d_mesh <- function(mask) {
  in_mask <- which(mask == 1, arr.ind = T)
  in_mask <- in_mask[,2:1]
  boundary <- INLA::inla.nonconvex.hull(in_mask, resolution = 100)
  mesh <- INLA::inla.mesh.2d(loc = in_mask, boundary = boundary, max.edge = c(2,4))
}
