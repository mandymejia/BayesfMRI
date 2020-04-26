#' Construct a triangular mesh from a 2D volumetric mask (need to generalize this to a 3D volume)
#'
#' @param mask A matrix of 0s and 1s representing a volumetric mask
#'
#' @return An inla.mesh object.  See `help(inla.mesh.2d)` for details.
#' @export
#' @importFrom INLA inla.nonconvex.hull inla.mesh.2d
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#'
vol2mesh <- function(mask){
  xy.in <- which(mask==1, arr.ind=TRUE)[,2:1]
  boundary <- inla.nonconvex.hull(xy.in, resolution = 100)
  mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))
}
