#' Construct a triangular mesh from a 2D volumetric mask (need to generalize this to a 3D volume)
#'
#' @inheritSection INLA_Description INLA Requirement
#' 
#' @param mask A matrix of 0s and 1s representing a volumetric mask
#'
#' @return An inla.mesh object.  See \code{\link{inla.mesh.2d}} for details.
#' 
#' @importFrom INLA inla.nonconvex.hull inla.mesh.2d
#' 
#' @export
#'
vol2mesh <- function(mask){
  xy.in <- which(mask==1, arr.ind=TRUE)[,2:1]
  boundary <- inla.nonconvex.hull(xy.in, resolution = 100)
  mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))
}
