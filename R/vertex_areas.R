#' Surface area of each vertex
#' 
#' Compute surface areas of each vertex in a triangular mesh.
#' 
#' @inheritSection INLA_Description INLA Requirement
#' 
#' @inheritParams mesh_Param_inla
#'
#' @return Vector of areas
#' 
#' @export
vertex_areas <- function(mesh) {
  if(missing(mesh)) { stop("`mesh` input is required.")}

  if (!inherits(mesh, "inla.mesh")) {
    stop("`mesh` needs to be of class `'inla.mesh'`.")
  }

  diag(INLA::inla.fmesher.smorg(
    mesh$loc,mesh$graph$tv, fem = 0, output = list("c0")
  )$c0)
}
