#' Compute surface areas of each vertex in a triangular mesh.
#' 
#' @inheritSection INLA_Description INLA Requirement
#' 
#' @inheritParams mesh_Param_inla
#'
#' @return Vector of areas
#' 
#' @export
compute_vertex_areas <- function(mesh)
{
  if(missing(mesh))
  { print("This needs a valid mesh! Input a mesh now")}

  if(class(mesh)=="inla.mesh")
    {
    areas <- diag(INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, fem = 0, output = list("c0"))$c0)
    return(areas)
  }  else {
    stop("Error in the class of mesh you input.It needs to be inla.mesh")

    }
}
