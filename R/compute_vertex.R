#' Compute surface areas of each vertex in a triangular mesh.
#'
#' @param mesh object of type mesh
#'
#' @return Vector of areas
#' @export
#' @importFrom INLA inla.fmesher.smorg
#' @examples

compute_vertex_areas <- function(mesh)
{
  if(!exists(mesh))
  { print("This needs a valid mesh! Input a mesh now")}

  if(class(mesh)==inla.mesh)
    {
    areas <- diag(inla.fmesher.smorg(mesh$loc,mesh$graph$tv, fem = 0, output = list("c0"))$c0)
    return(areas)
  }  else {
    stop("Error in the class of mesh you input.It needs to be inla.mesh")

    }

}
