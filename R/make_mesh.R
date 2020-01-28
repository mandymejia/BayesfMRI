#' Make Mesh
#'
#' Make triangular mesh from faces and vertices.
#'
#' @param vertices Matrix of vertices
#' @param faces Matrix of faces
#'
#' @return Triangular mesh from matrices and vertices
#' @export
#' @importFrom INLA inla.mesh.create
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
make_mesh <- function(vertices, faces){

  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  # Construct mesh
  mesh <- inla.mesh.create(loc = as.matrix(vertices), tv = as.matrix(faces))
  return(mesh)
}


