#' Make Mesh
#'
#' Make triangular mesh from faces and vertices.
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#'
#' @return Triangular mesh from matrices and vertices
#'
#' @export
make_mesh <- function(vertices, faces){

  check_INLA(FALSE)

  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  # Construct mesh
  mesh <- INLA::inla.mesh.create(loc = as.matrix(vertices), tv = as.matrix(faces))
  return(mesh)
}


