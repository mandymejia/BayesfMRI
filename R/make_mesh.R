#' Make Mesh
#'
#' Make INLA triangular mesh from faces and vertices
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#'
#' @return INLA triangular mesh
#'
#' @export
make_mesh <- function(vertices, faces) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "`make_mesh` requires the `INLA` package. Please install it.", 
      call. = FALSE
    )
  }

  # Check index of faces
  if (min(faces) == 0) { faces <- faces + 1 }

  # Construct mesh
  mesh <- INLA::inla.mesh.create(
    loc = as.matrix(vertices), tv = as.matrix(faces)
  )

  # # Trying to avoid using INLA.
  # gal_mesh <- galerkin_db(faces, vertices, surface = TRUE)
  # mesh <- list(
  #   n = nrow(gal_mesh$C),
  #   loc = vertices,
  #   graph = list(tv = faces),
  #   idx = list(loc = seq(nrow(gal_mesh$C)))
  # )
  # class(mesh) <- "inla.mesh"
  
  return(mesh)
}
