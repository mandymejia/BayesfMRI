#' Make Mesh
#'
#' Make INLA triangular mesh from faces and vertices
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @param use_INLA (logical) Use the INLA package to make the mesh? Default:
#'  \code{FALSE}. Otherwise, mesh construction is based on an internal function,
#'  \code{galerkin_db}.
#'
#' @return INLA triangular mesh
#'
#' @export
make_mesh <- function(vertices, faces, use_INLA = FALSE){

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("`make_mesh` requires the `INLA` package. Please install it.", call. = FALSE)
  }


  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  # Construct mesh
  if(use_INLA) {
    mesh <- INLA::inla.mesh.create(loc = as.matrix(vertices), tv = as.matrix(faces))
  } else {
    gal_mesh <- galerkin_db(faces, vertices, surface = TRUE)
    mesh <- list(
      n = nrow(gal_mesh$C),
      loc = vertices,
      graph = list(tv = faces),
      idx = list(loc = seq(nrow(gal_mesh$C)))
    )
    class(mesh) <- "inla.mesh"
  }
  return(mesh)
}

#' Remove part of a mesh without using INLA functions
#'
#' @param mask a 0-1 vector
#' @param mesh a mesh resulting from a call to \code{make_mesh} with
#'   \code{use_INLA = FALSE}
#'
#' @return a mesh object with fewer vertices than the original input mesh
#' @keywords internal
submesh <- function(mask, mesh) {
  t.count <- Matrix::rowSums(matrix((mask >= 0.5)[mesh$graph$tv], nrow(mesh$graph$tv),3))
  tri <- which(t.count == 3)
  tv <- mesh$graph$tv[tri, , drop = FALSE]
  v <- sort(unique(as.vector(tv)))
  idx <- rep(as.integer(NA), nrow(mesh$loc))
  idx[v] <- seq_len(length(v))
  tv <- matrix(idx[tv], nrow(tv), 3)
  loc <- mesh$loc[v, , drop = FALSE]
  # mesh <- INLA::inla.mesh.create(loc = loc, tv = tv, refine = FALSE)
  mesh <- make_mesh(vertices = loc,faces = tv,use_INLA = FALSE)
  idx <- rep(as.integer(NA), length(idx))
  idx[v] <- mesh$idx$loc
  mesh$idx$loc <- idx
  class(mesh) <- "inla.mesh"
  return(mesh)
}
