#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?is.session} for more details.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh.
#' @param mask A vector of 0s and 1s of length V, where locations with value 0 will be excluded from analysis.
#'
#' @return A list containing...
#' @export
#'
#' @examples
BayesGLMfMRI <- function(data, vertices, faces, mask=NULL){

  INLA:::inla.dynload.workaround() #avoid error on creating mesh

  mesh <- make_mesh(vertices, faces, mask)
  spde <- inla.mesh.create(mesh)
  areas <- compute_vertex_areas(mesh)

}
