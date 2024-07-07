#' Make A matrix
#'
#' Make A matrix
#'
#' @param spatial See \code{BayesGLM}
#' @return The A matrix
#' @keywords internal
make_A_mat <- function(spatial){
  nV <- get_nV(spatial)

  ### `A`. -----
	# [TO DO]
  # A_sparse <- make_A_mat_rs(
  #   spatial$surf,
  #   ciftiTools::mask_surf(spatial$surf, spatial$mask)
  # )

  A_sparse <- Matrix::Diagonal(nV$model)[spatial$maskMdat,]
}

#' Make A matrix with resampling framework
#'
#' Make the A matrix for downsampling surface mesh data to a lower resolution.
#'
#' @param surf The full-resolution \code{"surf"} object.
#' @param surf_rs The downsampled \code{"surf"} object.
#' @return The A matrix
#' @keywords internal
make_A_mat_rs <- function(surf, surf_rs){

  stop("[TO DO] refine")

  stopifnot(ciftiTools::is.surf(surf))
  stopifnot(ciftiTools::is.surf(surf_rs))

  #construct the mesh
  v1 <- surf_rs$vertices
  f1 <- surf_rs$faces
  v1.sums <- sqrt(rowSums(v1^2))
  v1.norm <- v1
  for (ii in seq(length(v1.sums))) {
    v1.norm[ii,] <- v1[ii,]/v1.sums[1]
  }
  v1.norm <- v1.norm*v1.sums[1]
  mesh <- INLA::inla.mesh.create(loc = v1.norm, tv = f1) #downsampled mesh

  #determine the data location coordinates
  v0 <- surf$vertices
  v0.sums <- sqrt(rowSums(v0^2))
  loc0 <- v0
  for (ii in seq(length(v0.sums))) {
    loc0[ii,] <- loc0[ii,]/v0.sums[1]
  }
  loc0 <- loc0*v0.sums[1]

  #compute the A matrix mapping between data locations and mesh locations
  INLA::inla.spde.make.A(mesh = mesh, loc = loc0)
}
