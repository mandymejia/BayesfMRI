#' Construct a triangular mesh from a 3D volumetric mask
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param mask An array of 0s and 1s representing a volumetric mask
#' @param res The spatial resolution in each direction, in mm. For example, c(2,2,2) indicates 2mm isotropic voxels.
#' @param buffer Size of extra voxels layers around the bounding box, in terms of voxels
#' @param neighbor_order What order of neighbors of data locations to keep? (1 = first-order neighbors, 2 = first- and second-order neighbors, etc.)
#'
#' @return An inla.spde2 object.
#' 
#' @export
#'
vol2spde <- function(mask, res = c(2,2,2), buffer=c(1,1,3,4,4), neighbor_order = 1){

  #[TO DO] check that mask is an array and contains 0s and 1s. If logical, make numeric.
  #[TO DO] check that buffer is increasing positive integers
  #[TO DO] check that the buffer argument is large enough to capture all the dependence layers, plus one.  Each dependence layer is 2 neighbor layers.

  #bounding box limits for the mask
  lim_x <- range(which(apply(mask, 1, sum) > 0)) #dim 1
  lim_y <- range(which(apply(mask, 2, sum) > 0)) #dim 2
  lim_z <- range(which(apply(mask, 3, sum) > 0)) #dim 3

  #indices (x0, y0 and z0 guaranteed to lay within original mask)
  x0 <- x <- (lim_x[1]:lim_x[2])
  y0 <- y <- (lim_y[1]:lim_y[2])
  z0 <- z <- (lim_z[1]:lim_z[2])

  #boundary layer(s)
  nB <- length(buffer)
  for(b in 1:nB){
    size_b <- sum(buffer[1:b])
    x <- c(min(x0)-size_b, x, max(x0) + size_b)
    y <- c(min(y0)-size_b, y, max(y0) + size_b)
    z <- c(min(z0)-size_b, z, max(z0) + size_b)
  }

  #define mask within expanded box
  mask_box0 <- mask[x0, y0, z0] #mask within the original bounding box
  mask_box <- array(0, dim=c(length(x), length(y), length(z)))
  mask_box[x0 + nB - min(x0) + 1, y0 + nB - min(y0) + 1, z0 + nB - min(z0) + 1] <- mask_box0 #mask within expanded box

  #convert indices to Euclidean coordinates based on voxel size
  x <- x*res[1]
  y <- y*res[2]
  z <- z*res[3]

  #x <- seq(from=0,to=1,length.out = length(x0))
  #y <- seq(from=0,to=1,length.out = length(y0))
  #z <- seq(from=0,to=1,length.out = length(z0))

  #dimensions
  nx <- length(x)
  ny <- length(y)
  nz <- length(z)

  #1D meshes
  mesh.x <- INLA::inla.mesh.1d(loc = x)
  mesh.y <- INLA::inla.mesh.1d(loc = y)
  mesh.z <- INLA::inla.mesh.1d(loc = z)

  fem.x <- INLA::inla.mesh.fem(mesh.x)
  fem.y <- INLA::inla.mesh.fem(mesh.y)
  fem.z <- INLA::inla.mesh.fem(mesh.z)

  #construct C and G for 3D mesh
  Cx <- fem.x$c0
  Gx <- fem.x$g1
  Cy <- fem.y$c0
  Gy <- fem.y$g1
  Cz <- fem.z$c0
  Gz <- fem.z$g1

  C3d <- kronecker(kronecker(Cz,Cy),Cx)
  G3d <- kronecker(kronecker(Cz,Cy),Gx) + kronecker(kronecker(Cz,Gy),Cx) + kronecker(kronecker(Gz,Cy),Cx)

  # spde <- INLA::inla.spde2.generic(M0 = C3d,
  #                            M1 = G3d,
  #                            M2 = G3d%*%solve(C3d, G3d),
  #                            theta.mu = c(Elog.kappa, Elog.tau),
  #                            theta.Q = c(Qlog.kappa, Qlog.tau),
  #                            B0 = matrix(c(0, 1, 0), 1, 3),
  #                            B1 = 2*matrix(c(0, 0, 1), 1, 3),
  #                            B2 = 1)
  #
  # #visualize marginal variances inside and outside of mask
  # Q_full <- INLA::inla.spde.precision(spde, theta = c(1,1)) #example Q with low spatial range (large kappa)
  # tmp_full <- diag(solve(Q_full)) #marginal variances
  # tmp_mask <- mask_box; tmp_mask[mask_box >= 0] <- tmp_full
  # tmp_mask[10,1,1] <- 0; tmp_mask[10,1,2] <- 0.0004 #to control the scale
  # image.plot(1000*tmp_mask[10,,], col=viridis(100))
  # mask_box_NAs <- mask_box; mask_box_NAs[mask_box==0] <- NA
  # tmp_mask_NAs <- tmp_mask*mask_box_NAs #mask out excluded locations
  # image.plot(1000*tmp_mask_NAs[10,,], col=viridis(100))

  #M2 has the same sparsity structure as Q
  M2 <- G3d%*%solve(C3d, G3d) #encodes first-order neighbors of data locations
  if(neighbor_order > 1){
    for(o in 2:neighbor_order){
      M2 <- M2 %*% M2 #encode the next order of neighbors
    }
  }

  #identify locations to remove (no dependence on any in-mask locations)
  idx <- which(mask_box==1) #indices of in-mask locations
  idx2 <- apply(M2[idx,], 1, function(x) which(x != 0)) #locations with any dependence with in-mask locations
  idx2 <- unique(as.vector(idx2)) #1252 locations for 257 original ROI locations
  mask_box2 <- mask_box; mask_box2[idx2] <- mask_box[idx2] + 1 #for visualization

  #recreate the SPDE after removing locations
  C3d <- C3d[idx2,idx2]
  G3d <- G3d[idx2,idx2]

  #[TO DO] Construct SPDE outside of this function, block diagonalizing over regions

  # spde <- INLA::inla.spde2.generic(M0 = C3d,
  #                            M1 = G3d,
  #                            M2 = G3d%*%solve(C3d, G3d),
  #                            theta.mu = c(Elog.kappa, Elog.tau),
  #                            theta.Q = c(Qlog.kappa, Qlog.tau),
  #                            B0 = matrix(c(0, 1, 0), 1, 3),
  #                            B1 = 2*matrix(c(0, 0, 1), 1, 3),
  #                            B2 = 1)
  #
  # #visualize marginal variances inside and outside of mask
  # Q <- INLA::inla.spde.precision(spde, theta = c(1,1)) #example Q with low spatial range (large kappa)
  # tmp <- diag(solve(Q)) #marginal variances
  # tmp_full <- rep(0, length(mask_box)); tmp_full[idx2] <- tmp #put back in excluded locations
  # tmp_mask <- mask_box; tmp_mask[mask_box >= 0] <- tmp_full; tmp_mask[mask_box2 == 0] <- NA
  # tmp_mask[10,1,1] <- 0; tmp_mask[10,1,2] <- 0.0004 #to control the scale
  # image.plot(1000*tmp_mask[10,,], col=viridis(100))
  # mask_box_NAs <- mask_box; mask_box_NAs[mask_box==0] <- NA
  # tmp_mask_NAs <- tmp_mask*mask_box_NAs #mask out excluded locations
  # image.plot(1000*tmp_mask_NAs[10,,], col=viridis(100))
  # par(mfrow=c(3,3))
  # summary(1000*tmp_mask_NAs, na.rm=TRUE)

  params <- list(res = res,
                 buffer = buffer)

  list(mats = list(C = C3d, G = G3d),
       idx2 = idx2, #idx of included locations in expanded box
       xyz0 = list(x0=x0, y0=y0, z0=z0), #bounding box
       xyz = list(x=x, y=y, z=z), #expanded bounding box
       mask_orig = mask, #original mask within original volume
       mask_box = mask_box, #original mask within expanded bounding box
       mask_box2 = mask_box2, #original mask + boundary layers
       params = params)
}
