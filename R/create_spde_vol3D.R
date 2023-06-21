# #' Create SPDE for 3D volumetric data
# #'
# #' @inheritSection INLA_Description INLA Requirement
# #'
# #' @param locs Locations of data points (Vx3 matrix)
# #' @param labs Region labels of data points (vector of length V). If NULL, treat observations as a single region.
# #' @param lab_set Only used if labs is not NULL. Vector of region labels for which to construct spde object. If NULL, use all region labels.
# #'
# #' @return SPDE object representing triangular mesh structure on data locations
# #'
# #' @importFrom Matrix sparseMatrix colSums Diagonal t solve
# #'
# #' @export
# create_spde_vol3D <- function(locs, labs, lab_set = NULL){

#   # Check to see that the `rdist` package is installed
#   if (!requireNamespace("rdist", quietly = TRUE)) {
#     stop("`create_spde_vol3D` requires the `rdist` package. Please install it.", call. = FALSE)
#   }

#   # Check to see that the `geometry` package is installed
#   if (!requireNamespace("geometry", quietly = TRUE)) {
#     stop("`create_spde_vol3D` requires the `geometry` package. Please install it.", call. = FALSE)
#   }

#   check_INLA(FALSE)

#   if(is.null(labs) & !is.null(lab_set)) stop('If labs is NULL, then lab_set must not be specified.')

#   # If no labels provided, construct mesh over all regions (treat as a single region)
#   if(is.null(labs)) labs <- rep(1, nrow(locs))

#   # If no set of labels specified, use all labels
#   if(is.null(lab_set)) lab_set <- unique(labs)

#   # For each of the labels specified, create a triangularization

#   G_all <- C_all <- vector('list', length=length(lab_set))
#   P_new_all <- FV_new_all <- vector('list', length=length(lab_set))
#   I_all <- A_all <- vector('list', length=length(lab_set))
#   for(value in lab_set){

#     ii <- which(lab_set==value)
#     ind <- (labs == value)
#     P <- locs[ind,] #(x,y,z) coordinates of selected locations

#     # In Triangulations for 3D data from original locations, note that there are big triangles that we don't want to keep.
#     # FV_orig <- geometry::delaunayn(P)
#     # open3d()
#     # rgl.viewpoint(60)
#     # rgl.light(120,60)
#     # tetramesh(FV_orig, P, alpha=0.7)

#     # Determine grid spacing (assumes isotropic voxels)
#     max_dist = get_spacing(P[,1])

#     # Create grid to get reasonable delaunay triangulations. Reasonable being triangles with similar size.
#     x = unique(P[,1])
#     y = unique(P[,2])
#     z = unique(P[,3])
#     #add vertices around the extreme of the boundary
#     x = c(x, min(x)-max_dist, max(x)+max_dist)
#     y = c(y, min(y)-max_dist, max(y)+max_dist)
#     z = c(z, min(z)-max_dist, max(z)+max_dist)

#     # # Subsample to a specified proportion
#     # if(!is.null(subsample)){
#     #   if(subsample < 0 | subsample > 1) stop('subsample must be a number between 0 and 1')
#     #   length.out.x <- length(x)*subsample
#     #   length.out.y <- length(y)*subsample
#     #   length.out.z <- length(z)*subsample
#     #   x = seq(min(x)-max_dist, max(x)+max_dist, length.out=length.out.x)
#     #   y = seq(min(y)-max_dist, max(y)+max_dist, length.out=length.out.y)
#     #   z = seq(min(z)-max_dist, max(z)+max_dist, length.out=length.out.z)
#     # }

#     PP <- expand.grid(x, y, z) #full lattice within a cube surrounding the data locations

#     # Create Triangulations for 3D data based on the grid
#     FV <- geometry::delaunayn(PP)
#     # open3d()
#     # rgl.viewpoint(60)
#     # rgl.light(120,60)
#     # tetramesh(FV, PP, alpha=0.7)

#     # Remove locations that are far from original data locations
#     D <- rdist::cdist(P,PP)
#     md <- apply(D, 2, min)
#     ind_keep <- md < max_dist
#     indices <- which(ind_keep == 1)

#     # Remove vertices from FV that are not associated with an original vertex
#     FV_new <- NULL
#     for (i in 1:nrow(FV)) {
#       if (sum(FV[i,] %in% indices) > 0) {
#         FV_new = rbind(FV_new, FV[i,])
#       }
#     }
#     fu <- unique(as.vector(FV_new))
#     v <- rep(0, max(fu))
#     v[fu] <- 1:length(fu)
#     for (i in 1:nrow(FV_new)){
#       FV_new[i,] <- v[FV_new[i,]]
#     }
#     P_new <- PP[fu,]
#     FV_new_all[[ii]] <- FV_new
#     P_new_all[[ii]] <- P_new
#     # Now P_new and FV contain points not in P originally
#       # # Visualize the mesh
#       # open3d()
#       # rgl.viewpoint(60)
#       # rgl.light(120,60)
#       # tetramesh(FV_new, P_new, alpha=0.7)

#     # Create observation matrix A (this assumes that all original data locations appear in mesh)
#     D <- rdist::cdist(P,P_new)
#     I_v <- apply(D, 1, function(x) {which(x == min(x))})
#     A_v <- diag(1, nrow = nrow(P_new))
#     A_all[[ii]] <- A_v[I_v,]
#     I_all[[ii]] <- I_v

#     # ij <- which(A!=0,arr.ind = T)
#     # i <- ij[,1]
#     # j <- ij[,2]

#     #construct G and C matrices that appear in SPDE precision
#     gal <- galerkin_db(FV = FV_new, P = P_new)
#     G_all[[ii]] <- gal$G
#     C_all[[ii]] <- gal$C
#   }

#   G <- bdiag(G_all)
#   C <- bdiag(C_all)

#   # Part 2
#   tG <- t(G)
#   M0 <- C
#   M1 <- G + tG
#   M2 <- tG %*% solve(C,G)

#   # Create the spde object. Note that the priors theta.mu and theta.Q have to be set reasonably here!!!
#   # spde <- INLA::inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2,
#   #                           B0 = matrix(c(0,1,0),1,3),
#   #                           B1 = matrix(c(0,0,1),1,3),
#   #                           B2 = 1,
#   #                           theta.mu = rep(0,2),
#   #                           theta.Q = Diagonal(2,c(1e-6,1e-6)),
#   #                           transform = "identity")

#   spde = list(M0 = M0, M1 = M1, M2 = M2,n.spde = nrow(M0))

#   out <- list(spde = spde,
#               vertices = P_new_all,
#               faces = FV_new_all,
#               idx = I_all,
#               Amat = bdiag(A_all))
#   class(out) <- 'BayesfMRI.spde'

#   return(out)
# }

# #' Determine grid spacing
# #'
# #' @param locations Vector of vertex locations in one dimension
# #'
# #' @return Value of minimum spacing between two locations
# #' @export
# get_spacing <- function(locations){
#   #locations is a vector of locations along one dimension
#   x = sort(unique(locations))
#   dx = diff(x)
#   return(min(dx))
# }


#' Create FEM matrices
#'
#' @param FV Matrix of faces in triangularization
#' @param P Matrix of vertex locations in triangularization
#' @param surface (logical) Will this create the SPDE matrices for a surface
#'   or not?
#'
#' @return A list of matrices C and G appearing in sparse SPDE precision
#'
#' @keywords internal
galerkin_db <- function(FV, P, surface=FALSE){
  d <- ncol(FV)-1
  if(surface){
    if(ncol(P) != (d + 1)){P <- t(P)}
    if(ncol(P) != (d + 1)){stop("Wrong dimension of P")}
  } else {
    if(ncol(P) != d){P <- t(P)}
    if(ncol(P) != d){stop("Wrong dimension of P")}
  }

  nV <- nrow(P)
  nF <- nrow(FV)
  Gi <- matrix(0, nrow = nF*(d+1), ncol = d+1)
  Gj = Gi; Gz = Gi; Ci = Gi; Cj = Gi; Cz = Gi;

  for( f in 1:nF){
    dd <- (d+1)*(f-1)+(1:(d+1))
    Gi[dd,] <- Ci[dd,] <- FV[f,] %*% t(rep(1,d+1))
    Gj[dd,] <- Cj[dd,] <- t(Gi[dd,])
    if(surface){
      r = t(P[FV[f,c(3,1,2)],]-P[FV[f,c(2,3,1)],])
      r1 = r[,1,drop=FALSE]
      r2 = r[,2,drop=FALSE]
      f_area <- as.double(sqrt((t(r1)%*%r1)*(t(r2)%*%r2)-(t(r1)%*%r2)^2)/2)
      Gz[dd,] = (t(r)%*%r)/(4*f_area)
      Cz[dd,] = (f_area/12)*(matrix(1,3,3)+diag(3))
    } else {
      m1 <- rbind(rep(1, d+1), t(P[FV[f,],]))
      m2 <- rbind(rep(0, d), diag(1, d))
      m <- solve(m1, m2)
      ddet <- abs(det(m1))
      Gz[dd,] <- ddet * (m %*% t(m)) / prod(1:d)
      Cz[dd,] <- ddet * (rep(1,d+1)+diag(d+1)) / prod(1:(d+2))
    }
  }

  G <- Matrix::sparseMatrix(i = as.vector(Gi), j = as.vector(Gj), x = as.vector(Gz), dims = c(nV,nV))
  Ce <- Matrix::sparseMatrix(i = as.vector(Ci), j = as.vector(Cj), x = as.vector(Cz), dims = c(nV,nV))
  # C <- Matrix::diag(Matrix::colSums(Ce), nrow = nV, ncol = nV)
  C <- Matrix::Diagonal(n = nV, x = Matrix::colSums(Ce))
  return(list(G = G, Ce = Ce, C = C))
}






