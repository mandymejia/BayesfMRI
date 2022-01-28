#' Create SPDE for 3D volumetric data
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param locs Locations of data points (Vx3 matrix)
#' @param labs Region labels of data points (vector of length V). If NULL, treat observations as a single region.
#' @param lab_set Only used if labs is not NULL. Vector of region labels for which to construct spde object. If NULL, use all region labels.
#' @return SPDE object representing triangular mesh structure on data locations
#'
#' @importFrom geometry delaunayn
#' @importFrom INLA inla.spde2.generic
#' @importFrom Matrix diag bdiag
#' @importFrom spam spam_rdist apply.spam diag bdiag.spam as.spam crossprod.spam solve.spam
#'
#' @export
create_spde_vol3D <- function(locs, labs, lab_set = NULL){
  if(is.null(labs) & !is.null(lab_set)) stop('If labs is NULL, then lab_set must not be specified.')

  # If no labels provided, construct mesh over all regions (treat as a single region)
  if(is.null(labs)) labs <- rep(1, nrow(locs))

  # If no set of labels specified, use all labels
  if(is.null(lab_set)) lab_set <- unique(labs)

  # For each of the labels specified, create a triangularization

  G_all <- C_all <- vector('list', length=length(lab_set))
  P_new_all <- FV_new_all <- vector('list', length=length(lab_set))
  I_all <- A_all <- vector('list', length=length(lab_set))
  for(value in lab_set){

    ii <- which(lab_set==value)
    ind <- (as.numeric(labs) == value)
    P <- locs[ind,] #(x,y,z) coordinates of selected locations

    # Determine grid spacing (assumes isotropic voxels)
    sorted_locs = sort(unique(P[,1]))
    diff_locs = diff(sorted_locs)
    max_dist <- min(diff_locs)

    # Create grid to get reasonable delaunay triangulations. Reasonable being triangles with similar size.
    x = unique(P[,1])
    y = unique(P[,2])
    z = unique(P[,3])
    #add vertices around the extreme of the boundary
    x = c(x, min(x)-max_dist, max(x)+max_dist)
    y = c(y, min(y)-max_dist, max(y)+max_dist)
    z = c(z, min(z)-max_dist, max(z)+max_dist)

    PP <- expand.grid(x, y, z) #full lattice within a cube surrounding the data locations

    # Create Triangulations for 3D data based on the grid
    FV <- geometry::delaunayn(PP)

    # Remove locations that are far from original data locations
    D <- rdist::cdist(P,PP)
    md <- apply(D, 2, min)
    ind_keep <- md < max_dist
    indices <- which(ind_keep == 1)

    # Remove vertices from FV that are not associated with an original vertex
    FV_new <- NULL
    for (i in 1:nrow(FV)) {
      if (sum(FV[i,] %in% indices) > 0) {
        FV_new = rbind(FV_new, FV[i,])
      }
    }
    fu <- unique(as.vector(FV_new))
    v <- rep(0, max(fu))
    v[fu] <- 1:length(fu)
    for (i in 1:nrow(FV_new)){
      FV_new[i,] <- v[FV_new[i,]]
    }
    P_new <- PP[fu,]
    FV_new_all[[ii]] <- FV_new
    P_new_all[[ii]] <- P_new

    # Create observation matrix A (this assumes that all original data locations appear in mesh)
    D <- rdist::cdist(P,P_new)
    I_v <- apply(D, 1, function(x) {which(x == min(x))})
    A_v <- Matrix::diag(1, nrow = nrow(P_new))
    A_all[[ii]] <- A_v[I_v,]
    I_all[[ii]] <- I_v

    #construct G and C matrices that appear in SPDE precision
    gal <- galerkin_db(FV = FV_new, P = P_new)
    G_all[[ii]] <- gal$G
    C_all[[ii]] <- gal$C
  }

  G <- Matrix::bdiag(G_all)
  C <- Matrix::bdiag(C_all)

  # Part 2
  tG <- Matrix::t(G)
  M0 = C
  M1 = G + tG
  M2 = tG %*% Matrix::solve(C,G)

  # Create the spde object. Note that the priors theta.mu and theta.Q have to be set reasonably here!!!
  # spde = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2,
  #                           B0 = matrix(c(0,1,0),1,3),
  #                           B1 = matrix(c(0,0,1),1,3),
  #                           B2 = 1,
  #                           theta.mu = rep(0,2),
  #                           theta.Q = Matrix::Diagonal(2,c(1e-6,1e-6)),
  #                           transform = "identity")

  spde <- list(
    M0 = M0,
    M1 = M1,
    M2 = M2,
    n.spde = nrow(M0)
  )

  out <- list(spde = spde,
              vertices = P_new_all,
              faces = FV_new_all,
              idx = I_all,
              Amat = Matrix::bdiag(A_all))
  class(out) <- 'BayesfMRI.spde'

  return(out)
}

#' Create FEM matrices
#'
#' @param FV Matrix of faces in triangularization
#' @param P Matrix of vertex locations in triangularization
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom spam as.spam.dgCMatrix diag.spam colSums
#' @return A list of matrices C and G appearing in sparse SPDE precision
#'
#' @export
galerkin_db <- function(FV, P){
  d <- ncol(FV)-1
  if(ncol(P) != d){P <- t(P)}
  nV <- nrow(P)
  nF <- nrow(FV)
  Gi <- matrix(0, nrow = nF*(d+1), ncol = d+1)
  Gj = Gi; Gz = Gi; Ci = Gi; Cj = Gi; Cz = Gi;

  for( f in 1:nF){
    m1 <- rbind(rep(1, d+1), t(P[FV[f,],]))
    m2 <- rbind(rep(0, d), diag(1, d))
    m <- solve(m1, m2)
    ddet <- abs(det(m1))

    dd <- (d+1)*(f-1)+(1:(d+1))
    Gi[dd,] <- FV[f,] %*% t(rep(1,d+1))
    Gj[dd,] <- t(Gi[dd,])
    Gz[dd,] <- ddet * (m %*% t(m)) / prod(1:d)
    Ci[dd,] <- Gi[dd,]
    Cj[dd,] <- Gj[dd,]
    Cz[dd,] <- ddet * (rep(1,d+1)+diag(d+1)) / prod(1:(d+2))
  }

  G <- Matrix::sparseMatrix(i = as.vector(Gi), j = as.vector(Gj), x = as.vector(Gz), dims = c(nV,nV))
  Ce <- Matrix::sparseMatrix(i = as.vector(Ci), j = as.vector(Cj), x = as.vector(Cz), dims = c(nV,nV))
  C <- Matrix::diag(Matrix::colSums(Ce), nrow = nV, ncol = nV)
  return(list(G = G, Ce = Ce, C = C))
}
