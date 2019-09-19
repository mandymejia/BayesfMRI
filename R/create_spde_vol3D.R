#' Create SPDE for 3D volumetric data
#'
#' Create SPDE for 3D volumetric data
#'
#' @param loc Locations of data points
#' @param lab Labels of data points
#' @param value Label value to analyze
#' @param max_dist Maximum distance of new vertex locations from the original vertex locations
#'
#' @return SPDE
#' @export
#' @importFrom INLA inla.mesh.create inla.spde2.generic
#' @importFrom geometry delaunayn
#' @importFrom Matrix sparseMatrix colSums Diagonal t solve
#' @importFrom rdist cdist
#'
# require(INLA)
# library(geometry)
# library(rgl)
# # Simulate locations
# locs <- matrix(rnorm(10*3, mean =0, sd =3),10,3)
# lab <- NULL # sample(1:3, 10, replace = TRUE)
# value <- NULL
# max_dist <- 1
#
# # Read in locations from actual data
# locs <- read.csv('/Users/sarahryan/Desktop/loc.100307.csv')
# lab <- read.csv('/Users/sarahryan/Desktop/labels.100307.csv')
# value <- 21
# max_dist <- 1
#
# # Create SPDE
# spde <- make_mesh_vol3D(locs, lab, value = 1)
#
# # For example, sample the field:
# Q = inla.spde2.precision(spde,theta = c(log(0.1),log(1)))
# x = inla.qsample(Q=Q)
create_spde_vol3D <- function(locs, lab, value = NULL, max_dist = NULL){

  # Which part to analyze
  if(!is.null(value)){
    ind <- lab == value
  }else{
    ind <- rep(TRUE, nrow(locs))
  }
  P <- locs[ind,]

  # Create Triangulations for 3D data from original locations. Note that there are big triangles that we don't want to keep.
  # FV_orig <- delaunayn(P)
    # open3d()
    # rgl.viewpoint(60)
    # rgl.light(120,60)
    # tetramesh(FV_orig, P, alpha=0.7)


  # Create grid to get reasonable delaunay triangulations.. Reasonable being triangles with similar size.
  x = unique(P[,1])
  y = unique(P[,2])
  z = unique(P[,3])
  PP <- expand.grid(x, y, z)

  # Create Triangulations for 3D data based on the grid
  FV <- delaunayn(PP)
    # open3d()
    # rgl.viewpoint(60)
    # rgl.light(120,60)
    # tetramesh(FV, PP, alpha=0.7)


  # Find extra edges to remove (i.e. edges that don't connect our original vertices, but are in our triangualations since we created a grid)
  D <- cdist(P,PP)
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
  # Now P_new and FV contain points not in P originally
    # # Visualize the mesh
    # open3d()
    # rgl.viewpoint(60)
    # rgl.light(120,60)
    # tetramesh(FV_new, P_new, alpha=0.7)

  # Create observation matrix A
  D <- cdist(P,P_new)
  I <- apply(D, 1, function(x) {which(x == min(x))})
  A <- diag(1, nrow = nrow(P_new))
  A <- A[I,]

  # ij <- which(A!=0,arr.ind = T)
  # i <- ij[,1]
  # j <- ij[,2]

  # Create FEM matrices
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

    G <- sparseMatrix(i = as.vector(Gi), j = as.vector(Gj), x = as.vector(Gz), dims = c(nV,nV))
    Ce <- sparseMatrix(i = as.vector(Ci), j = as.vector(Cj), x = as.vector(Cz), dims = c(nV,nV))
    C <- diag(colSums(Ce), nrow = nV, ncol = nV)
    return(list(G = G, Ce = Ce, C = C))
  }

  gal <- galerkin_db(FV = FV_new, P = P_new)
  G <- gal$G
  Ce <- gal$Ce
  C <- gal$C

  # Part 2
  tG <- t(G)
  M0 = C
  M1 = G + tG
  M2 = tG %*% solve(C,G)


  # Create the spde object. Note that the priors theta.mu and theta.Q have to be set reasonably here!!!
  spde = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2,
                            B0 = matrix(c(0,1,0),1,3),
                            B1 = matrix(c(0,0,1),1,3),
                            B2 = 1,
                            theta.mu = rep(0,2),
                            theta.Q = Diagonal(2,c(1e-6,1e-6)),
                            transform = "identity")

  return(spde)
}




