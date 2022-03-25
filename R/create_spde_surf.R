#' Create SPDE for 2D surface data
#'
#' @param mesh A list with elements \code{n} (the number of vertices),
#'   \code{loc} (the n x 3 matrix of vertices), \code{graph$tv} (the matrix of
#'   faces with 3 columns), and \code{idx$loc} (a vector of the indices ,
#'   usually 1:n)
#'
#' @return SPDE object representing triangular mesh structure on data locations
#' @export
create_spde_surf <- function(mesh) {
  gal <- galerkin_db(mesh$graph$tv, mesh$loc,surface = TRUE)

  G <- gal$G
  M0 <- gal$C

  tG <- Matrix::t(G)
  M1 = G + tG
  M2 = tG %*% Matrix::solve(M0,G)

  spde = list(M0 = M0, M1 = M1, M2 = M2,n.spde = nrow(M0))

  out <- list(spde = spde,
              vertices = mesh$loc,
              faces = mesh$graph$tv,
              idx = seq(nrow(M0)),
              Amat = Matrix::Diagonal(n = nrow(M0),x = 1))
  return(out)
}
