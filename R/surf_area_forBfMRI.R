# "forBfMRI" because `ciftiTools` has a script already,
#  but here we include the INLA-based routine too and
#  prioritize that over our routine in `ciftiTools`.

#' @rdname surf_area
#' @importFrom INLA inla.mesh.create inla.fmesher.smorg
#' @keywords internal
surf_area_INLA <- function(surf, by=c("vertex", "face")) {  
  mesh <- inla.mesh.create(
    loc = as.matrix(surf$vertices), 
    tv = as.matrix(surf$faces)
  )
  diag(inla.fmesher.smorg(
    mesh$loc, mesh$graph$tv, fem = 0, output = list("c0")
  )$c0)
}

#' Surface area calculation
#' 
#' Calculate surface area of a \code{"surf"} object by vertex or face.
#' 
#' @param mesh The \code{"BfMRI.mesh"} object.
#' @param by \code{"vertex"} or \code{"face"}. For \code{"vertex"}, the result
#'  is the area associated with each vertex: the sum the area of each triangular
#'  face it is a part of, divided by three. For \code{"face"}, the result is
#'  the surface area of each face. 
#' 
#' @return Vector of surface areas \code{by} vertex or face, in the same order
#'  as the corresponding component of \code{surf}. The units are the square of
#'  the units of \code{surf$vertices}. 
#' @importFrom ciftiTools surf_area
#' @keywords internal
surf_area_forBfMRImesh <- function(mesh, by=c("vertex", "face")) {
  surf <- mesh_to_surf(mesh)
  if (requireNamespace("INLA", quietly = TRUE)) {
    surf_area_INLA(surf, by)
  } else {
    ciftiTools::surf_area(surf, by)
  }
}