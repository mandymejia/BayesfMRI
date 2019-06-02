#' Make Mesh
#'
#' Make triangular mesh from faces and vertices. If mask is given, then mesh is masked.
#'
#' @param faces Matrix of faces
#' @param verts Matrix of vertices
#' @param mask Mask to be applied to mesh. If none given, full mesh is used.
#'
#' @return Triangular mesh from matrices and vertices
#' @export
#' @importFrom INLA inla.mesh.create
make_mesh <- function(faces, verts, mask = NULL){

  # Number of vertices
  V <- nrow(verts)

  # Check the dimension of mask
  if(is.null(mask)){
    mask <- rep(1, V)
  }else{
    mask <- is.numeric(mask)
  }
  if(length(mask) != V | !is.vector(mask)){
    stop("Mask should be a vector of length V")
  }

  # Check 0s and 1s
  values <- sort(unique(mask))
  if(values[1] != 0 | values[2] !=1 | length(values) != 2){
    stop("Mask should be composed of only 0s and 1s")
  }

  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  # Apply mask to vertices
  inmask <- which(mask)
  verts <- verts[mask,]

  # Identify and remove any triangles where at least one vertex is not included in the motor mask
  faces <- faces[(faces[,1] %in% inmask) & (faces[,2] %in% inmask) & (faces[,3] %in% inmask),]

  # Re-number faces
  faces_new <- faces*0
  for(ii in 1:nrow(faces_sh)){
    faces_new[ii,1] <- which(inmask == faces[ii,1])
    faces_new[ii,2] <- which(inmask == faces[ii,2])
    faces_new[ii,3] <- which(inmask == faces[ii,3])
  }

  # Construct mesh
  mesh <- inla.mesh.create(loc = as.matrix(verts), tv = as.matrix(faces_new))
  return(mesh)
}


