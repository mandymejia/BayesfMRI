#' Apply Mask to Vertices and Faces
#'
#' Apply a binary mask to a set of vertices and faces.  Vertices not in the mask are removed,
#' and faces (triangles) with any vertices outside of the mask are removed.  Finally,
#' vertex numbering in the masked faces matrix is corrected.
#'
#' @param vertices Vx3 matrix of vertices
#' @param faces matrix of faces
#' @param mask 0/1 vector of length V indicating which vertices are to be retained (1) and which are to be excluded (0)
#'
#' @return List containing masked vertices and faces matrices
#' @export
mask_vertices_faces <- function(vertices, faces, mask){

  # Number of vertices
  V <- nrow(vertices)

  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  mask <- as.numeric(mask)
  if(length(mask) != V | !is.vector(mask)){
    stop("Mask should be a vector of length V")
  }
  # Check only 0s and 1s
  values <- sort(unique(mask))
  if(! (min(values %in% 0:1)) ) stop("Mask should be composed of only 0s and 1s")

  inmask <- which(mask==1)

  # Apply mask to vertices
  vertices_new <- vertices[inmask,]

  ### Apply mask to faces (triangles)

  # Identify triangles where any vertex is outside of the mask
  faces <- faces[(faces[,1] %in% inmask) & (faces[,2] %in% inmask) & (faces[,3] %in% inmask),]

  # Re-number faces
  faces_new <- faces*0
  for(ii in 1:nrow(faces)){
    faces_new[ii,1] <- which(inmask == faces[ii,1])
    faces_new[ii,2] <- which(inmask == faces[ii,2])
    faces_new[ii,3] <- which(inmask == faces[ii,3])
  }

  # Return updated vertices and faces
  result <- list(vertices=vertices_new, faces=faces_new)
  return(result)
}


