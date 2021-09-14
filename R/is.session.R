#' Is session?
#'
#' Verify that the \code{"session"} object is valid.
#'
#' A valid "session" object is a list with the following named fields
#' - BOLD: \eqn{T x V} matrix of BOLD responses, rows are time points, columns are voxels
#' - design: \eqn{T x K} matrix containing the K task regressors
#' - nuisance (optional): \eqn{T x L} matrix containing the L nuisance regressors
#'
#' @param sess A list representing a task fMRI session (see Details)
#'
#' @return True if all checks pass, or an (hopefully useful) error message
#'
#' @export
is.session <- function(sess){

  ## check number of fields
  num_fields <- length(sess)

  ## BOLD, design, nuisance
  if(num_fields == 3){

    ## check identities of fields
    fields = c('BOLD','design','nuisance')
    if(! all.equal(names(sess),fields)){stop(
      paste0('You are missing the following fields',setdiff(names(sess),fields)))}

    ## check each field's type
    if(! (is.numeric(sess$BOLD))){stop('I expected BOLD to be numeric, but it is not')}
    if(! (is.matrix(sess$BOLD))){stop('I expected BOLD to be a matrix, but it is not')}

    if(! (is.matrix(sess$design))){stop('I expected design to be a matrix, but it is not')}
    if(! (is.matrix(sess$design))){stop('I expected design to be a matrix, but it is not')}

    if(! (is.matrix(sess$nuisance))){stop('I expected nuisance to be a matrix, but it is not')}
    if(! (is.matrix(sess$nuisance))){stop('I expected nuisance to be a matrix, but it is not')}

    ## check the dimensions of each field: T
    if(nrow(sess$BOLD) != nrow(sess$design)){stop("BOLD and design don't have the same number of rows (time points)")}
    if(nrow(sess$BOLD) != nrow(sess$nuisance)){stop("BOLD and nuisance don't have the same number of rows (time points)")}
    if(nrow(sess$design) != nrow(sess$nuisance)){stop("design and nuisance don't have the same number of rows (time points)")}

  ## BOLD, design
  } else if (num_fields==2){

    ## check identities of fields
    fields = c('BOLD','design')
    if(! all.equal(names(sess),fields)){
      stop(paste0('You are missing the following fields',setdiff(names(sess),fields)))
    }

    ## check each field's type
    if(! (is.numeric(sess$BOLD))){stop('I expected BOLD to be numeric, but it is not')}
    if(! (is.matrix(sess$BOLD))){stop('I expected BOLD to be a matrix, but it is not')}

    if(!('matrix' %in% class(sess$design) | "dgCMatrix" %in% class(sess$design))){stop('I expected design to be a matrix, but it is not')}
    is_pw <- !is.matrix(sess$design) #if prewhitening has been done, design is a large sparse matrix (class dgCMatrix)

    if(!is_pw){
      if(nrow(sess$BOLD) != nrow(sess$design)){stop("BOLD and design don't have the same number of rows (time points)")}
    } else {
      if(class(sess$design) != "dgCMatrix"){stop(paste0('I expected the class of design to be dgCMatrix (for prewhitened data) or matrix (for non-prewhitened data), but it is of class ', class(sess$design)))}
      is_missing <- is.na(sess$BOLD[1,])
      nvox <- sum(!is_missing)
      if(nrow(sess$BOLD) != nrow(sess$design)/nvox){stop("If prewhitening has been performed, the number of rows in the design matrix should be T*V, where T=nrow(BOLD) and V is the number of columns of BOLD that are non-NA.")}
    }

  } else {
    stop('I expected the session to have 2 or 3 fields, but it does not.')
  }

  return(TRUE)
}

