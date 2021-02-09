#' Fit classical massive univariate GLM to task fMRI data after prewhitening
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#'
#' @return A list of classical GLM task activation estimates, where each element represents a session.
#'
#' @export
classicalGLM_pw <- function(data, scale_BOLD=TRUE, scale_design = TRUE){

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  V <- ncol(data[[1]]$BOLD) #number of data locations
  is_missing <- is.na(data[[1]]$BOLD[1,])
  K <- ncol(data[[1]]$design) / sum(!is_missing) #number of tasks

  # Come back to these checks ----
  # for(s in 1:n_sess){
  #   if(! is.session_pw(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
  #   if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
  #   if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  # }
  # >> /Come back to these checks ----

  GLM_result <- vector('list', length=n_sess)
  names(GLM_result) <- session_names
  for(s in 1:n_sess){

    #scale data to represent % signal change (or just center to remove baseline if scale=FALSE)
    BOLD_s <- data[[s]]$BOLD
    if(sum(is_missing) > 0) BOLD_s <- BOLD_s[,!is_missing]
    # BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD, transpose = FALSE)
    if(scale_design) {
      design_s <- scale_design_mat(as.matrix(data[[s]]$design))
    } else {
      # design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline
    }

    #regress nuisance parameters from BOLD data and design matrix
    if('nuisance' %in% names(data[[s]])){
      design_s <- data[[s]]$design
      nuisance_s <- data[[s]]$nuisance
      y_reg <- nuisance_regress(BOLD_s, nuisance_s)
      X_reg <- nuisance_regress(design_s, nuisance_s)
    } else {
      y_reg <- BOLD_s
      X_reg <- data[[s]]$design
    }
    y_reg <- c(y_reg)
    # X_reg <- bdiag(rep(list(X_reg),V))
    # XTX_inv <- try(solve(t(X_reg) %*% X_reg))
    XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_reg)))
    if("try-error" %in% class(XTX_inv)) {
      stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
    } else {
      if(sum(is_missing) > 0) {
        GLM_result[[s]] <- matrix(NA, K, V)
        # cat("V=",V,"K=",K,"\n")
        coef_output <- matrix(as.numeric(XTX_inv %*% t(X_reg) %*% y_reg), byrow = F)
        # cat("length(coef_output)=",length(coef_output),"\n")
        # print(table(is_missing))
        # GLM_result[[s]][!is_missing,] <- coef_output
        GLM_result[[s]][,!is_missing] <- coef_output
        GLM_result[[s]] <- t(GLM_result[[s]])
      } else {
        GLM_result[[s]] <- matrix(as.numeric(t(XTX_inv %*% t(X_reg) %*% y_reg)),V,K, byrow = T)
      }
    }
  }

  return(GLM_result)

}

#' Is session? (Prewhitened data)
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
is.session_pw <- function(sess){

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

    if(! (is.matrix(sess$design) | is.list(sess$design))){stop('I expected design to be a matrix or a list, but it is not')}
    if(is.list(sess$design)) {
      all_list_classes <- sapply(sess$design,class)
      if(!all(all_list_classes == 'matrix' | all_list_classes == 'dgeMatrix')) stop("If the design is a list, every list element should be a matrix.")
    }

    if(! (is.matrix(sess$nuisance))){stop('I expected nuisance to be a matrix, but it is not')}


    ## check the dimensions of each field: T
    if(nrow(sess$BOLD) != nrow(sess$design)){stop("BOLD and design don't have the same number of rows (time points)")}
    if(nrow(sess$BOLD) != nrow(sess$nuisance)){stop("BOLD and nuisance don't have the same number of rows (time points)")}
    if(nrow(sess$design) != nrow(sess$nuisance)){stop("design and nuisance don't have the same number of rows (time points)")}

    ## BOLD, design
  } else if (num_fields==2){

    ## check identities of fields
    fields = c('BOLD','design')
    if(! all.equal(names(sess),fields)){stop(
      paste0('You are missing the following fields',setdiff(names(sess),fields)))}

    ## check each field's type
    if(! (is.numeric(sess$BOLD))){stop('I expected BOLD to be numeric, but it is not')}
    if(! (is.matrix(sess$BOLD))){stop('I expected BOLD to be a matrix, but it is not')}

    if(! class(sess$design) %in% c("matrix","dgCMatrix") ){stop('I expected design to be a matrix, but it is not')}

    ## check the dimensions of each field: T
    if(is.matrix(sess$design)) {
      if(nrow(sess$BOLD) != nrow(sess$design)){stop("BOLD and design don't have the same number of rows (time points)")}
    } else {
      if(! (nrow(sess$BOLD) * sum(!is.na(sess$BOLD[1,]))) == nrow(sess$design))
        stop("If prewhitening has been performed, the number of rows in the design matrix is not equal to the product of the number of data locations and the number of time points.")
    }


  } else {
    stop('I expected the session to have 2 or 3 fields, but it does not.')
  }

  return(TRUE)
}

