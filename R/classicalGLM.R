#' Fit classical massive univariate GLM to task fMRI data
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
classicalGLM <- function(data, scale_BOLD=TRUE, scale_design = TRUE){

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  V <- ncol(data[[1]]$BOLD) #number of data locations
  K <- ncol(data[[1]]$design) #number of tasks
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  GLM_result <- vector('list', length=n_sess)
  names(GLM_result) <- session_names
  for(s in 1:n_sess){

    #scale data to represent % signal change (or just center to remove baseline if scale=FALSE)
    BOLD_s <- data[[s]]$BOLD
    BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD, transpose = FALSE)
    if(scale_design) {
      design_s <- scale_design_mat(data[[s]]$design)
    } else {
      design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline
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
    XTX_inv <- try(solve(t(X_reg) %*% X_reg))
    if("try-error" %in% class(XTX_inv)) {
      stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
    } else {
      GLM_result[[s]] <- t(XTX_inv %*% t(X_reg) %*% y_reg)
    }
  }

  return(GLM_result)

}

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
        GLM_result[[s]] <- matrix(NA, V, K)
        # cat("V=",V,"K=",K,"\n")
        coef_output <- matrix(as.numeric(t(XTX_inv %*% t(X_reg) %*% y_reg)))
        # cat("length(coef_output)=",length(coef_output),"\n")
        # print(table(is_missing))
        GLM_result[[s]][!is_missing,] <- coef_output
      } else {
        GLM_result[[s]][is] <- matrix(as.numeric(t(XTX_inv %*% t(X_reg) %*% y_reg)),V,K, byrow = T)
      }
    }
  }

  return(GLM_result)

}
