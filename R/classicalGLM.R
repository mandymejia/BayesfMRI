#' Fit classical massive univariate GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams avg_sessions_Param
#'
#' @return A list of classical GLM task activation estimates, where each element represents a session.
#'
#' @importFrom Matrix solve crossprod
#' @export
classicalGLM <- function(data,
                         scale_BOLD=TRUE,
                         scale_design = TRUE,
                         avg_sessions = TRUE){

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  # Name sessions and check compatibility of multi-session arguments
  session_names <- names(data)
  n_sess <- length(session_names)

  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE

  V <- ncol(data[[1]]$BOLD) #number of data locations
  is_missing <- is.na(data[[1]]$BOLD[1,])
  V2 <- sum(!is_missing)

  is_pw <- !is.matrix(session_data$single_session$design) #if prewhitening has been done, design is a large sparse matrix (class dgCMatrix)
  if(is_pw){
    if(scale_BOLD | scale_design) stop('If data is prewhitened, scale_BOLD and scale_design should be FALSE. Setting both to FALSE.')
    scale_BOLD <- FALSE
    scale_design <- FALSE
  }

  #checks
  if(!is_pw){
    K <- ncol(data[[1]]$design) #number of tasks
    for(s in 1:n_sess){
      if(!is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
      if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
      if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }
  }
  if(is_pw){
    K <- ncol(data[[1]]$design) / sum(!is_missing) #number of tasks
    for(s in 1:n_sess){
      if(!is.session_pw(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
      if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
      if(ncol(data[[s]]$design) != K*V2) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }
  }

  if(avg_sessions) num_GLM <- n_sess + 1
  if(!avg_sessions) num_GLM <- n_sess

  GLM_result <- vector('list', length=num_GLM)
  names(GLM_result) <- session_names

  if(avg_sessions){
    y_reg_all <- X_reg_all <- NULL
  }

  for(s in 1:num_GLM){

    if(s <= n_sess){

      #scale data to represent % signal change (or just center to remove baseline if scale=FALSE)
      BOLD_s <- data[[s]]$BOLD
      if(sum(is_missing) > 0) BOLD_s <- BOLD_s[,!is_missing]
      BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD, transpose = FALSE)
      if(scale_design) {
        design_s <- scale_design_mat(data[[s]]$design)
      } else {
        #only center if no prewhitening (prewhitening centers data)
        if(!is_pw) design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline
        if(is_pw) design_s <- data[[s]]$design
      }

      #regress nuisance parameters from BOLD data and design matrix (only for non-pw data, since pw data has already been nuisance regressed)
      do_nuisance <- ('nuisance' %in% names(data[[s]]))
      if(is_pw & do_nuisance) stop('Prewhitened data should not include nuisance. Contact developer.')
      if(do_nuisance){
        nuisance_s <- data[[s]]$nuisance
        y_reg <- nuisance_regress(BOLD_s, nuisance_s)
        X_reg <- nuisance_regress(design_s, nuisance_s)
      } else {
        y_reg <- BOLD_s
        X_reg <- design_s
      }

      if(avg_sessions){
        y_reg_all <- rbind(y_reg_all, y_reg)
        X_reg_all <- rbind(X_reg_all, X_reg)
      }
    }

    if(s == n_sess + 1){
      y_reg <- y_reg_all
      X_reg <- X_reg_all
    }

    # ESTIMATE MODEL COEFFICIENTS

    beta_hat_s <- SE_beta_hat_s <- matrix(NA, K, V)
    if(is_pw){
      y_reg <- c(y_reg) #make y a vector
      XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_reg)))
      if("try-error" %in% class(XTX_inv)) {
        stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      }
      coef_s <- as.matrix(XTX_inv %*% t(X_reg) %*% y_reg) #a vector of (beta1, beta2, .. for location 1, beta 1, beta 2, ... for location 2, ...)
      beta_hat_s[,!is_missing] <- coef_s #RHS is a vector
      resid_s <- matrix(y_reg - X_reg %*% coef_s, ncol = V2)
    }

    if(!is_pw){
      XTX_inv <- try(solve(t(X_reg) %*% X_reg))
      if("try-error" %in% class(XTX_inv)) stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      coef_s <- t(XTX_inv %*% t(X_reg) %*% y_reg)
      beta_hat_s[,!is_missing] <- coef_s #RHS is a matrix
      resid_s <- y_reg - X_reg %*% coef_s
    }



    # ESTIMATE STANDARD ERRORS OF ESTIIMATES

    ntime <- length(y_reg) / V2
    if(round(ntime) != ntime) stop('Error, contact developer')

    #compute residual SD
    var_error <- colVars(resid_s) * (ntime - 1) / (ntime - K - 1)
    if(is_pw) var_error <- rep(mean(var_error), length(var_error))
    sd_error <- sqrt(var_error)

    #compute SE of betas
    SE_beta_s <- sqrt(diag(XtX_inv)) * rep(sd_error, each = K)
    SE_beta_hat_s[,!is_missing] <- SE_beta_s

    #DoF
    DOF <- ntime - K - 1




    GLM_result[[s]] <- list(estimates = beta_hat_s,
                            SE_estimates = SE_beta_hat_s,
                            DOF = DOF)
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
