#' Fit classical massive univariate GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param mask (Optional) A length \eqn{V} logical vector indicating if each
#'  vertex is to be included.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams avg_sessions_Param
#' @param num_permute The number of permutations that should be performed
#'   in order to allow for permutation testing to determine activations. A value
#'   of 0 will not perform any permutations, and permutation testing will not
#'   be available.
#'
#' @return A list of lists containing classical GLM task activation estimates, standard error estimates, and degrees of freedom. Each list element represents a session.
#'
#' @importFrom matrixStats colVars
#' @export
classicalGLM <- function(data, mask = NULL, scale_BOLD=TRUE, scale_design = TRUE, avg_sessions = TRUE, num_permute = 0){
  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  # Name sessions and check compatibility of multi-session arguments
  session_names <- names(data)
  n_sess <- length(session_names)
  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE

  do_permute <- num_permute > 0

  is_pw <- !is.matrix(data[[1]]$design) #if prewhitening has been done, design is a large sparse matrix (class dgCMatrix)
  if(is_pw){
    if(scale_BOLD | scale_design) warning('If data is prewhitened, scale_BOLD and scale_design should be FALSE. Setting both to FALSE.')
    scale_BOLD <- FALSE
    scale_design <- FALSE
  }
  #check dimensions
  V <- ncol(data[[1]]$BOLD)
  for(s in 1:n_sess){
    if(!is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
  }

  #ID any zero-variance voxels and remove from analysis
  zero_var <- sapply(data, function(x){
    x$BOLD[is.na(x$BOLD)] <- 0 #to detect medial wall locations coded as NA
    x$BOLD[is.nan(x$BOLD)] <- 0 #to detect medial wall locations coded as NaN
    vars <- matrixStats::colVars(x$BOLD)
    return(vars < 1e-6)
  })
  zero_var <- (rowSums(zero_var) > 0) #check whether any vertices have zero variance in any session

  #remove zero var locations from mask
  if(sum(zero_var) > 0){
    if(is.null(mask)) num_flat <- sum(zero_var) else num_flat <- sum(zero_var[mask==1])
    #if(num_flat > 1) warning(paste0('I detected ', num_flat, ' vertices that are flat (zero variance), NA or NaN in at least one session. Removing these from analysis. See mask returned with function output.'))
    #if(num_flat == 1) warning(paste0('I detected 1 vertex that is flat (zero variance), NA or NaN in at least one session. Removing it from analysis. See mask returned with function output.'))
    mask_orig <- mask
    if(!is.null(mask)) mask[zero_var==TRUE] <- 0
    if(is.null(mask)) mask <- !zero_var
  } else {
    mask_orig <- NULL
  }

  #apply mask to data
  if(is.null(mask)){
    mask_use <- rep(TRUE, V)
  } else {
    mask_use <- as.logical(mask)
  }
  V_all <- length(mask_use)
  V <- sum(mask_use)
  for(s in 1:n_sess){
    data[[s]]$BOLD <- data[[s]]$BOLD[,mask_use]
  }

  #check dimensions
  if(!is_pw){
    K <- ncol(data[[1]]$design) #number of tasks
    for(s in 1:n_sess){
      if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }
  }
  if(is_pw){
    K <- ncol(data[[1]]$design) / V #number of tasks
    for(s in 1:n_sess){
      if(ncol(data[[s]]$design) != K*V) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
    }
  }

  if(avg_sessions) {
    num_GLM <- n_sess + 1
    session_names <- append(session_names,"avg")
  }
  if(!avg_sessions) num_GLM <- n_sess
  GLM_result <- vector('list', length=num_GLM)
  names(GLM_result) <- session_names
  if(avg_sessions){
    y_reg_all <- X_reg_all <- NULL
    if(do_permute) permuted_y_all <- NULL
  }
  for(s in 1:num_GLM){
    if(s <= n_sess){
      #scale data to represent % signal change (or just center to remove baseline if scale=FALSE)
      BOLD_s <- data[[s]]$BOLD #[,!is_missing]
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
      if(do_permute) {
        permuted_y <- sapply(seq(num_permute), function(m) {
          ntime_m <- nrow(BOLD_s)
          perm_idx <- sample(seq(ntime_m),size = ntime_m,replace = F)
          perm_y <- BOLD_s[perm_idx,]
          # shift_size <- sample(-ntime_m:ntime_m, size = 1)
          # perm_y <- BOLD_s[((seq(ntime_m) + shift_size) %% ntime_m) + 1,]
          return(perm_y)
        }, simplify = F)
      }
      if(avg_sessions){
        y_reg_all <- cbind(y_reg_all, y_reg)
        X_reg_all <- rbind(X_reg_all, X_reg)
        if(do_permute) {
          if(!is.null(permuted_y_all))
            permuted_y_all <- mapply(rbind, permuted_y_all, permuted_y, SIMPLIFY = F)
          if(is.null(permuted_y_all)) permuted_y_all <- permuted_y
        }
      }
    }
    #for average, analyze time-concatenated data and design
    if(s == n_sess + 1){
      y_reg <- y_reg_all
      X_reg <- X_reg_all
      if(do_permute) permuted_y <- permuted_y_all
    }
    # ESTIMATE MODEL COEFFICIENTS
    beta_hat_s <- SE_beta_hat_s <- matrix(NA, K, V_all)
    beta_hat_null <- SE_beta_hat_null <- NULL
    if(do_permute) {
      beta_hat_null <- SE_beta_hat_null <- array(NA, dim = c(K,V_all,num_permute))
    }
    ntime <- length(y_reg) / V
    if(round(ntime) != ntime) stop('Error, contact developer')
    #compute residual SD
    DOF <- ntime - K - 1
    if(is_pw){
      y_reg <- c(y_reg) #make y a vector (grouped by location)
      XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_reg)))
      if("try-error" %in% class(XTX_inv)) {
        stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      }
      XTX_invXT <- Matrix::tcrossprod(XTX_inv,X_reg)
      coef_s <- as.matrix(XTX_invXT %*% y_reg) #a vector of (estimates for location 1, estimates for location 2, ...)
      beta_hat_s[,mask_use==TRUE] <- coef_s #RHS is a vector
      resid_s <- matrix(y_reg - X_reg %*% coef_s, ncol = V)
      if(do_permute) {
        cat("Performing permutations for session",session_names[s],"\n")
        for(m in seq(num_permute)) {
          y_perm <- c(permuted_y[[m]])
          beta_hat_null[,mask_use,m] <- as.matrix(XTX_invXT %*% y_perm)
          fit_perm <- X_reg %*% c(beta_hat_null[,mask_use,m])
          resid_perm <- matrix(y_perm - fit_perm, ncol = V)
          var_error_perm <- matrixStats::colVars(resid_perm) * (ntime - 1) / DOF #correct for DOF
          var_error_perm <- rep(mean(var_error_perm), length(var_error_perm))
          sd_error_perm <- sqrt(var_error_perm)
          #compute SE of betas
          SE_beta_hat_null[,mask_use,m] <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error_perm, each = K)
        }
      }
    }
    if(!is_pw){
      XTX_inv <- try(solve(t(X_reg) %*% X_reg))
      if("try-error" %in% class(XTX_inv)) stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      XTX_invXT <- tcrossprod(XTX_inv,X_reg)
      coef_s <-XTX_invXT %*% y_reg
      beta_hat_s[,mask_use==TRUE] <- coef_s   #RHS is a matrix
      resid_s <- y_reg - X_reg %*% coef_s
      if(do_permute) {
        cat("Performing permutations for session",session_names[s],"\n")
        for(m in seq(num_permute)) {
          y_perm <- c(permuted_y[[m]])
          beta_hat_null[,mask_use,m] <- as.matrix(XTX_invXT %*% y_perm)
          fit_perm <- X_reg %*% c(beta_hat_null[,mask_use,m])
          resid_perm <- matrix(y_perm - fit_perm, ncol = V)
          var_error_perm <- matrixStats::colVars(resid_perm) * (ntime - 1) / DOF #correct for DOF
          var_error_perm <- rep(mean(var_error_perm), length(var_error_perm))
          sd_error_perm <- sqrt(var_error_perm)
          #compute SE of betas
          SE_beta_hat_null[,mask_use,m] <- matrix(sqrt(diag(XTX_inv)), nrow=K, ncol=V) * matrix(sd_error_perm, nrow=K, ncol=V, byrow = TRUE)
        }
      }
    }
    # ESTIMATE STANDARD ERRORS OF ESTIMATES
    var_error <- matrixStats::colVars(resid_s) * (ntime - 1) / DOF #correct for DOF
    if(is_pw) var_error <- rep(mean(var_error), length(var_error))
    sd_error <- sqrt(var_error)
    #compute SE of betas
    if(is_pw) SE_beta_s <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, each = K)
    if(!is_pw) SE_beta_s <- matrix(sqrt(diag(XTX_inv)), nrow=K, ncol=V) * matrix(sd_error, nrow=K, ncol=V, byrow = TRUE)
    SE_beta_hat_s[,mask_use==TRUE] <- SE_beta_s

    GLM_result[[s]] <- list(estimates = t(beta_hat_s),
                            SE_estimates = t(SE_beta_hat_s),
                            null_estimates = aperm(beta_hat_null, perm = c(2,1,3)),
                            null_SE_estimates = aperm(SE_beta_hat_null, perm = c(2,1,3)),
                            DOF = DOF,
                            mask = mask,
                            mask_orig = mask_orig)
  }
  class(GLM_result) <- 'classicalGLM'
  return(GLM_result)
}

