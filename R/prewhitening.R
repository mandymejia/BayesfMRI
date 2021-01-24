#' Compute inverse covariance matrix for AR process (up to a constant scaling factor)
#'
#' @param ar vector of p AR parameters
#' @param ntime number of time points in timeseries
#'
#' @return inverse covariance matrix for AR process (up to a constant scaling factor)
#' @export
getInvCovAR <- function(ar, ntime){
  require(Matrix)
  Inv0 <- diag(ntime)
  incr0 <- matrix(0, nrow=ntime, ncol=ntime)
  offs <- row(Inv0) - col(Inv0) #identifies each off-diagonal
  p <- length(ar)
  for(k in 1:p){
    incr <- incr0 #matrix of zeros
    incr[offs==k] <- -1*ar[k]
    Inv0 <- Inv0 + incr
  }
  Inv <- Inv0 %*% t(Inv0)
  return(Inv)
}


#' Compute square root of a symmetric matrix (such as inverse covariance mat) using SVD
#'
#' @param Inv Symmetric matrix
#'
#' @return Matrix square root of \code{Inv}
#' @keywords internal
getSqrtInv <- function(Inv){
  require(Matrix)
  #
  # Cov = UD^2U', then Inv = UD^(-2)U', and Inv_sqrt = UD^(-1)U'
  ei <- eigen(Inv)
  d2inv <- ei$values #diag(D^(-2))
  if(sum(d2inv<0)>0) print('negative eigenvalues')
  sqrtInv <- ei$vectors %*% diag(sqrt(d2inv)) %*% t(ei$vectors)
  return(sqrtInv)
}

# data <- session_data
# ar_order <- 6
# surface_FWHM <- 5
# surface_sigma <- surface_FWHM / (2*sqrt(2*log(2)))
# cifti_data <- cifti_ss
# hemisphere <- hem

#' Prewhiten cifti session data
#'
#' @param data List of sessions (see \code{is.session} and \code{is.session_pw})
#' @param scale_BOLD (logical) Should the BOLD response be scaled? (Default is TRUE)
#' @param scale_design (logical) Should the design matrix be scaled? (Default is TRUE)
#' @param ar_order Order of the AR used to prewhiten the data at each location
#' @param surface_sigma range parameter for smoothing. Remember that
#'   sigma = \code{FWHM / (2*sqrt(2*log(2))})
#' @param cifti_data A \code{xifti} object used to map the AR coefficient
#'   estimates onto the surface mesh for smoothing.
#' @param wb_path Relative path to the Human Connectome Project workbench
#' @param hemisphere 'left' or 'right'
#' @importFrom Matrix bandSparse sparseMatrix bdiag
#' @importFrom ciftiTools smooth_cifti
#' @importFrom stats ar.yw
#'
#' @return The prewhitened data (in a list), the smoothed, averaged AR
#'   coefficient estimates used in the prewhitening, the smoothed, average
#'   residual variance after prewhitening, and the value given for \code{ar_order}.
#' @export
prewhiten_cifti <- function(data, scale_BOLD = TRUE, scale_design = TRUE, ar_order = 6, surface_sigma = NULL, cifti_data, wb_path = NULL, hemisphere = NULL) {
  require(Matrix)
  require(ciftiTools)
  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  is_missing <- is.na(data[[1]]$BOLD[1,])
  if(length(is_missing) > 0) {
    cat("Some response locations are missing. Prewhitening will be done for",
        sum(!is_missing),"of", length(is_missing),"data locations.\n")
    for(s in 1:n_sess) {
      data[[s]]$BOLD <- data[[s]]$BOLD[,!is_missing]
    }
  }

  V <- ncol(data[[1]]$BOLD) #number of data locations
  K <- ncol(data[[1]]$design) #number of tasks
  ntime <- nrow(data[[1]]$BOLD) # Number of time steps
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }


  GLM_result <- vector('list', length=n_sess)
  names(GLM_result) <- session_names
  AR_coeffs <- array(dim = c(V,ar_order,n_sess))
  AR_resid_var <- array(dim = c(V,n_sess))
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
    XTX_inv <- try(solve(crossprod(X_reg)))
    if("try-error" %in% class(XTX_inv)) {
      stop("There is some numerical instability in your design matrix
                    (due to very large or very small values). Scaling the design
                    matrix is suggested.")
    } else {
      # GLM_result[[s]] <- t(XTX_inv %*% t(X_reg) %*% y_reg)
      betas <- t(XTX_inv %*% t(X_reg) %*% y_reg)
      resids <- y_reg - tcrossprod(X_reg,betas)
      for(v in 1:V) {
        if(is.na(resids[1,v])) next
        ar_v <- ar.yw(resids[,v],aic = FALSE,order.max = ar_order)
        # aic_order <- ar.yw(resids[,v])$order # This should be the order
        # of the time series if the AIC is used to select the best order
        AR_coeffs[v,,s] <- ar_v$ar # The AR parameter estimates
        AR_resid_var[v,s] <- ar_v$var.pred # Resulting variance
      }
    }
    data[[s]]$BOLD <- y_reg
    data[[s]]$design <- X_reg
  }

  avg_AR <- apply(AR_coeffs,1:2, mean)
  avg_var <- apply(as.matrix(AR_resid_var),1,mean)
  if(!is.null(surface_sigma) & !is.null(cifti_data) & !is.null(hemisphere) & !is.null(wb_path)) {
    cat("Smoothing AR coefficients and residual variance...")
    rows.keep <- which(!is.na(avg_AR[,1]))
    avg_xifti <- cifti_data
    avg_xifti$data[[paste0("cortex_",hemisphere)]] <- avg_AR[rows.keep,]
    smooth_avg_xifti <- smooth_cifti(avg_xifti,surface_sigma = surface_sigma,
                                     volume_sigma = surface_sigma,wb_path = wb_path)
    avg_AR[rows.keep,] <- smooth_avg_xifti$data[[paste0("cortex_",hemisphere)]]
    var_xifti <- cifti_data
    var_xifti$data[[paste0("cortex_",hemisphere)]] <- as.matrix(avg_var[rows.keep])
    smooth_var_xifti <- smooth_cifti(var_xifti,surface_sigma = surface_sigma,
                                     volume_sigma = surface_sigma,wb_path = wb_path)
    avg_var[rows.keep] <- smooth_var_xifti$data[[paste0("cortex_",hemisphere)]]
    cat("done!\n")
  }
  # Create the sparse pre-whitening matrix
  off_diags <- row(diag(ntime)) - col(diag(ntime))

  # Initialize the block diagonal covariance matrix
  template_pw <- Matrix::bandSparse(n = ntime,
                                    k = 0:(ar_order + 1),
                                    symmetric = TRUE)
  template_pw_list <- rep(list(template_pw),V)
  rows.rm <- numeric()
  cat("Prewhitening...\n")
  for(v in 1:V) {
    if(v %% 100 == 0) cat("Location",v,"of",V,"\n")
    rows <- (1:ntime) + ntime*(v-1)
    rows.rm <- c(rows.rm,rows[1:(ar_order + 1)], rows[(ntime - ar_order):ntime])
    ar.v <- avg_AR[v,]
    var.v <- avg_var[v]
    if(is.na(ar.v[1])) {
      template_pw_list[[v]] <- Matrix::sparseMatrix(i=NULL, j=NULL, dims=c(ntime, ntime))
    } else {
      Inv.v <- getInvCovAR(ar.v, ntime) #inverse correlation (sort of)
      Dinv.v <- diag(rep(1/sqrt(var.v), ntime)) #inverse diagonal SD matrix
      Inv.v <- Dinv.v %*% Inv.v %*% Dinv.v #inverse covariance
      sqrtInv.v <- getSqrtInv(Inv.v) #find V^(1) and take sqrt
      #after the (p+1)th off-diagonal, values of sqrtInv.v very close to zero
      diags <- list(); length(diags) <- ar_order+2
      for(k in 0:(ar_order+1)){
        diags[[k+1]] <- sqrtInv.v[off_diags==k]
      }
      matband <- bandSparse(ntime, k=0:(ar_order+1), symm=TRUE, diag=diags)
      template_pw_list[[v]] <- matband
    }
  }
  cat("done!\n")

  sqrtInv_all <- Matrix::bdiag(template_pw_list)


  # prewhite_mat <- Reduce(rbind,template_pw_list)
  pw_data <- sapply(data,function(data_s) {
    bold_out <- matrix(NA,ntime, length(is_missing))
    pw_BOLD <- as.vector(sqrtInv_all %*% c(data_s$BOLD))
    bold_out[,!is_missing] <- pw_BOLD
    # pw_bold <- matrix(pw_BOLD,ntime,V)
    # all_design <- Reduce(rbind,rep(list(data_s$design),V))
    all_design <- bdiag(rep(list(data_s$design),V))
    pw_design <- sqrtInv_all %*% all_design
    return(list(BOLD = bold_out, design = pw_design))
  }, simplify = F)

  # return(GLM_result)
  return(list(data = pw_data,
              AR_coeffs = avg_AR,
              AR_var = avg_var,
              ar_order = ar_order))
}
