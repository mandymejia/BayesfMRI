#' Fit classical massive univariate GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param beta_names (Optional) Names of tasks represented in design matrix
#' @param vertices (Only required for AR smoothing) A \eqn{V x 3} matrix, where
#' each row contains the Euclidean coordinates at which a given vertex in the
#' mesh is located. \eqn{V} is the number of vertices in the mesh. Should provide
#' vertices and faces or mesh, but not both.
#' @param faces (Only required for AR smoothing) An \eqn{F x 3} matrix, where
#' each row contains the vertex indices for a given triangular face in the mesh.
#' \eqn{F} is the number of faces in the mesh. Should provide vertices and faces
#' or mesh, but not both.
#' @param mesh (Only required for AR smoothing) An \code{"inla.mesh"} object
#' (see \code{\link{make_mesh}} for surface data). Should provide mesh or vertices
#' and faces, but not both.
#' @param mask (Optional) A length \eqn{V} logical vector indicating if each
#'  vertex is to be included.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#' @inheritParams num.threads_Param
#' @inheritParams avg_sessions_Param
#'
#' @return A list of lists containing classical GLM task activation estimates, standard error estimates, and degrees of freedom. Each list element represents a session.
#'
#' @importFrom matrixStats colVars
#' @export
classicalGLM <- function(data,
                         beta_names = NULL,
                         vertices = NULL,
                         faces = NULL,
                         mesh = NULL,
                         mask = NULL,
                         scale_BOLD=TRUE,
                         scale_design = TRUE,
                         ar_order = 6,
                         ar_smooth = 5,
                         num.threads=4,
                         avg_sessions = TRUE){

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)
  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE

  is_pw <- FALSE # This will need to be removed throughout
  mask_orig <- mask # This will also need to be removed throughout

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')


  # is_pw <- !is.matrix(data[[1]]$design) #if prewhitening has been done, design is a large sparse matrix (class dgCMatrix)
  # if(is_pw){
  #   if(scale_BOLD | scale_design) warning('If data is prewhitened, scale_BOLD and scale_design should be FALSE. Setting both to FALSE.')
  #   scale_BOLD <- FALSE
  #   scale_design <- FALSE
  # }

  V <- ncol(data[[1]]$BOLD) #number of data locations
  ntime <- nrow(data[[1]]$BOLD)
  K <- ncol(data[[1]]$design) #number of tasks

  if(!is.null(beta_names)){
    if(length(beta_names) != K) stop(paste0('I detect ', K, ' task based on the design matrix, but the length of beta_names is ', length(beta_names), '.  Please fix beta_names.'))
  }

  if(is.null(beta_names)){
    beta_names_maybe <- colnames(data[[1]]$design) #grab beta names from design (if not provided)
    if(!is.null(beta_names_maybe)) beta_names <- beta_names_maybe
    if(is.null(beta_names_maybe)) beta_names <- paste0('beta',1:K)
  }

  #check dimensions
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  #apply mask to data
  if(!is.null(mask)) {
    mask <- as.logical(mask)
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
  } else {
    mask <- rep(TRUE, V)
  }
  V <- sum(mask)
  V_all <- length(mask)

  #make sure no invalid vertices exist after masking
  tmp <- make_mask(data)
  if(!is.null(tmp)) stop('Flat/NA/NaN vertices exist outside the mask. Update mask to exclude these.')

  #collect data and design matrices
  design <- vector('list', length=n_sess)

  #compute AR coefficients and average over sessions
  if(is.null(ar_order)) ar_order <- 0
  ar_order <- as.numeric(ar_order)
  do_pw <- (ar_order > 0)
  if(do_pw){
    AR_coeffs <- array(dim = c(V,ar_order,n_sess))
    AR_resid_var <- array(dim = c(V,n_sess))
  }
  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

    #scale data to represent % signal change (or just center if scale=FALSE)
    BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale_BOLD)
    if(scale_design) design_s <- scale_design_mat(data[[s]]$design)
    if(!scale_design) design_s <- scale(data[[s]]$design, scale = F)
    design[[s]] <- design_s #after scaling but before nuisance regression

    #regress nuisance parameters from BOLD data and design matrix
    if('nuisance' %in% names(data[[s]])){
      nuisance_s <- data[[s]]$nuisance
      data[[s]]$BOLD <- nuisance_regression(BOLD_s, nuisance_s)
      data[[s]]$design <- nuisance_regression(design_s, nuisance_s)
    } else {
      data[[s]]$BOLD <- BOLD_s
      data[[s]]$design <- design_s
    }

    #estimate prewhitening parameters
    if(do_pw){
      resids <- nuisance_regression(data[[s]]$BOLD, data[[s]]$design)
      AR_est <- pw_estimate(resids, ar_order)
      AR_coeffs[,,s] <- AR_est$phi
      AR_resid_var[,s] <- AR_est$sigma_sq
    }
  }

  #average prewhitening parameters across sessions
  if(do_pw){
    avg_AR <- apply(AR_coeffs, 1:2, mean)
    avg_var <- apply(as.matrix(AR_resid_var), 1, mean)
  }

  #smooth prewhitening parameters
  if(is.null(ar_smooth)) ar_smooth <- 0
  if(ar_smooth > 0){

    #check that only mesh OR vertices+faces supplied
    has_mesh <- !is.null(mesh)
    has_verts_faces <- !is.null(vertices) & !is.null(faces)
    has_howmany <- has_mesh + has_verts_faces
    if(has_howmany != 1) stop('Must supply either mesh or vertices and faces for AR smoothing.')

    #if mesh provided, grab vertices and faces
    if(has_mesh){
      vertices <- mesh$loc
      faces <- mesh$graph$tv
    }

    #smooth prewhitening parameters
    if(is.null(ar_smooth)) ar_smooth <- 0
    if(ar_smooth > 0) {
      AR_smoothed_list <- pw_smooth(vertices=vertices,
                                    faces=faces,
                                    AR=avg_AR,
                                    var=avg_var,
                                    FWHM=ar_smooth)
      avg_AR <- AR_smoothed_list$AR
      avg_var <- AR_smoothed_list$var
    }
  }

  if(avg_sessions) {
    num_GLM <- n_sess + 1
    session_names <- append(session_names,"avg")
  }
  if(!avg_sessions) num_GLM <- n_sess
  GLM_result <- vector('list', length=num_GLM)
  names(GLM_result) <- session_names


  #If not prewhitening, do matrix operations to analyze all vertices simultaneously
  if(!do_pw){

    if(avg_sessions) y_reg_all <- X_reg_all <- NULL

    #if analyzing a session
    for(s in 1:num_GLM){
      if(s <= n_sess){
        y_reg <- data[[s]]$BOLD
        X_reg <- data[[s]]$design
        if(avg_sessions){
          y_reg_all <- rbind(y_reg_all, y_reg)
          X_reg_all <- rbind(X_reg_all, X_reg)
        }
      }
      #for average, analyze time-concatenated data and design
      if(s == n_sess + 1){ y_reg <- y_reg_all; X_reg <- X_reg_all }
      ntime_s <- nrow(X_reg)

      # ESTIMATE MODEL COEFFICIENTS
      beta_hat_s <- SE_beta_hat_s <- matrix(NA, K, V_all)
      XTX_inv <- try(solve(t(X_reg) %*% X_reg))
      if("try-error" %in% class(XTX_inv)) stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      coef_s <- XTX_inv %*% t(X_reg) %*% y_reg
      beta_hat_s[,mask==TRUE] <- coef_s   #RHS is a matrix
      resid_s <- y_reg - X_reg %*% coef_s

      # ESTIMATE STANDARD ERRORS OF ESTIIMATES
      #compute residual SD
      DOF <- ntime_s - K - 1
      var_error <- matrixStats::colVars(resid_s) * (ntime_s - 1) / DOF #correct for DOF
      #if(is_pw) var_error <- rep(mean(var_error), length(var_error))
      sd_error <- sqrt(var_error)
      #compute SE of betas
      #if(is_pw) SE_beta_s <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, each = K)
      #if(!is_pw)
      SE_beta_s <- matrix(sqrt(diag(XTX_inv)), nrow=K, ncol=V) * matrix(sd_error, nrow=K, ncol=V, byrow = TRUE)
      SE_beta_hat_s[,mask==TRUE] <- SE_beta_s
    }
  }

  #HERE ----

  if(do_pw){

    #loop over voxels, do prewhitening, compute beta-hat and SE-beta-hat
    #do concatenation for averaging over sessions within each voxel

  }





    GLM_result[[s]] <- list(estimates = t(beta_hat_s),
                            SE_estimates = t(SE_beta_hat_s),
                            DOF = DOF,
                            mask = mask,
                            mask_orig = mask_orig)
  }
  class(GLM_result) <- 'classicalGLM'
  return(GLM_result)
}

