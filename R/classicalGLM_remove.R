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


  #apply mask to data
  if(!is.null(mask)) {
    mask <- as.logical(mask)
    for(s in 1:n_sess){
      data[[s]]$BOLD <- data[[s]]$BOLD[,mask]
    }
  } else {
    mask <- rep(TRUE, V)
  }







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

