#' BayesGLM
#'
#' Performs spatial Bayesian GLM for task fMRI activation
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param data A list of sessions in the \code{"BfMRI.sess"} object format. Each
#'  session is a list with elements \code{"BOLD"}, \code{"design"}, and
#'  optionally \code{"nuisance"}. Each element should be a numeric matrix with
#'  \eqn{T} rows. The name of each element in \code{data} is the name of that
#'  session. Colnames of design matrix represent field names and must match across
#'  sessions. See \code{?is.BfMRI.sess} for details.
#'
#' @param design_multiple (Optional) A list of \eqn{T \times K \times D} arrays of \eqn{D}
#' different design matrices for model comparison.
#' @param vertices,faces For cortical surface data, the geometry is based on
#'  the \code{vertices} and \code{faces} arguments (must provide both).
#'
#'  \code{vertices} is a \eqn{V \times 3} matrix, where each row contains the
#'  Euclidean coordinates at which a given vertex in the mesh is located.
#'  \eqn{V} is the number of vertices in the mesh.
#'
#'  \code{faces} is a \eqn{F \times 3} matrix, where each row contains the
#'  vertex indices for a given triangular face in the mesh. \eqn{F} is the
#'  number of faces in the mesh.
#' @param labels For volumetric data, geometry is provided through a 3D
#' array containing label values for spatially distinct ROIs to be analyzed.
#' @param nbhd_order For volumetric data, what order neighborhood around data
#' locations to keep? (0 = no neighbors, 1 = 1st-order neighbors, 2 = 1st- and
#' 2nd-order neighbors, etc.). Smaller values will provide greater computational
#' efficiency at the cost of higher variance around the edge of the data.
#' @param buffer For volumetric data, size of extra voxels layers around the
#' bounding box, in terms of voxels. Set to NULL for no buffer.
#' @param mask (Optional) A length \eqn{V} logical vector indicating the
#'  vertices to include in analysis. (Currently only for surface-based analysis)
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams Bayes_Param
# @inheritParams EM_Param
#' @inheritParams ar_order_Param
#' @inheritParams ar_smooth_Param
#' @inheritParams aic_Param
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_Param
#' @inheritParams verbose_Param
# @inheritParams combine_sessions_Param
#' @param meanTol,varTol,snrTol Tolerance for mean, variance and SNR of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for mean and variance, \code{50} for SNR. Note: \code{snrTol}
#'  currently not in use, but SNR maps are returned for visualization.
# @inheritParams emTol_Param
#'
#' @return A \code{"BayesGLM"} object: a list with elements
#'  \describe{
#'    \item{INLA_model_obj}{The full result of the call to \code{INLA::inla}.}
#'    \item{field_estimates}{The estimated coefficients for the Bayesian model.}
#'    \item{result_classical}{Results from the classical model: field estimates, field standard error estimates, residuals, degrees of freedom, and the mask.}
#'    \item{mesh}{The model mesh including only the locations analyzed, i.e. within \code{mask}, without missing values, and meeting \code{meanTol} and \code{varTol}.}
#'    \item{mesh_orig}{The original mesh provided.}
#'    \item{mask}{A mask of \code{mesh_orig} indicating the locations inside \code{mesh}.}
#'    \item{design}{The design matrix, after centering and scaling, but before any nuisance regression or prewhitening.}
#'    \item{field_names}{The names of the fields.}
#'    \item{session_names}{The names of the sessions.}
#'    \item{hyperpar_posteriors}{Hyperparameter posterior densities.}
#'    \item{theta_estimates}{Theta estimates from the Bayesian model.}
#'    \item{posterior_Sig_inv}{For joint group modeling.}
#'    \item{mu_theta}{For joint group modeling.}
#'    \item{Q_theta}{For joint group modeling.}
#'    \item{y}{For joint group modeling: The BOLD data after any centering, scaling, nuisance regression, or prewhitening.}
#'    \item{X}{For joint group modeling: The design matrix after any centering, scaling, nuisance regression, or prewhitening.}
#'    \item{prewhiten_info}{Vectors of values across locations: \code{phi} (AR coefficients averaged across sessions), \code{sigma_sq} (residual variance averaged across sessions), and AIC (the maximum across sessions).}
#'    \item{call}{match.call() for this function call.}
#'  }
#'
#' @importFrom matrixStats colVars
#' @importFrom Matrix bandSparse bdiag crossprod solve Diagonal
#' @importFrom parallel detectCores makeCluster clusterMap stopCluster
#' @importFrom stats as.formula var
#' @importFrom fMRItools is_1 nuisance_regression scale_timeseries
#'
#' @importFrom utils tail
#'
#' @importFrom methods as
#' @export
BayesGLM <- function(
  data,
  design_multiple = NULL,
  vertices = NULL,
  faces = NULL,
  mask = NULL,
  labels = NULL,
  buffer = c(1,1,3,4,4),
  nbhd_order = 1,
  # Below arguments shared with `BayesGLM_cifti`.
  #combine_sessions = TRUE,
  scale_BOLD = c("auto", "mean", "sd", "none"),
  scale_design = TRUE, #[TO DO] Delete this?  It is done by BayesGLM_cifti and could be done by the user
  Bayes = TRUE,
  #EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  num.threads = 4,
  return_INLA = c("trimmed", "full", "minimal"),
  verbose = 1,
  meanTol = 1e-6,
  varTol = 1e-6,
  snrTol = 50#, emTol = 1e-3,
  ){

  # Preliminary steps. ---------------------------------------------------------
  ## Check simple arguments.
  ## These checks are in a separate function because they are shared with
  ## `BayesGLM_cifti`.
  argChecks <- BayesGLM_argChecks(
    scale_BOLD = scale_BOLD,
    scale_design = scale_design,
    Bayes = Bayes,
    #EM = EM,
    ar_order = ar_order,
    ar_smooth = ar_smooth,
    aic = aic,
    num.threads = num.threads,
    return_INLA = return_INLA,
    verbose = verbose,
    #combine_sessions = combine_sessions,
    varTol = varTol,
    meanTol = meanTol
    #emTol = emTol
  )
  scale_BOLD <- argChecks$scale_BOLD
  do_Bayesian <- argChecks$do_Bayesian
  do_EM <- argChecks$do_EM
  do_pw <- argChecks$do_pw
  return_INLA <- argChecks$return_INLA

  do_EM <- FALSE; #emTol <- 1e-3

  ## Sessions and data dimensions. ---------------------------------------------
  if (!is.BfMRI.sess(data)) {
    stop("`data` must be a list of sessions, as described in `?is.BfMRI.sess`.")
  }
  session_names <- names(data)
  if(is.null(session_names)) stop('Session names should be provided via names(data).')

  ## Geometric inputs: check ------------------------------------------------
  is_surf <- !is.null(vertices) && !is.null(faces)
  is_vol <- !is.null(labels)
  if(is_surf + is_vol != 1) stop('Must provide vertices + faces OR labels, but not both. Please check inputs')

  #TO DO:
  #add "(ignored if classicalGLM_only = TRUE) to some params"
  #add if statements for some of code if classicalGLM_only = TRUE

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of fields

  nS <- nS_orig <- length(session_names) # number of sessions
  nK <- ncol(data[[1]]$design) # number of fields
  nV <- ncol(data[[1]]$BOLD) # number of data locations (before any masking)
  nT <- vapply(data, function(x){ nrow(x$BOLD) }, 0) # numbers of timepoints

 # if (nS == 1 && combine_sessions) combine_sessions <- FALSE

  ## Field names: check or make. -----------------------------------------------
  field_names <- colnames(data[[1]]$design)
  if(is.null(field_names)) stop('Field names should be provided through `design` element of data. See documentation for details.')

  ## Define a few return variables that may or may not be calculated. ----------
  INLA_model_obj <- hyperpar_posteriors <- Q_theta <- NULL
  field_estimates <- hyperpar_posteriors <- mu_theta <- y_all <- XA_all_list <- NULL
  theta_estimates <- Sig_inv <- mesh <- mesh_orig <- NULL

  ## Set up spatial data -------------------------------------------------------

  # [TO DO] Implement masking for volumetric analysis.  Currently within is_surf.

  mask_orig <- mask

  if(is_surf){

    # [TO DO] Add a check that if(is_surf), both vertices and faces supplied and in correct format

    ## Mesh: Create INLA mesh  -----------------------------------------------------
    mesh <- make_mesh(vertices, faces)
    if (mesh$n != nV) { stop("Mesh has ", mesh$n, " locations, but the data has ", nV, " locations.") }

    ## Mask: check or make.  -----------------------------------------------------
    # Get `mask` based on intersection of input mask and `make_mask` checks.
    if (is.null(mask)) { mask <- rep(TRUE, ncol(data[[1]]$BOLD)) } #start with provided mask or mask of all 1's
    masks_quality <- make_mask(data, meanTol=meanTol, varTol=varTol) #, snrTol=snrTol)
    mask <- mask & masks_quality$mask #remove locations with bad data
    if (!any(mask)) { stop("No in-mask data locations.") }

    # If any masked locations, apply to `mesh` and `data`.
    mesh_orig <- NULL #for output only. initialize to NULL, only update if applying a mask to the mesh
    mask <- as.logical(mask) # [TO DO] check that provided mask is logical
    if (!all(mask)) {
      # `mesh`
      mesh_orig <- mesh #for later plotting
      mesh <- submesh(mask, mesh)
      mask <- !is.na(mesh$idx$loc) #update mask (sometimes vertices not excluded by mask will be excluded in mesh)
      mesh$idx$loc <- mesh$idx$loc[mask]
      # `data`
      for (ss in 1:nS) { data[[ss]]$BOLD <- data[[ss]]$BOLD[,mask,drop=FALSE] }
    }
    # Update number of locations after masking
    nV <- sum(mask); nV_all <- length(mask)

    data_loc <- rep(TRUE, nV) #because we are not including boundary vertices in surface case

    spde <- INLA::inla.spde2.matern(mesh)
  }

  if(is_vol){

    #[TO DO] Allow the user to additionally specify a mask input excluding certain within-ROI locations
    #Simply remove those locations (in addition to bad data locations, as below) from the labels array before creating SPDE

    ## SPDE: Create  -----------------------------------------------------

    #labels object here is an array of ROI labels  #[TO DO] Check that labels contains positive integers only, and zeros for background
    mask <- (labels != 0) # mask of the ROI within the full brain volume
    masks_quality <- make_mask(data, meanTol=meanTol, varTol=varTol)
    mask2 <-  masks_quality$mask #mask within the mask (remove locations with low mean or variance)
    mask[mask==TRUE] <- mask2 #remove bad locations from mask
    labels[mask==FALSE] <- 0 #remove bad locations from labels
    ROIs <- unique(labels[mask])
    nR <- length(ROIs)

    #construct the C and G for the SPDE by block-diagonalizing over ROIs
    C_list <- G_list <- spde_list <- vector('list', length=nR)
    for(r in 1:nR){
      mask_r <- (labels == ROIs[r])
      spde_list[[r]] <- vol2spde(mask_r, nbhd_order = nbhd_order, buffer=buffer)
      C_list[[r]] <- spde_list[[r]]$mats$C
      G_list[[r]] <- spde_list[[r]]$mats$G
    }
    C_sub <- Matrix::bdiag(C_list)
    G_sub <- Matrix::bdiag(G_list)

    #construct the SPDE
    Elog.kappa <- Elog.tau <- 0 #prior means for log(kappa) and log(tau)
    Qlog.kappa <- Qlog.tau <- 0.1 #prior precisions for log(kappa) and log(tau)
    spde <- INLA::inla.spde2.generic(M0 = C_sub,
                                   M1 = G_sub,
                                   M2 = G_sub%*%solve(C_sub, G_sub),
                                   theta.mu = c(Elog.kappa, Elog.tau),
                                   theta.Q = diag(c(Qlog.kappa, Qlog.tau)),
                                   B0 = matrix(c(0, 1, 0), 1, 3),
                                   B1 = 2*matrix(c(0, 0, 1), 1, 3),
                                   B2 = 1)

    #[TO DO] test this code for multiple regions, it might break
    #get indices of data locations
    data_loc_list <- lapply(spde_list, function(x) which(x$idx2 %in% x$idx))
    data_loc <- data_loc_list[[1]]
    if(nR > 1){
      before <- 0 #cumulative sum of previous regions
      for(r in 2:nR){
        before <- 0 + nrow(spde_list[[r-1]]$mats$C)
        data_loc_r <- data_loc_list[[r]]
        data_loc_r <- data_loc_r + before
        data_loc <- c(data_loc, data_loc_r)
      }
    }

    # Update number of locations after masking
    nV <- sum(mask2); nV_all <- length(mask2);
    mask <- mask2 #mask within the ROIs

    # Mask BOLD data
    if (!all(mask)) {
      for (ss in 1:nS) { data[[ss]]$BOLD <- data[[ss]]$BOLD[,mask,drop=FALSE] }
    }

  }

  if(verbose==1) cat(paste0('\tNumber of data locations: ',length(data_loc),'\n'))
  if(verbose==1) cat(paste0('\tNumber of data + boundary locations: ',spde$n.spde,'\n'))

  # ------------------------------------------------------------------------------
  # Case 1: Fit a models including prewhitening, inference and spatial Bayesian GLM

  if(is.null(design_multiple)){

    # Identify any missing fields across sessions for bookkeeping -----------------

    valid_cols <- matrix(NA, nrow=nS, ncol=ncol(data[[1]]$design))
    for (ss in 1:nS) {
      cols_ss <- (colSums(abs(data[[ss]]$design)) > 0)
      cols_ss[is.na(cols_ss)] <- FALSE
      valid_cols[ss,] <- cols_ss
    }

    ## Scale, nuisance regress, and/or concatenate session data. -----------------
    #collect data and design matrices
    design <- vector('list', length=nS)
    nK2 <- if (is.null(data[[1]]$nuisance)) { 0 } else { ncol(data[[1]]$nuisance) } #number of nuisance regressors
    for (ss in seq(nS)) {
      # Scale data.
      # (`scale_timeseries` expects VxT data, so transpose before and after.)
      data[[ss]]$BOLD <- t(fMRItools::scale_timeseries(
        t(data[[ss]]$BOLD),
        scale=scale_BOLD,
        transpose=FALSE
      ))

      # Remove any missing fields from design matrix for classical GLM
      cols_ss <- valid_cols[ss,]
      if(!all(cols_ss)) warning(paste0('For session ',ss,', ignoring ',sum(!cols_ss),' design matrix columns of zeros for classical GLM.'))

      # Scale design matrix (ignore columns of zeros)
        data[[ss]]$design[,cols_ss] <- if (scale_design) {
        scale_design_mat(data[[ss]]$design[,cols_ss])
      } else {
        scale(data[[ss]]$design[,cols_ss], scale = FALSE)
      }
      design[[ss]] <- data[[ss]]$design #after scaling but before nuisance regression

      #regress nuisance parameters from BOLD data and design matrix
      if ('nuisance' %in% names(data[[ss]])) {
        nuis_ss <- data[[ss]]$nuisance
        stopifnot((is.matrix(nuis_ss) | is.data.frame(nuis_ss)) && is.numeric(nuis_ss))
        nuis_ss <- scale(nuis_ss, scale=FALSE)
        data[[ss]]$BOLD <- nuisance_regression(data[[ss]]$BOLD, nuis_ss)
        data[[ss]]$design[,cols_ss] <- nuisance_regression(data[[ss]]$design[,cols_ss], nuis_ss) #[TO DO] if design matrix varies spatially, need to adapt this. Design matrix will start as TxKxV and continue in that format after this step.
        data[[ss]]$nuisance <- NULL
      }
    }

    # #concatenate sessions if combine_sessions=TRUE
    # if(combine_sessions){
    #   #concatenate BOLD data across all sessions
    #   data <- list(
    #     session_combined = list(
    #       BOLD = do.call(rbind, lapply(data, function(sess){ sess$BOLD })),
    #       design = do.call(rbind, lapply(data, function(sess){ sess$design }))
    #     )
    #   )
    #
    #   # Update nT, nS, session_names
    #   nT <- nrow(data$session_combined$BOLD)
    #   sess_names_orig <- session_names
    #   session_names <- 'session_combined'
    #   nS <- 1
    #   valid_cols <- matrix(TRUE, nrow=1, ncol=nK)
    # } else {
    #   # [TO DO]: allow different `nT`.
    #   # Is this only problematic when `do_pw`?
    #   if (length(unique(nT)) > 1) {
    #     stop("Not supported yet: different BOLD time durations while `combine_sessions=FALSE`.")
    #   }
    #   nT <- nT[1]
    # }

    # [TO DO] "Always prewhiten" even if we do not want to prewhiten so that the data is in a consistent format. Just skip the actual PW steps.  Or at least call organize_data.

    # Prewhitening. --------------------------------------------------------------
    if (do_pw) {
      if (verbose>0) cat("\tPrewhitening...")
      ## Estimate prewhitening parameters. ---------------------------------------
      AR_coeffs <- array(dim = c(nV,ar_order,nS))
      AR_resid_var <- array(dim = c(nV,nS))
      AR_AIC <- if (aic) { array(dim = c(nV,nS)) } else { NULL }

      #estimate prewhitening parameters for each session
      for (ss in 1:nS) {
        cols_ss <- valid_cols[ss,]
        resids <- nuisance_regression(data[[ss]]$BOLD, data[[ss]]$design[,cols_ss]) #[TO DO] if design matrix varies spatially, need to adapt this
        AR_est <- pw_estimate(resids, ar_order, aic=aic)
        AR_coeffs[,,ss] <- AR_est$phi
        AR_resid_var[,ss] <- AR_est$sigma_sq
        if (aic) { AR_AIC[,ss] <- AR_est$aic }
      }

      #average prewhitening parameters across sessions
      avg_AR <- apply(AR_coeffs, 1:2, mean)
      avg_var <- apply(as.matrix(AR_resid_var), 1, mean)
      if (aic) { max_AIC <- apply(AR_AIC, 1, max) } else { max_AIC <- NULL }

      #smooth prewhitening parameters
      if (ar_smooth > 0) {
        AR_smoothed_list <- pw_smooth(
          vertices=mesh$loc, faces=mesh$graph$tv,
          #mask=mask,
          AR=avg_AR, var=avg_var, FWHM=ar_smooth
        )
        avg_AR <- AR_smoothed_list$AR
        avg_var <- AR_smoothed_list$var
      }

      ## Create the sparse pre-whitening matrix. ---------------------------------
      if (is.null(num.threads) | num.threads < 2) {
        # Initialize the block diagonal covariance matrix
        template_pw <- Matrix::bandSparse(
          n = nT, k = 0:(ar_order + 1), symmetric = TRUE #[TO DO]: Check that nT is correct here for multi-session analysis
        )
        template_pw_list <- rep(list(template_pw), nV)
        for (vv in 1:nV) {
          if(vv %% 100 == 0) if (verbose>0) cat("\tLocation",vv,"of",nV,"\n")
          template_pw_list[[vv]] <- .getSqrtInvCpp(
            AR_coeffs = avg_AR[vv,],
            nTime = nT,  #[TO DO]: Check that nT is correct here for multi-session analysis
            avg_var = avg_var[vv]
          )
        }
      } else {
        if (!requireNamespace("parallel", quietly = TRUE)) {
          stop("Prewhitening in parallel requires the `parallel` package. Please install it.", call. = FALSE)
        }
        max_threads <- max(parallel::detectCores(), 25)
        num_threads <- min(max_threads,num.threads)
        cl <- parallel::makeCluster(num_threads)
        template_pw_list <- parallel::clusterMap(
          cl,
          .getSqrtInvCpp,
          AR_coeffs = split(avg_AR, row(avg_AR)),
          nTime = nT,   #[TO DO]: Check that nT is correct here for multi-session analysis
          avg_var = avg_var,
          SIMPLIFY = FALSE
        )
        parallel::stopCluster(cl)
      }
      if (verbose>0) cat("\tDone!\n")

      #consider using a variant of bdiag_m if this is very slow.  See help(Matrix::bdiag)
      sqrtInv_all <- Matrix::bdiag(template_pw_list)

    } else {

      # Standardize variance (in lieu of PW) -----------------------------------

      if (verbose>0) cat("\tStandardizing residual variance...\n")

      ## Estimate residual variance ---------------------------------------
      AR_resid_var <- array(dim = c(nV,nS))

      #estimate resid var for each session
      for (ss in 1:nS) {
        cols_ss <- valid_cols[ss,]
        resids <- nuisance_regression(data[[ss]]$BOLD, data[[ss]]$design[,cols_ss]) #[TO DO] if design matrix varies spatially, need to adapt this
        resid_var_ss <- matrixStats::colVars(resids)
        AR_resid_var[,ss] <- resid_var_ss
      }

      #average across sessions
      avg_var <- apply(as.matrix(AR_resid_var), 1, mean)
      avg_var <- avg_var/mean(avg_var, na.rm=TRUE)

      #[TO DO] Smoothing of variance estimates

      #set up to pre-multiply the data and design by 1/sqrt(var) for each voxel (happens within organize_data below)
      diag_values <- rep(1/sqrt(avg_var), each = nT)
      sqrtInv_all <- Diagonal(x = diag_values)

    } # END PREWHITENING

    # Classical GLM. -------------------------------------------------------------
    #organize data
    y_all <- c()
    XA_all_list <- NULL
    # Classical GLM
    result_classical <- vector('list', length=nS)
    for (ss in seq(nS)) {
      if (verbose>0) {
        if (nS==1) {
          cat('\tFitting classical GLM.\n')
        } else {
          if (ss==1) { cat('\tFitting classical GLM for session ') }
          if (ss!=nS) { cat(paste0(ss, ", ")) } else { cat(paste0(ss, ".\n")) }
        }
      }

      cols_ss <- valid_cols[ss,] #classical GLM will ignore
      nK_ss <- sum(cols_ss) #in case some fields missing

      # Set up vectorized data and big sparse design matrix.
      # Apply prewhitening, if applicable.
      data_ss <- organize_data(
        data[[ss]]$BOLD, data[[ss]]$design,
        n_mesh = spde$n.spde, inds = data_loc,
        sqrtInv_all = sqrtInv_all
      )

      #setup for Bayesian GLM
      y_all <- c(y_all, data_ss$BOLD)
      #XA_ss <- data_ss$design
      XA_ss <- lapply(data_ss$design, function(x) { return(x %*% data_ss$A) }) #post-multiply each design matrix by A (n_data x n_mesh) for Bayesian GLM
      XA_list <- do.call(cbind, XA_ss) #cbind (expanded) design matrices for each field
      XA_list <- list(XA_list)
      names(XA_list) <- session_names[ss]
      XA_all_list <- c(XA_all_list, XA_list)

      #XA_list <- list(data_ss$bigX)
      #names(XA_list) <- session_names[ss]
      #XA_all_list <- c(XA_all_list, XA_list)

      #setup for classical GLM
      y_ss <- data_ss$BOLD #a vector (grouped by location)
      X_ss <- do.call(cbind, data_ss$design) #cbind non-expanded design matrices for each field
      valid_cols_bigX <- (colSums(abs(X_ss)) > 0) #because X_ss is a big sparse matrix, any missing fields will manifest as a block of zeros
      valid_cols_bigX[is.na(valid_cols_bigX)] <- FALSE
      X_ss <- X_ss[,valid_cols_bigX]

      #perform classical GLM after any prewhitening
      beta_hat_s <- SE_beta_hat_s <- matrix(NA, nV_all, nK)
      XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_ss)))
      if (inherits(XTX_inv, "try-error")) {
        stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
      }
      coef_s <- as.matrix(XTX_inv %*% t(X_ss) %*% y_ss) #a vector of (estimates for location 1, estimates for location 2, ...)
      coef_s_mat <- matrix(coef_s, nrow = nV, ncol = nK_ss) #re-form into a VxK matrix
      beta_hat_s[mask==TRUE,cols_ss] <- coef_s_mat
      resid_s <- t(matrix(y_ss - X_ss %*% coef_s, nrow = nT))

      # ESTIMATE STANDARD ERRORS OF ESTIMATES
      #compute residual SD
      #using length(y_ss)/nV instead of nT here because we want nT for single session case and nT*nS for multi-session case
      DOF_true <- (length(y_ss)/nV) - nK_ss - nK2 - 1
      DOF_false <- (length(y_ss)/nV - 1)
      var_error <- matrixStats::rowVars(resid_s) * DOF_false / DOF_true #correct DOF
      if(do_pw) var_error <- rep(mean(var_error), length(var_error)) #if prewhitening has been done, use same estimate of residual SD everywhere
      sd_error <- sqrt(var_error) #length = nV
      #compute SE of betas
      SE_beta_s <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, each = nK_ss) #blocks of XTX_inv represent fields, so we should repeat each location-specific SD K times
      SE_beta_hat_s[mask==TRUE,cols_ss] <- SE_beta_s

      colnames(beta_hat_s) <- colnames(SE_beta_hat_s) <- field_names

      result_classical[[ss]] <- list(
        estimates = beta_hat_s,
        SE_estimates = SE_beta_hat_s,
        resids = resid_s,
        DOF = DOF_true,
        mask = mask
      )
    }
    names(result_classical) <- session_names

    # FIR Model ------------------------------------------------------------------
    # result_FIR <- vector('list', length=nS)
    # for (ss in seq(nS)) {
    #
    #   #check whether to proceed with FIR modeling
    #   FIR_ss <- data[[ss]]$design_FIR
    #   if(is.null(FIR_ss)) next()
    #   nFIR <- ncol(FIR_ss)
    #   print(paste0('Number of FIR regressors: ', nFIR))
    #   print(paste0('Number of volumes: ', nrow(FIR_ss)))
    #   if(nFIR > nrow(FIR_ss)){
    #     warning('More FIR regressors than volumes. Consider reducing FIR_nsec.')
    #     next()
    #   }
    #
    #   if (verbose>0) cat("\tFitting FIR model.\n")
    #
    #   #set up vectorized data and big sparse design matrix
    #   #if(!do_pw) data_ss <- organize_data(data[[ss]]$BOLD, data[[ss]]$design)
    #   #if(do_pw) data_ss <- data[[ss]] #data has already been "organized" (big sparse design) in prewhitening step above
    #
    #   #y_ss <- matrix(data_ss$BOLD, nrow=nT) #a vector (grouped by location)
    #   y_ss <- data[[ss]]$BOLD #[TO DO] implement prewhitening case (may not be needed if we are not doing inference)
    #   X_ss <- cbind(FIR_ss, 1, data[[ss]]$nuisance) #need the intercept since FIR bases are not centered
    #
    #   #fit model
    #   beta_hat_s <- matrix(NA, nV_all, nFIR)
    #   XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_ss)))
    #   if (inherits(XTX_inv, "try-error")) {
    #     stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
    #   }
    #   coef_s <- as.matrix(XTX_inv %*% t(X_ss) %*% y_ss) #a vector of (estimates for location 1, estimates for location 2, ...)
    #   beta_hat_s[mask==TRUE,] <- t(coef_s[1:nFIR,]) #drop the intercept and nuisance, transpose to V x nFIR
    #
    #   colnames(beta_hat_s) <- colnames(FIR_ss)
    #   result_FIR[[ss]] <- beta_hat_s
    # }
    # names(result_FIR) <- session_names
    #

    # Bayesian GLM. --------------------------------------------------------------
    if (do_Bayesian) {
      #construct betas and repls objects
      if(is_surf) {
        n_mesh <- length(mesh$idx$loc)
        if(!all.equal(mesh$idx$loc, 1:n_mesh)) stop('Developer: Check definition of `spatial` in organize_replicates')
        #data_loc <- 1:n_mesh
      }
      if(is_vol) {
        n_mesh <- spde$n.spde
        #data_loc <- data_loc
      }

      replicates_list <- organize_replicates(n_sess=nS,
                                             field_names=field_names,
                                             n_mesh=n_mesh)
                                             #data_loc=data_loc) #indices of original data locations
      betas <- replicates_list$betas
      repls <- replicates_list$repls
      model_data <- make_data_list(y=y_all, X=XA_all_list,
                                   betas=betas, repls=repls)
      Amat <- model_data$X
      model_data$XA_all_list <- NULL

      ## EM Model. ---------------------------------------------------------------
      #if(do_EM) {
        #stop()
        # if (!requireNamespace("MatrixModels", quietly = TRUE)) {
        #   stop("EM requires the `MatrixModels` package. Please install it.", call. = FALSE)
        # }
        # if (verbose>0) cat('\tEstimating Bayesian model with EM.\n')
        # Psi_k <- spde$Amat
        # Psi <- Matrix::bdiag(rep(list(Psi_k),nK))
        # A <- Matrix::crossprod(model_data$X %*% Psi)
        # # Initial values for kappa and tau
        # kappa2 <- 4
        # phi <- 1 / (4*pi*kappa2*4)
        # # Using values based on the classical GLM
        # if (verbose>0) cat("\t\tFinding best guess initial values.\n")
        # beta_hat <- MatrixModels:::lm.fit.sparse(model_data$X, model_data$y)
        # res_y <- (model_data$y - model_data$X %*% beta_hat)@x
        # sigma2 <- stats::var(res_y)
        # beta_hat <- matrix(beta_hat, ncol = nK*nS)
        # rcpp_spde <- create_listRcpp(spde$spde)
        # if(nS > 1) {
        #   field_cols <- sapply(seq(nS), function(ss) seq(nK) + nK *(ss - 1))
        #   beta_hat <- apply(field_cols,1,function(x) beta_hat[,x])
        # }
        # n_threads <- parallel::detectCores()
        # n_threads <- min(n_threads,nK,num.threads)
        # cl <- parallel::makeCluster(n_threads)
        # kappa2_phi_rcpp <- parallel::parApply(
        #   cl = cl,
        #   beta_hat,
        #   2,
        #   .initialKP,
        #   theta = c(kappa2, phi),
        #   spde = rcpp_spde,
        #   n_sess = nS,
        #   tol = emTol,
        #   verbose = FALSE
        # )
        # parallel::stopCluster(cl)
        # if (verbose>0) cat("\t\tDone!\n")
        # theta <- c(t(kappa2_phi_rcpp), sigma2)
        # theta_init <- theta
        # Ns <- 50 # This is a level of approximation used for the Hutchinson trace estimator
        # if(verbose>0) cat("\t\tStarting EM algorithm.\n")
        # em_output <-
        #   .findTheta(
        #     theta = theta,
        #     spde = rcpp_spde,
        #     y = model_data$y,
        #     X = model_data$X,
        #     QK = make_Q(theta, rcpp_spde, nS),
        #     Psi = as(Psi, "dgCMatrix"),
        #     A = as(A, "dgCMatrix"),
        #     Ns = 50,
        #     tol = emTol,
        #     verbose = verbose>0
        #   )
        # if(verbose>0) cat("\t\tEM algorithm complete!\n")
        # kappa2_new <- phi_new <- sigma2_new <- mu_theta <- NULL
        # list2env(em_output, envir = environment())
        # Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
        #                  MoreArgs = list(spde=rcpp_spde), SIMPLIFY = F)
        # Q_theta <- Matrix::bdiag(Qk_new)
        # if(nS > 1) Q_theta <- Matrix::bdiag(lapply(seq(nS), function(x) Q_theta))
        # Sig_inv <- Q_theta + A/sigma2_new
        # m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
        # mu_theta <- Matrix::solve(Sig_inv, m)
        # # Prepare results
        # field_estimates <- matrix(NA, nrow = length(mask), ncol = nK*nS)
        # field_estimates[mask == 1,] <- matrix(mu_theta,nrow = nV, ncol = nK*nS)
        # colnames(field_estimates) <- rep(field_names, nS)
        # field_estimates <- lapply(seq(nS), function(ss) field_estimates[,(seq(nK) + nK * (ss - 1))])
        # names(field_estimates) <- session_names
        # avg_field_estimates <- NULL
        # if(combine_sessions) avg_field_estimates <- Reduce(`+`,field_estimates) / nS
        # theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
        # names(theta_estimates) <- c("sigma2",paste0("phi_",seq(nK)),paste0("kappa2_",seq(nK)))
        # #extract stuff needed for group analysis
        # tau2_init <- 1 / (4*pi*theta_init[seq(nK)]*theta_init[(seq(nK) + nK)])
        # mu_init <- c(log(1/tail(theta_init,1)), c(rbind(log(sqrt(tau2_init)),log(sqrt(theta_init[seq(nK)])))))
        # tau2 <- 1 / (4*pi*kappa2_new*phi_new)
        # mu_theta <- c(log(1/sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new)))))
        # if (verbose>0) cat("\t\tDone!\n")

      ## INLA Model. -------------------------------------------------------------
      #} else {
        #estimate model using INLA
        if (verbose>0) cat('\tEstimating Bayesian model with INLA...')
        #organize the formula and data objects
        repl_names <- names(repls)
        hyper_initial <- c(-2,2)
        hyper_initial <- rep(list(hyper_initial), nK)
        hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

        formula_vec <- paste0('f(',field_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
        formula_vec <- c('y ~ -1', formula_vec)
        formula_str <- paste(formula_vec, collapse=' + ')
        formula <- as.formula(formula_str)

        INLA_model_obj <- INLA::inla(
          formula,
          data=model_data,
          #data=INLA::inla.stack.data(model_data, spde=spde),
          control.predictor=list(A=Amat, compute = TRUE),
          verbose = verbose>1, keep = FALSE, num.threads = num.threads,
          control.inla = list(strategy = "gaussian", int.strategy = "eb"),
          control.family=list(hyper=list(prec=list(initial=1))),
          control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
        )
        if (verbose>0) cat("\tDone!\n")

        #extract useful stuff from INLA model result
        field_estimates <- extract_estimates(
          INLA_model_obj=INLA_model_obj,
          session_names=session_names,
          mask=mask,
          inds=data_loc
        ) #posterior means of latent field
        theta_estimates <- INLA_model_obj$summary.hyperpar$mean
        hyperpar_posteriors <- get_posterior_densities2(
          INLA_model_obj=INLA_model_obj, #spde,
          field_names
        ) #hyperparameter posterior densities

        #construct object to be returned
        INLA_model_obj <- switch(return_INLA,
          trimmed=trim_INLA_model_obj(INLA_model_obj, minimal=FALSE),
          full=INLA_model_obj,
          minimal=trim_INLA_model_obj(INLA_model_obj, minimal=TRUE)
        )
        attr(INLA_model_obj, "format") <- return_INLA
      #}
    } else {
      field_estimates <- lapply(result_classical, function(x){ x$estimates })
      attr(field_estimates, "GLM_type") <- "classical"
    }

    # Clean up and return. -------------------------------------------------------
    prewhiten_info <- if (do_pw) {
      list(phi = avg_AR, sigma_sq = avg_var, AIC = max_AIC)
    } else {
      NULL
    }

    result_multiple <- NULL

  } else {

    # ------------------------------------------------------------------------------
    # Case 2: Fit and compare multiple models

    y_ss <- data[[1]]$BOLD
    X_ss <- design_multiple[[1]]
    X2_ss <- data[[1]]$nuisance

    # Step 1: Identify no-signal locations. Compare null vs. canonical model using out-of-sample prediction error.

    #ingredients for prediction
    nT2 <- round(nT/2)
    inds1 <- 1:nT2
    inds2 <- setdiff(1:nT, inds1)
    y_list <- list(y_ss[inds1,], y_ss[inds2,])
    X2_list <- list(cbind(1, X2_ss[inds1,]), cbind(1, X2_ss[inds2,])) #for nuisance regression of the "held out" data
    X1_list <- list(design_can[inds1,], design_can[inds2,]) #canonical HRF task regressors

    RSS_OS <- matrix(0, nrow=nV_all, ncol=2)
    RSS_OS[mask==FALSE,] <- NA
    for(pred in 1:2){

      train <- pred
      test <- setdiff(1:2, train)

      # (i) estimate task coefficients based on training set (for canonical model only, not necessary for null model)
      y_train <- y_list[[train]]
      X2_train <- X2_list[[train]]
      X1_train <- X1_list[[train]]
      X_train_can <- cbind(X1_train, X2_train)
      nK_can <- ncol(X1_train)

      XtX_inv <- try(Matrix::solve(Matrix::crossprod(X_train_can))) #this includes nuisance regressors
      if (inherits(XtX_inv, "try-error")) {
        warning(paste0("Numerical instability in design matrix"))
      }
      coefs_can <- (XtX_inv %*% t(X_train_can) %*% y_train)[1:nK_can,] #save task coefficients only (KxV)

      # (ii) do nuisance regression on test set
      y_test <- y_list[[test]]
      X2_test <- X2_list[[test]]
      X1_test <- X1_list[[test]]

      Imat <- diag(1, nrow(X2_test))
      Hmat <- X2_test %*% try(Matrix::solve(Matrix::crossprod(X2_test))) %*% t(X2_test)
      y_test_nuis <- (Imat - Hmat) %*% y_test #regress out nuisance from y
      X1_test_nuis <- (Imat - Hmat) %*% X1_test #regress out nuisance from task regressors

      # (iii) apply coefficients to generate prediction errors
      print(dim(y_test_nuis))
      print(dim(X1_test_nuis))
      print(dim(coefs_can))
      resid_list <- list(canonical = y_test_nuis - X1_test_nuis %*% coefs_can, null = y_test_nuis)
      RSS_OS_pred <- sapply(resid_list, function(x) colSums(x^2)) #Vx2
      RSS_OS[mask==TRUE,] <- RSS_OS[mask==TRUE,] + RSS_OS_pred #sum over both directions of prediction
    }

    noHRF <- (RSS_OS[,2] < RSS_OS[,1]) #for which locations is the null model RSS less than the canonical error RSS

    #loop over models
    nP <- dim(X_ss)[3]
    RSS <- matrix(NA, nrow=nV_all, ncol=nP) #keep track of residual sum of squares (proxy for R^2 or AIC)

    if(verbose > 0) cat('\tFitting models: Model ')
    for(pp in 1:nP){

      if(verbose > 0) cat(paste0(pp,'\t'))

      X_sp <- cbind(X_ss[,,pp], rep(1, nT), X2_ss) #combine design with intercept and nuisance

      #1. Compute residual SD

      XtX_inv_pp <- try(Matrix::solve(Matrix::crossprod(X_sp)))
      if (inherits(XtX_inv_pp, "try-error")) {
        warning(paste0("Numerical instability in design matrix for model ",pp))
      }
      coef_pp <- XtX_inv_pp %*% t(X_sp) %*% y_ss
      resid_pp <- y_ss - X_sp %*% coef_pp #TxV matrix
      RSS[mask==TRUE,pp] <- sqrt(colSums(resid_pp^2)/(nT - ncol(X_sp)))

      #determine best model (minimum residual squared error)
      bestmodel <- apply(RSS, 1, function(x){
        wm <- which.min(x)
        varx <- var(x, na.rm=TRUE)
        if(is.na(varx)) varx <- 0
        if(varx==0) wm <- NA
        wm
      })

      result_multiple <- list(bestmodel = bestmodel,
                                    noHRF = noHRF,
                                    RSS = RSS,
                                    RSS_OS = RSS_OS)
    }

    result_classical <- prewhiten_info <- design <- NULL

  }

  result <- list(
    field_estimates = field_estimates,
    INLA_model_obj = INLA_model_obj,
    hyperpar_posteriors = hyperpar_posteriors,
    theta_estimates = theta_estimates,
    result_classical = result_classical,
    result_multiple = result_multiple,
    #result_FIR = result_FIR,
    mesh = mesh,
    spde = spde,
    data_loc = data_loc,
    mesh_orig = mesh_orig,
    mask = mask,
    mask_orig = mask_orig, #[TO DO] return the params passed into the function instead?
    masks_quality = masks_quality,
    design = design,
    field_names = field_names,
    session_names = session_names,
    n_sess_orig = nS_orig,
    # For joint group model ~~~~~~~~~~~~~
    #posterior_Sig_inv = Sig_inv,
    y = y_all,
    X = XA_all_list,
    prewhiten_info = prewhiten_info,
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call = match.call()
  )
  class(result) <- "BayesGLM"

  result
}


#' Bayes GLM arg checks
#'
#' Checks arguments for \code{BayesGLM} and \code{BayesGLM_cifti}
#'
#' Avoids duplicated code between \code{BayesGLM} and \code{BayesGLM_cifti}
#'
#' @param scale_BOLD,scale_design See \code{\link{BayesGLM}}.
#' @param Bayes,EM See \code{\link{BayesGLM}}.
#' @param ar_order,ar_smooth,aic See \code{\link{BayesGLM}}.
#' @param num.threads See \code{\link{BayesGLM}}.
#' @param return_INLA See \code{\link{BayesGLM}}.
#' @param verbose See \code{\link{BayesGLM}}.
# @param combine_sessions See \code{\link{BayesGLM}}.
#' @param meanTol,varTol,emTol See \code{\link{BayesGLM}}.
#'
#' @return The arguments that may have changed, in a list: \code{scale_BOLD},
#'  \code{do_Bayesian}, \code{do_EM}, and \code{do_pw}.
#'
#' @keywords internal
BayesGLM_argChecks <- function(
    #combine_sessions = FALSE,
    scale_BOLD = c("auto", "mean", "sd", "none"),
    scale_design = TRUE,
    Bayes = TRUE,
    EM = FALSE,
    ar_order = 6,
    ar_smooth = 5,
    aic = FALSE,
    num.threads = 4,
    return_INLA = c("trimmed", "full", "minimal"),
    verbose=1,
    meanTol=1e-6,
    varTol=1e-6,
    emTol=1e-3
){

  #stopifnot(is_1(combine_sessions, "logical"))

  if (isTRUE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'auto'"); scale_BOLD <- "auto"
  }
  if (isFALSE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'none'"); scale_BOLD <- "none"
  }
  scale_BOLD <- match.arg(scale_BOLD, c("auto", "mean", "sd", "none"))
  stopifnot(is_1(scale_design, "logical"))

  stopifnot(is_1(Bayes, "logical"))
  stopifnot(is_1(EM, "logical"))
  if (EM && !Bayes) {
    warning("EM is a Bayesian method: setting `Bayes` to `TRUE`.")
    Bayes <- TRUE
  }
  if (Bayes) {
    if (!EM) { check_INLA(require_PARDISO=FALSE) }
  }

  if (isTRUE(return_INLA)) {
    message("Setting `return_INLA` to 'trimmed'"); return_INLA <- "trimmed"
  }
  if (isFALSE(return_INLA)) {
    message("Setting `return_INLA` to 'minimal'"); return_INLA <- "minimal"
  }
  return_INLA <- match.arg(return_INLA, c("trimmed", "full", "minimal"))

  # Rename these arguments.
  do_Bayesian <- Bayes; rm(Bayes)
  do_EM <- EM; rm(EM)

  if(do_EM) stop() #not currently available

  if (is.null(ar_order)) ar_order <- 0
  stopifnot(is_1(ar_order, "numeric"))
  do_pw <- ar_order > 0
  if (is.null(ar_smooth)) ar_smooth <- 0
  stopifnot(is_1(ar_smooth, "numeric"))
  stopifnot(is_1(aic, "logical"))
  stopifnot(is_1(num.threads, "numeric"))
  stopifnot(num.threads <= parallel::detectCores())
  if (isTRUE(verbose)) { verbose <- 2 }
  if (isFALSE(verbose)) { verbose <- 0 }
  stopifnot(is_posNum(verbose, zero_ok=TRUE))
  stopifnot(is_posNum(meanTol))
  stopifnot(is_posNum(varTol))
  stopifnot(is_posNum(emTol))

  list(
    scale_BOLD=scale_BOLD,
    do_Bayesian=do_Bayesian,
    do_EM = do_EM,
    do_pw = do_pw,
    return_INLA=return_INLA
  )
}
