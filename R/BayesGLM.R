#' BayesGLM
#'
#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance. See \code{?create.session} and \code{?is.session}
#'  for more details.
#' List element names represent session names.
#' @param beta_names (Optional) Names of tasks represented in design matrix
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mesh_Param_inla
#' @param mask (Optional) A length \eqn{V} logical vector indicating if each
#'  vertex is to be included.
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams Bayes_Param
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#' @param ar_smooth FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#' @param aic Use the AIC to select AR model order between \code{0} and \code{ar_order}? Default: \code{FALSE}.
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_result_Param_TRUE
#' @param outfile File name where results will be written (for use by
#'  \code{BayesGLM2}).
#' @inheritParams verbose_Param_inla
#' @inheritParams contrasts_Param_inla
#' @inheritParams avg_sessions_Param
#' @param meanTol,varTol Tolerance for mean and variance of each data location. Locations which
#'  do not meet these thresholds are masked out of the analysis. Defaults: \code{1e-6}.
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use `id_activations()`? Default: `TRUE`.
#'
#' @return A list containing...
#'
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @importFrom Matrix bandSparse bdiag crossprod solve
#' @importFrom parallel detectCores makeCluster clusterMap stopCluster
#' @importFrom stats as.formula
#'
#' @export
BayesGLM <- function(
  data,
  beta_names = NULL,
  vertices = NULL,
  faces = NULL,
  mesh = NULL,
  mask = NULL,
  scale_BOLD=c("auto", "mean", "sd", "none"),
  scale_design = TRUE,
  Bayes=TRUE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  num.threads=4,
  return_INLA_result=TRUE,
  outfile = NULL,
  verbose=FALSE,
  contrasts = NULL,
  avg_sessions = FALSE,
  meanTol=1e-6,
  varTol=1e-6,
  trim_INLA = TRUE){

  #TO DO:
  #add "(ignored if classicalGLM_only = TRUE) to some params"
  #add if statements for some of code if classicalGLM_only = TRUE
  #make need_mesh object for if classicalGLM_only == TRUE and no AR smoothing
  #do masking only (no involvement of mesh) if need_mesh == FALSE <-- copy from classicalGLM()

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  # Rename and coerce to logical
  do_Bayesian <- as.logical(Bayes)
  if (do_Bayesian) { check_INLA(require_PARDISO=TRUE) }

  # Check prewhitening arguments.
  if (is.null(ar_order)) ar_order <- 0
  ar_order <- as.numeric(ar_order)
  do_pw <- (ar_order > 0)
  if (is.null(ar_smooth)) ar_smooth <- 0
  ar_smooth <- as.numeric(ar_smooth)

  #we need a mesh if doing spatial Bayesian modeling OR any AR smoothing
  need_mesh <- (do_Bayesian || (do_pw && ar_smooth > 0))
  #check that only mesh OR vertices+faces supplied
  if (need_mesh) {
    has_mesh <- !is.null(mesh)
    has_verts_faces <- !is.null(vertices) && !is.null(faces)
    if (!xor(has_mesh, has_verts_faces)) {
      stop('Must supply EITHER mesh OR vertices and faces.')
    }
    if (is.null(mesh)) mesh <- make_mesh(vertices, faces)
  } else {
    mesh <- NULL
  }

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)
  if (n_sess == 1 && avg_sessions) avg_sessions <- FALSE

  if (!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if (!all.equal(unique(data_classes),'list')) {
    stop('I expect data to be a list of lists (sessions), but it is not')
  }

  #check dimensions
  V <- ncol(data[[1]]$BOLD) #number of data locations
  ntime <- vapply(data, function(x){ nrow(x$BOLD) }, 0)
  K <- ncol(data[[1]]$design) #number of tasks
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }
  if (need_mesh) {
    if (mesh$n != V) { stop("Mesh has ", mesh$n, " locations, but the data has ", V, " locations.") }
  }

  if (!is.null(beta_names)) {
    if(length(beta_names) != K) {
      stop(paste0(
        'I detect ', K,
        ' task based on the design matrix, but the length of beta_names is ',
        length(beta_names), '.  Please fix beta_names.'
      ))
    }
  } else {
    # Grab beta names from design (if provided)
    bn_maybe <- colnames(data[[1]]$design)
    beta_names <- if (is.null(bn_maybe)) { paste0("beta", seq(K)) } else { bn_maybe }
  }

  if (do_Bayesian && is.null(outfile)) {
    message('No value supplied for `outfile`, which is required for post-hoc Bayesian group modeling.')
  }

  # Scale
  if (isTRUE(scale_BOLD)) { cat("Setting `scale_BOLD <- 'auto'`"); scale_BOLD <- "auto" }
  if (isFALSE(scale_BOLD)) { cat("Setting `scale_BOLD <- 'none'`"); scale_BOLD <- "none" }
  scale_BOLD <- match.arg(scale_BOLD, c("auto", "mean", "sd", "none"))

  # `mask`  -----------------------------
  # Get `mask` based on intersection of input mask and `make_mask` checks.
  if (is.null(mask)) { mask <- rep(TRUE, ncol(data[[1]]$BOLD)) }
  mask <- mask & make_mask(data, meanTol=meanTol, varTol=varTol)
  if (!any(mask)) { stop("No in-mask data locations.") }

  # If any masked locations, apply to `mesh` and `data`.
  mesh_orig <- NULL #for output only. initialize to NULL, only update if applying a mask to the mesh
  if (!all(mask)) {
    # `mask2` is in case `need_mesh==FALSE`
    mask <- mask2 <- as.logical(mask)

    # `mesh`
    if (need_mesh) {
      mesh_orig <- mesh #for later plotting
      mesh <- excursions::submesh.mesh(mask, mesh)
      mask2 <- !is.na(mesh$idx$loc) #update mask (sometimes vertices not excluded by mask will be excluded in mesh)
      mesh$idx$loc <- mesh$idx$loc[mask2]
    }
    # `data`
    for (ss in 1:n_sess) {
      data[[ss]]$BOLD <- data[[ss]]$BOLD[,mask2,drop=FALSE]
    }
  }
  if (do_Bayesian) {spde <- INLA::inla.spde2.matern(mesh)}

  V <- sum(mask2)
  V_all <- length(mask2)

  # ---------------------------------
  #collect data and design matrices
  design <- vector('list', length=n_sess)
  for(s in 1:n_sess){
    #scale data to represent % signal change (or just center if scale=FALSE)
    data[[s]]$BOLD <- scale_timeseries(t(data[[s]]$BOLD), scale=scale_BOLD)
    if (scale_design) {
      data[[s]]$design <- scale_design_mat(data[[s]]$design)
    } else {
      data[[s]]$design <- scale(data[[s]]$design, scale = FALSE)
    }
    design[[s]] <- data[[s]]$design #after scaling but before nuisance regression

    #regress nuisance parameters from BOLD data and design matrix
    if ('nuisance' %in% names(data[[s]])) {
      nuisance_s <- scale(data[[s]]$nuisance, scale=FALSE)
      data[[s]]$BOLD <- nuisance_regression(data[[s]]$BOLD, nuisance_s)
      data[[s]]$design <- nuisance_regression(data[[s]]$design, nuisance_s)
      data[[s]]$nuisance <- NULL
    }
  }

  #concatenate sessions if avg_sessions=TRUE
  if(avg_sessions){
    #concatenate BOLD data across all sessions
    data <- list(
      session_avg = list(
        BOLD = do.call(rbind, lapply(data, function(sess){ sess$BOLD })),
        design = do.call(rbind, lapply(data, function(sess){ sess$design }))
      )
    )

    #update ntime, n_sess, session_names
    ntime <- nrow(data$session_avg$BOLD)
    session_names <- 'session_avg'
    n_sess_orig <- n_sess
    n_sess <- 1
  } else {
    # [TO DO]: allow different `ntime`.
    # Is this only problematic when `do_pw`?
    if (length(unique(ntime)) > 1) {
      stop("Not supported yet: different BOLD time durations while `avg_sessions=FALSE`.")
    }
    ntime <- ntime[1]
  }

  ## ESTIMATE PREWHITENING PARAMETERS -----

  #compute AR coefficients and average over sessions
  if (do_pw) {
    AR_coeffs <- array(dim = c(V,ar_order,n_sess))
    AR_resid_var <- array(dim = c(V,n_sess))

    #estimate prewhitening parameters for each session
    for (ss in 1:n_sess) {
      resids <- nuisance_regression(data[[ss]]$BOLD, data[[ss]]$design)
      AR_est <- pw_estimate(resids, ar_order, aic=aic)
      AR_coeffs[,,ss] <- AR_est$phi
      AR_resid_var[,ss] <- AR_est$sigma_sq
    }

    #average prewhitening parameters across sessions
    avg_AR <- apply(AR_coeffs, 1:2, mean)
    avg_var <- apply(as.matrix(AR_resid_var), 1, mean)

    #smooth prewhitening parameters
    if (ar_smooth > 0) {
      AR_smoothed_list <- pw_smooth(
        vertices=mesh$loc, faces=mesh$graph$tv,
        AR=avg_AR, var=avg_var, FWHM=ar_smooth
      )
      avg_AR <- AR_smoothed_list$AR
      avg_var <- AR_smoothed_list$var
    }
  }

  ## APPLY PREWHITENING -----
  if (do_pw) {
    # Create the sparse pre-whitening matrix
    cat(" .... prewhitening... ")
    if (is.null(num.threads) | num.threads < 2) {
      # Initialize the block diagonal covariance matrix
      template_pw <- Matrix::bandSparse(
        n = ntime, k = 0:(ar_order + 1), ymmetric = TRUE
      )
      template_pw_list <- rep(list(template_pw), V)
      for(vv in 1:V) {
        if(vv %% 100 == 0) cat("\n Location",vv,"of",V,"")
        template_pw_list[[vv]] <- prewhiten.v(AR_coeffs = avg_AR[vv,],
                                             ntime = ntime,
                                             AR_var = avg_var[vv])
      }
    } else {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Prewhitening in parallel requires the `parallel` package. Please install it.", call. = FALSE)
      }
      max_threads <- max(parallel::detectCores(), 25)
      num_threads <- min(max_threads,num.threads)
      cl <- parallel::makeCluster(num_threads)
      template_pw_list <-
        parallel::clusterMap(
          cl,
          prewhiten.v,
          AR_coeffs = split(avg_AR, row(avg_AR)),
          ntime = ntime,
          AR_var = avg_var,
          SIMPLIFY = FALSE
        )
      parallel::stopCluster(cl)
    }
    cat("done!\n")

    #consider using a variant of bdiag_m if this is very slow.  See help(Matrix::bdiag)
    sqrtInv_all <- Matrix::bdiag(template_pw_list)

    #apply prewhitening matrix to BOLD and design for each session
    data <- sapply(data, function(data_s) {
      #bold_out <- matrix(NA, ntime, V)
      bold_out <- as.vector(sqrtInv_all %*% c(data_s$BOLD))
      #bold_out[,mask] <- pw_BOLD
      all_design <- organize_data(data_s$BOLD, data_s$design)$design
      pw_design <- sqrtInv_all %*% all_design
      return(list(BOLD = bold_out, design = pw_design))
    }, simplify = F)
  }

  ## ESTIMATE CLASSICAL GLM -----

  #organize data
  y_all <- c()
  X_all_list <- NULL
  # Classical GLM
  num_GLM <- n_sess
  classical_session_names <- session_names
  result_classical <- vector('list', length=num_GLM)
  for(s in 1:num_GLM){
    #set up vectorized data and big sparse design matrix
    if(!do_pw) data_s <- organize_data(data[[s]]$BOLD, data[[s]]$design)
    if(do_pw) data_s <- data[[s]] #data has already been "organized" (big sparse design) in prewhitening step above

    y_all <- c(y_all, data_s$BOLD)
    X_list <- list(data_s$design)
    names(X_list) <- session_names[s]
    X_all_list <- c(X_all_list, X_list)
    y_reg <- data_s$BOLD #a vector (grouped by location)
    X_reg <- data_s$design

    #perform classical GLM after any prewhitening
    beta_hat_s <- SE_beta_hat_s <- matrix(NA, V_all, K)
    XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_reg)))
    if (inherits(XTX_inv, "try-error")) {
      stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
    }
    coef_s <- as.matrix(XTX_inv %*% t(X_reg) %*% y_reg) #a vector of (estimates for location 1, estimates for location 2, ...)
    coef_s_mat <- matrix(coef_s, nrow = V, ncol = K)
    beta_hat_s[mask2==TRUE,] <- coef_s_mat
    resid_s <- t(matrix(y_reg - X_reg %*% coef_s, nrow = ntime))

    # ESTIMATE STANDARD ERRORS OF ESTIIMATES
    #compute residual SD
    #using length(y_reg)/V instead of ntime here because we want ntime for single session case and ntime*n_sess for multi-session case
    DOF_true <- (length(y_reg)/V) - K - 1
    DOF_false <- (length(y_reg)/V - 1)
    var_error <- matrixStats::rowVars(resid_s) * DOF_false / DOF_true #correct DOF
    if(do_pw) var_error <- rep(mean(var_error), length(var_error)) #if prewhitening has been done, use same estimate of residual SD everywhere
    sd_error <- sqrt(var_error)
    #compute SE of betas
    SE_beta_s <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, times = K) #each?
    SE_beta_hat_s[mask2==TRUE,] <- SE_beta_s

    result_classical[[s]] <- list(estimates = beta_hat_s,
                                  SE_estimates = SE_beta_hat_s,
                                  resids = resid_s,
                                  DOF = DOF_true,
                                  mask = mask2)
  }
  names(result_classical) <- classical_session_names

  ## ESTIMATE BAYESIAN GLM -----

  if(do_Bayesian){

    ### THIS PART WILL CHANGE WITH EM VERSION

    #construct betas and repls objects
    replicates_list <- organize_replicates(n_sess=n_sess, beta_names=beta_names, mesh=mesh)
    betas <- replicates_list$betas
    repls <- replicates_list$repls

    #organize the formula and data objects
    repl_names <- names(repls)
    hyper_initial <- c(-2,2)
    hyper_initial <- rep(list(hyper_initial), K)
    hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

    formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
    formula_vec <- c('y ~ -1', formula_vec)
    formula_str <- paste(formula_vec, collapse=' + ')
    formula <- as.formula(formula_str, env = globalenv())

    model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

    #estimate model using INLA
    cat('\n .... estimating model with INLA')
    check_INLA(require_PARDISO=FALSE)

    INLA_result <- INLA::inla(
      formula, data=model_data, control.predictor=list(A=model_data$X, compute = TRUE),
      verbose = verbose, keep = FALSE, num.threads = num.threads,
      control.inla = list(strategy = "gaussian", int.strategy = "eb"),
      control.family=list(hyper=list(prec=list(initial=1))),
      control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
    )
    cat("done!\n")

    #extract useful stuff from INLA model result
    beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names, mask=mask2) #posterior means of latent task field
    theta_posteriors <- get_posterior_densities(object=INLA_result, spde, beta_names) #hyperparameter posterior densities

    #extract stuff needed for group analysis
    mu.theta <- INLA_result$misc$theta.mode
    Q.theta <- solve(INLA_result$misc$cov.intern) #not needed for EM version

    #construct object to be returned
    if(!return_INLA_result){
      INLA_result <- NULL
    } else {
      if(trim_INLA) INLA_result <- trim_INLA_result(INLA_result)
    }

    ### END OF PART THAT WILL CHANGE WITH EM VERSION

  }

  if(do_pw) prewhiten_info <- list(phi = avg_AR, sigma_sq = avg_var)
  if(!do_pw) prewhiten_info <- NULL
  if(!do_Bayesian) INLA_result <- beta_estimates <- theta_posteriors <- mu.theta <- Q.theta <- y_all <- X_all_list <- NULL

  result <- list(INLA_result = INLA_result,
                 beta_estimates = beta_estimates,
                 result_classical = result_classical,
                 mesh = mesh,
                 mesh_orig = mesh_orig,
                 mask = mask,
                 design = design,
                 session_names = session_names,
                 beta_names = beta_names,
                 theta_posteriors = theta_posteriors,
                 mu.theta = mu.theta, #for joint group model
                 Q.theta = Q.theta, #for joint group model
                 y = y_all, #for joint group model
                 X = X_all_list, #for joint group model
                 prewhiten_info = prewhiten_info,
                 call = match.call())

  class(result) <- "BayesGLM"


  if(!is.null(outfile)){
    saveRDS(result, file=outfile)
  }

  return(result)
}
