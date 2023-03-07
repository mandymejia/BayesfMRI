#' Bayes GLM arg checks
#' 
#' Checks arguments for BayesGLM and BayesGLM_cifti
#' 
#' Avoids duplicated code between BayesGLM and BayesGLM_cifti
#' 
#' @param scale_BOLD,scale_design See \code{\link{BayesGLM}}.
#' @param Bayes,EM See \code{\link{BayesGLM}}.
#' @param ar_order,ar_smooth,aic See \code{\link{BayesGLM}}.
#' @param num.threads See \code{\link{BayesGLM}}.
#' @param return_INLA_result See \code{\link{BayesGLM}}.
#' @param outfile See \code{\link{BayesGLM}}.
#' @param verbose See \code{\link{BayesGLM}}.
#' @param avg_sessions See \code{\link{BayesGLM}}.
#' @param meanTol,varTol,emTol See \code{\link{BayesGLM}}.
#' @param trim_INLA See \code{\link{BayesGLM}}.
#'
#' @return The arguments that may have changed, in a list: \code{scale_BOLD},
#'  \code{do_Bayesian}, \code{do_EM}, and \code{do_pw}.
#' 
#' @keywords internal
BayesGLM_argChecks <- function(
  scale_BOLD = c("auto", "mean", "sd", "none"),
  scale_design = TRUE,
  Bayes = TRUE,
  EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  num.threads = 4,
  return_INLA_result = TRUE,
  outfile = NULL,
  verbose=FALSE,
  avg_sessions = FALSE,
  meanTol=1e-6,
  varTol=1e-6,
  emTol=1e-3,
  trim_INLA = TRUE){

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
    if (!EM) { check_INLA(require_PARDISO=TRUE) }
    if (is.null(outfile)) {
      warning('No value supplied for `outfile`, which is required for post-hoc Bayesian group modeling.')
    }
  }
  # Rename these arguments.
  do_Bayesian <- Bayes; rm(Bayes)
  do_EM <- EM; rm(EM)

  if (is.null(ar_order)) ar_order <- 0
  stopifnot(is_1(ar_order, "numeric"))
  do_pw <- ar_order > 0
  if (is.null(ar_smooth)) ar_smooth <- 0
  stopifnot(is_1(ar_smooth, "numeric"))
  stopifnot(is_1(aic, "logical"))
  stopifnot(is_1(num.threads, "numeric"))
  stopifnot(num.threads <= parallel::detectCores())
  stopifnot(is_1(return_INLA_result, "logical"))
  stopifnot(is_1(outfile, "character"))
  stopifnot(is_1(verbose, "logical"))
  stopifnot(is_1(avg_sessions, "logical"))
  stopifnot(is_posNum(meanTol))
  stopifnot(is_posNum(varTol))
  stopifnot(is_posNum(emTol))
  stopifnot(is_1(trim_INLA, "logical"))

  list(
    scale_BOLD=scale_BOLD,
    do_Bayesian=do_Bayesian,
    do_EM = do_EM, 
    do_pw = do_pw
  )
}

#' BayesGLM
#'
#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param data A list of sessions in the \code{"BfMRI.sess"} object format. Each
#'  session is a list with elements "BOLD", "design", and optionally "nuisance".
#'  The name of each element in \code{data} is the name of that session. See 
#'  \code{?is.BfMRI.sess} for details.
#' 
#'  Note that the argument \code{session_names} can be used instead of providing
#'  the session names by the list element names of \code{data}.
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mesh_Param_inla
#' @param mask (Optional) A length \eqn{V} logical vector indicating the
#'  vertices to be included.
#' @inheritParams task_names_Param
#' @inheritParams session_names_Param
#' @inheritParams contrasts_Param
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams Bayes_Param
#' @inheritParams EM_Param
#' @inheritParams ar_order_Param
#' @inheritParams ar_smooth_Param
#' @inheritParams aic_Param
#' @inheritParams num.threads_Param
#' @inheritParams return_INLA_result_Param
#' @inheritParams outfile_Param
#' @inheritParams verbose_Param_inla
#' @inheritParams avg_sessions_Param
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for both.
#' @inheritParams emTol_Param
#' @inheritParams trim_INLA_Param
#'
#' @return A list containing...
#'
#' @importFrom excursions submesh.mesh
#' @importFrom matrixStats colVars
#' @importFrom Matrix bandSparse bdiag crossprod solve
#' @importFrom parallel detectCores makeCluster clusterMap stopCluster
#' @importFrom stats as.formula
#' @importFrom fMRItools is_1
#'
#' @importFrom utils tail
#'
#' @export
BayesGLM <- function(
  data,
  vertices = NULL,
  faces = NULL,
  mesh = NULL,
  mask = NULL,
  # Below arguments shared with `BayesGLM_cifti`
  task_names = NULL,
  session_names = NULL,
  contrasts = NULL,
  scale_BOLD = c("auto", "mean", "sd", "none"),
  scale_design = TRUE,
  Bayes = TRUE,
  EM = FALSE,
  ar_order = 6,
  ar_smooth = 5,
  aic = FALSE,
  num.threads = 4,
  return_INLA_result = TRUE,
  outfile = NULL,
  verbose = FALSE,
  avg_sessions = TRUE,
  meanTol = 1e-6,
  varTol = 1e-6,
  emTol = 1e-3,
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

  # Preliminary steps. ---------------------------------------------------------
  ## Check simple arguments.
  ## These checks are in a separate function because they are shared with
  ## `BayesGLM_cifti`.
  argChecks <- BayesGLM_argChecks(
    scale_BOLD = scale_BOLD,
    scale_design = scale_design,
    Bayes = Bayes,
    EM = EM,
    ar_order = ar_order,
    ar_smooth = ar_smooth,
    aic = aic,
    num.threads = num.threads,
    return_INLA_result = return_INLA_result,
    outfile = outfile,
    verbose = verbose,
    avg_sessions = avg_sessions,
    varTol = varTol,
    meanTol = meanTol,
    emTol = emTol,
    trim_INLA = trim_INLA
  )
  scale_BOLD <- argChecks$scale_BOLD
  do_Bayesian <- argChecks$do_Bayesian
  do_EM <- argChecks$do_EM
  do_pw <- argChecks$do_pw
  need_mesh <- do_Bayesian || (do_pw && ar_smooth > 0)
  
  ## Sessions and data dimensions. ---------------------------------------------
  if (!is.BfMRI.sess(data)) {
    stop("`data` must be a list of sessions, as described in `?is.BfMRI.sess`.")
  }

  if (is.null(session_names)) {
    session_names <- names(data)
  } else {
    if (!is.null(names(data)) && !identical(session_names, names(data))) {
      warning("Using `session_names` rather than `names(data)`.")
      names(data) <- session_names
    }
  }
  n_sess <- length(session_names)
  if (n_sess == 1 && avg_sessions) avg_sessions <- FALSE
  V <- ncol(data[[1]]$BOLD) #number of data locations
  ntime <- vapply(data, function(x){ nrow(x$BOLD) }, 0)
  K <- ncol(data[[1]]$design) #number of tasks

  ## Mesh: check or make. ------------------------------------------------------
  # We need a mesh if doing spatial Bayesian modeling or any AR smoothing.
  if (need_mesh) {
    # Check that only mesh OR vertices+faces supplied
    has_mesh <- !is.null(mesh)
    has_verts_faces <- !is.null(vertices) && !is.null(faces)
    if (!xor(has_mesh, has_verts_faces)) {
      stop('Must supply EITHER mesh OR vertices and faces.')
    }
    if (is.null(mesh)) mesh <- make_mesh(vertices, faces) # This function has been modified to no longer require INLA (2022-03-24)
    if (mesh$n != V) { stop("Mesh has ", mesh$n, " locations, but the data has ", V, " locations.") }
  } else {
    mesh <- NULL
  }

  ## Mask: check or make.  -----------------------------------------------------
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
      # mesh <- excursions::submesh.mesh(mask, mesh) # This is commented out because we now have our own submesh function!
      mesh <- submesh(mask, mesh)
      mask2 <- !is.na(mesh$idx$loc) #update mask (sometimes vertices not excluded by mask will be excluded in mesh)
      mesh$idx$loc <- mesh$idx$loc[mask2]
    }
    # `data`
    for (ss in 1:n_sess) {
      data[[ss]]$BOLD <- data[[ss]]$BOLD[,mask2,drop=FALSE]
    }
  }
  if (do_Bayesian && !do_EM) {spde <- INLA::inla.spde2.matern(mesh)}
  if (do_EM) {spde <- create_spde_surf(mesh)}

  V <- sum(mask2)
  V_all <- length(mask2)

  ## Beta names: check or make. ------------------------------------------------
  if (!is.null(task_names)) {
    if (length(task_names) != K) {
      stop(
        'I detect ', K,
        ' task based on the design matrix, but the length of task_names is ',
        length(task_names), '.  Please fix task_names.'
      )
    }
  } else {
    # Grab beta names from design (if provided)
    task_names <- colnames(data[[1]]$design)
    if (is.null(task_names)) { task_names <- paste0("beta", seq(K)) }
  }

  ## Scale, nuisance regress, and/or concatenate session data. -----------------
  #collect data and design matrices
  design <- vector('list', length=n_sess)
  K2 <- if (is.null(data[[1]]$nuisance)) { 0 } else { ncol(data[[1]]$nuisance) }
  for (ss in seq(n_sess)) {
    #scale data to represent % signal change (or just center if scale=FALSE)
    data[[ss]]$BOLD <- scale_timeseries(t(data[[ss]]$BOLD), scale=scale_BOLD)
    if (scale_design) {
      data[[ss]]$design <- scale_design_mat(data[[ss]]$design)
    } else {
      data[[ss]]$design <- scale(data[[ss]]$design, scale = FALSE)
    }
    design[[ss]] <- data[[ss]]$design #after scaling but before nuisance regression

    #regress nuisance parameters from BOLD data and design matrix
    if ('nuisance' %in% names(data[[ss]])) {
      nuisance_s <- scale(data[[ss]]$nuisance, scale=FALSE)
      data[[ss]]$BOLD <- nuisance_regression(data[[ss]]$BOLD, nuisance_s)
      data[[ss]]$design <- nuisance_regression(data[[ss]]$design, nuisance_s)
      data[[ss]]$nuisance <- NULL
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

  # Prewhitening. --------------------------------------------------------------
  ## Estimate prewhitening parameters. -----------------------------------------

  #compute AR coefficients and average over sessions
  if (do_pw) {
    AR_coeffs <- array(dim = c(V,ar_order,n_sess))
    AR_resid_var <- array(dim = c(V,n_sess))
    AR_AIC <- if (aic) { array(dim = c(V,n_sess)) } else { NULL }

    #estimate prewhitening parameters for each session
    for (ss in 1:n_sess) {
      resids <- nuisance_regression(data[[ss]]$BOLD, data[[ss]]$design)
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
  }

  ## Apply prewhitening. -------------------------------------------------------
  if (do_pw) {
    # Create the sparse pre-whitening matrix
    cat(".... prewhitening... ")
    if (is.null(num.threads) | num.threads < 2) {
      # Initialize the block diagonal covariance matrix
      template_pw <- Matrix::bandSparse(
        n = ntime, k = 0:(ar_order + 1), symmetric = TRUE
      )
      template_pw_list <- rep(list(template_pw), V)
      for(vv in 1:V) {
        if(vv %% 100 == 0) cat("\n Location",vv,"of",V,"")
        # template_pw_list[[vv]] <- prewhiten.v(AR_coeffs = avg_AR[vv,],
        #                                      ntime = ntime,
        #                                      AR_var = avg_var[vv])
        # This is a more efficient prewhitening function written in C++
        template_pw_list[[vv]] <- getSqrtInvCpp(AR_coeffs = avg_AR[vv,],
                                              nTime = ntime,
                                              avg_var = avg_var[vv])
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
          # prewhiten.v,
          getSqrtInvCpp,
          AR_coeffs = split(avg_AR, row(avg_AR)),
          # ntime = ntime,
          nTime = ntime,
          # AR_var = avg_var,
          avg_var = avg_var,
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

  # Classical GLM. -------------------------------------------------------------
  #organize data
  y_all <- c()
  X_all_list <- NULL
  # Classical GLM
  num_GLM <- n_sess
  classical_session_names <- session_names
  result_classical <- vector('list', length=num_GLM)
  for (ss in seq(num_GLM)) {
    #set up vectorized data and big sparse design matrix
    if(!do_pw) data_s <- organize_data(data[[ss]]$BOLD, data[[ss]]$design)
    if(do_pw) data_s <- data[[ss]] #data has already been "organized" (big sparse design) in prewhitening step above

    y_all <- c(y_all, data_s$BOLD)
    X_list <- list(data_s$design)
    names(X_list) <- session_names[ss]
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

    # ESTIMATE STANDARD ERRORS OF ESTIMATES
    #compute residual SD
    #using length(y_reg)/V instead of ntime here because we want ntime for single session case and ntime*n_sess for multi-session case
    DOF_true <- (length(y_reg)/V) - K - K2 - 1
    DOF_false <- (length(y_reg)/V - 1)
    var_error <- matrixStats::rowVars(resid_s) * DOF_false / DOF_true #correct DOF
    if(do_pw) var_error <- rep(mean(var_error), length(var_error)) #if prewhitening has been done, use same estimate of residual SD everywhere
    sd_error <- sqrt(var_error)
    #compute SE of betas
    SE_beta_s <- sqrt(Matrix::diag(XTX_inv)) * rep(sd_error, times = K) #each?
    SE_beta_hat_s[mask2==TRUE,] <- SE_beta_s

    result_classical[[ss]] <- list(
      estimates = beta_hat_s,
      SE_estimates = SE_beta_hat_s,
      resids = resid_s,
      DOF = DOF_true,
      mask = mask2
    )
  }
  names(result_classical) <- classical_session_names

  # Bayesian GLM. --------------------------------------------------------------
  if (do_Bayesian) {
    #construct betas and repls objects
    replicates_list <- organize_replicates(n_sess=n_sess, task_names=task_names, mesh=mesh)
    betas <- replicates_list$betas
    repls <- replicates_list$repls
    model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

    ## EM Model. ---------------------------------------------------------------
    if(do_EM) {
      if (!requireNamespace("MatrixModels", quietly = TRUE)) {
        stop("EM requires the `MatrixModels` package. Please install it.", call. = FALSE)
      }
      mesh$Amat <- spde$Amat
      mesh$spde <- spde$spde
      cat('\n.... estimating model with EM')
      Psi_k <- spde$Amat
      Psi <- Matrix::bdiag(rep(list(Psi_k),K))
      A <- Matrix::crossprod(model_data$X %*% Psi)
      # Initial values for kappa and tau
      kappa2 <- 4
      phi <- 1 / (4*pi*kappa2*4)
      # Using values based on the classical GLM
      if(verbose) cat("... FINDING BEST GUESS INITIAL VALUES\n")
      beta_hat <- MatrixModels:::lm.fit.sparse(model_data$X, model_data$y)
      res_y <- (model_data$y - model_data$X %*% beta_hat)@x
      sigma2 <- stats::var(res_y)
      beta_hat <- matrix(beta_hat, ncol = K*n_sess)
      rcpp_spde <- create_listRcpp(spde$spde)
      if(n_sess > 1) {
        task_cols <- sapply(seq(n_sess), function(j) seq(K) + K *(j - 1))
        beta_hat <- apply(task_cols,1,function(x) beta_hat[,x])
      }
      n_threads <- parallel::detectCores()
      n_threads <- min(n_threads,K,num.threads)
      cl <- parallel::makeCluster(n_threads)
      kappa2_phi_rcpp <- parallel::parApply(
        cl = cl,
        beta_hat,
        2,
        initialKP,
        theta = c(kappa2, phi),
        spde = rcpp_spde,
        n_sess = n_sess,
        tol = emTol,
        verbose = FALSE
      )
      parallel::stopCluster(cl)
      if(verbose) cat("...... DONE!\n")
      theta <- c(t(kappa2_phi_rcpp), sigma2)
      theta_init <- theta
      Ns <- 50 # This is a level of approximation used for the Hutchinson trace estimator
      if(verbose) cat("... STARTING EM ALGORITHM\n")
      em_output <-
        findTheta(
          theta = theta,
          spde = rcpp_spde,
          y = model_data$y,
          X = model_data$X,
          QK = make_Q(theta, rcpp_spde, n_sess),
          Psi = as(Psi, "dgCMatrix"),
          A = as(A, "dgCMatrix"),
          Ns = 50,
          tol = emTol,
          verbose = verbose
        )
      if(verbose) cat(".... EM algorithm complete!\n")
      kappa2_new <- phi_new <- sigma2_new <- mu <- NULL
      list2env(em_output, envir = environment())
      Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
                       MoreArgs = list(spde=rcpp_spde), SIMPLIFY = F)
      Q <- Matrix::bdiag(Qk_new)
      if(n_sess > 1) Q <- Matrix::bdiag(lapply(seq(n_sess), function(x) Q))
      Sig_inv <- Q + A/sigma2_new
      m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
      mu <- Matrix::solve(Sig_inv, m)
      # Prepare results
      task_estimates <- matrix(NA, nrow = length(mask2), ncol = K*n_sess)
      task_estimates[mask2 == 1,] <- matrix(mu,nrow = V, ncol = K*n_sess)
      colnames(task_estimates) <- rep(task_names, n_sess)
      task_estimates <- lapply(seq(n_sess), function(ns) task_estimates[,(seq(K) + K * (ns - 1))])
      names(task_estimates) <- session_names
      avg_task_estimates <- NULL
      if(avg_sessions) avg_task_estimates <- Reduce(`+`,task_estimates) / n_sess
      theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
      names(theta_estimates) <- c("sigma2",paste0("phi_",seq(K)),paste0("kappa2_",seq(K)))
      #extract stuff needed for group analysis
      tau2_init <- 1 / (4*pi*theta_init[seq(K)]*theta_init[(seq(K) + K)])
      mu.theta_init <- c(log(1/tail(theta_init,1)), c(rbind(log(sqrt(tau2_init)),log(sqrt(theta_init[seq(K)])))))
      tau2 <- 1 / (4*pi*kappa2_new*phi_new)
      mu.theta <- c(log(1/sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new)))))
      cat("... done!\n")

    ## INLA Model. -------------------------------------------------------------
    } else {
      #estimate model using INLA
      cat('\n .... estimating model with INLA')
      #organize the formula and data objects
      repl_names <- names(repls)
      hyper_initial <- c(-2,2)
      hyper_initial <- rep(list(hyper_initial), K)
      hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

      formula_vec <- paste0('f(',task_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
      formula_vec <- c('y ~ -1', formula_vec)
      formula_str <- paste(formula_vec, collapse=' + ')
      formula <- as.formula(formula_str, env = globalenv())

      INLA_result <- INLA::inla(
        formula, data=model_data, control.predictor=list(A=model_data$X, compute = TRUE),
        verbose = verbose, keep = FALSE, num.threads = num.threads,
        control.inla = list(strategy = "gaussian", int.strategy = "eb"),
        control.family=list(hyper=list(prec=list(initial=1))),
        control.compute=list(config=TRUE), contrasts = NULL, lincomb = NULL #required for excursions
      )
      if(verbose) cat("done!\n")

      #extract useful stuff from INLA model result
      task_estimates <- extract_estimates(object=INLA_result, session_names=session_names, mask=mask2) #posterior means of latent task field
      hyperpar_posteriors <- get_posterior_densities(object=INLA_result, spde, task_names) #hyperparameter posterior densities

      #extract stuff needed for group analysis
      mu.theta <- INLA_result$misc$theta.mode
      Q.theta <- solve(INLA_result$misc$cov.intern) #not needed for EM version

      #construct object to be returned
      if(!return_INLA_result){
        INLA_result <- NULL
      } else {
        if(trim_INLA) INLA_result <- trim_INLA_result(INLA_result)
      }
    }
  }

  # Clean up and return. -------------------------------------------------------
  prewhiten_info <- NULL
  if (do_pw) prewhiten_info <- list(phi = avg_AR, sigma_sq = avg_var, AIC=max_AIC)

  if (do_Bayesian) {
    if (do_EM) {
      INLA_result <- hyperpar_posteriors <- Q.theta <- NULL
    } else {
      theta_estimates <- Sig_inv <- NULL
    }
  } else {
    INLA_result <- hyperpar_posteriors <- Q.theta <- NULL
    task_estimates <- hyperpar_posteriors <- mu.theta <- y_all <- X_all_list <- NULL
  }

  result <- list(
    INLA_result = INLA_result,
    task_estimates = task_estimates,
    result_classical = result_classical,
    mesh = mesh,
    mesh_orig = mesh_orig,
    mask = mask,
    design = design,
    task_names = task_names,
    session_names = session_names,
    hyperpar_posteriors = hyperpar_posteriors,
    theta_estimates = theta_estimates,
    # For joint group model ~~~~~~~~~~~~~
    posterior_Sig_inv = Sig_inv, 
    mu.theta = mu.theta, 
    Q.theta = Q.theta, 
    y = y_all, 
    X = X_all_list, 
    prewhiten_info = prewhiten_info,
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call = match.call()
  )

  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    saveRDS(result, file=outfile)
  }

  result
}
