#' Group-level Bayesian GLM
#'
#' Performs group-level Bayesian GLM estimation and inference using the joint
#'  approach described in Mejia et al. (2020).
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param results Either (1) a length \eqn{N} list of \code{"BayesGLM"} objects,
#'  or (2) a length \eqn{N} character vector of files storing \code{"BayesGLM"}
#'  objects saved with \code{\link{saveRDS}}.
#' @param contrasts (Optional) A list of contrast vectors that specify the
#'  group-level summaries of interest. If \code{NULL}, use contrasts that compute
#'  the average of each field (task HRF) across subjects and sessions.
#'
#'  Each contrast vector is length \eqn{K * S * N} vector specifying a
#'  group-level summary of interest, where \eqn{K} is the number
#'  of fields (task HRFs), \eqn{S} is the number of sessions, and \eqn{N} is the
#'  number of subjects. For a single subject-session the contrast
#'  for the first field would be:
#'
#'  \code{contrast1 <- c(1, rep(0, K-1))}
#'
#'  and so the full contrast vector representing the group average across
#'  sessions and subjects for the first task would be:
#'
#'  \code{rep(rep(contrast1, S), N) /S /N}.
#'
#'  To obtain the group average for the first task, for just the first sessions
#'  from each subject:
#'
#'  \code{rep(c(contrast1, rep(0, K*(S-1))), N) /N}.
#'
#'  To obtain the mean difference between the first and second sessions, for the
#'  first task:
#'
#'  \code{rep(c(contrast1, -contrast1, rep(0, K-2)), N) /N}.
#'
#'  To obtain the mean across sessions of the first task, just for the first
#'  subject:
#'
#'  \code{c(rep(contrast1, S-1), rep(0, K*S*(N-1)) /S}.
#'
#' @param quantiles (Optional) Vector of posterior quantiles to return in
#'  addition to the posterior mean.
#' @param excursion_type (For inference only) The type of excursion function for
#'  the contrast (">", "<", "!="), or a vector thereof (each element
#'  corresponding to one contrast).  If \code{NULL}, no inference performed.
#' @param contrast_names (Optional) Names of contrasts.
#' @param gamma (For inference only) Activation threshold for the excursion set,
#'  or a vector thereof (each element corresponding to one contrast). Default:
#'  \code{0}.
#' @param alpha (For inference only) Significance level for activation for the
#'  excursion set, or a vector thereof (each element corresponding to one
#'  contrast). Default: \code{.05}.
#' @param nsamp_theta Number of theta values to sample from posterior. Default:
#'  \code{50}.
#' @param nsamp_beta Number of beta vectors to sample conditional on each theta
#'  value sampled. Default: \code{100}.
#' @param num_cores The number of cores to use for sampling betas in parallel. If
#'  \code{NULL} (default), do not run in parallel.
#' @inheritParams verbose_Param
#'
#' @return A list containing the estimates, PPMs and areas of activation for each contrast.
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag crossprod
#'
#' @export
BayesGLM2 <- function(
  results,
  contrasts = NULL,
  quantiles = NULL,
  excursion_type=NULL,
  contrast_names = NULL,
  gamma = 0,
  alpha = 0.05,
  nsamp_theta = 50,
  nsamp_beta = 100,
  num_cores = NULL,
  verbose = 1){

  use_INLA <- TRUE # alternative: use EM model, but it's been removed.

  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("`BayesGLM2` requires the `abind` package. Please install it.", call. = FALSE)
  }

  # Check `results`, reading in the files if needed.
  results_ok <- FALSE
  if (is.character(results)) {
    if (!all(file.exists(results))) {
      stop("`results` is a character vector, but not all elements are existing files.")
    }
    results <- lapply(results, readRDS) # [TO DO]: delete w/ each read-in, stuff not needed
  }
  if (!is.list(results)) {
    stop("`results` must be a list of `'BayesGLM'` objects, or a character vector of files with `'BayesGLM'` results saved.")
  }
  is_BayesGLM <- all(vapply(results, inherits, FALSE, "BayesGLM"))
  is_cifti <- all(vapply(results, inherits, FALSE, "BayesGLM_cifti"))
  if (!is_BayesGLM && !is_cifti) {
    stop("`results` must be a list of all `'BayesGLM'` or all `'BayesGLM_cifti'` objects, or a character vector of files with `'BayesGLM(_cifti)'` results.")
  }
  rm(is_BayesGLM) # use `is_cifti`

  model_names <- if (is_cifti) {
    names(results[[1]]$BayesGLM_results)[!vapply(results[[1]]$BayesGLM_results, is.null, FALSE)]
  } else {
    "BayesGLM"
  }

  nM <- length(model_names)                 # models
  nN <- length(results)                     # subjects
  nS <- length(results[[1]]$session_names)  # sessions
  nK <- length(results[[1]]$task_names)     # fields

  session_names <- results[[1]]$session_names
  task_names <- results[[1]]$task_names

  # Check that every subject has the same models, sessions, tasks
  for (nn in seq(nN)) {
    sub_nn <- results[[nn]]
    stopifnot(identical(
      model_names,
      names(sub_nn$BayesGLM_results)[!vapply(sub_nn$BayesGLM_results, is.null, FALSE)]
    ))
    stopifnot(identical(session_names, sub_nn$session_names))
    stopifnot(identical(task_names, sub_nn$task_names))
  }

  # Check `contrasts`.
  # `contrasts` should be fields * sessions * subjects
  if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
  if(is.null(contrasts)) {
    if (verbose>0) cat('Using a contrast that computes the average across subjects for each task. If other contrasts are desired, provide `contrasts`.\n')
    contrasts <- vector('list', length=nK)
    names(contrasts) <- paste0(task_names, '_avg')
    for (kk in 1:nK) {
      # (1/J, 0, 0, ..., 0) for k=1,
      # (0, 1/J, 0, ..., 0) for k=2,
      # ...,
      # (0, 0, ..., 0, 1/J) for k=K
      # for each session, for each subject
      # where J == S * N
      contrast_1 <- c(rep(0, kk-1), 1/(nS*nN), rep(0, nK-kk)) # length nK
      contrasts[[kk]] <- rep(rep(contrast_1, nS), nN)         # length nK*nS*nN
    }
  } else {
    #Check that each contrast vector is numeric and length J*K
    if(any(sapply(contrasts, length) != nK*nS*nN)) {
      stop('Each contrast vector must be of length K*S*N (fields times sessions times subjects).')
    }
    if(any(!sapply(contrasts, is.numeric))) {
      stop('Each contrast vector must be numeric, but at least one is not.')
    }
    if (is.null(names(contrasts))) {
      names(contrasts) <- paste0("contrast_", seq(length(contrasts)))
    }
  }
  # Override `names(contrasts)` with `contrast_names` if provided.
  if (!is.null(contrast_names)) {
    stopifnot(length(contrast_names) == length(contrasts))
    names(contrasts) <- contrast_names
  }
  nC <- length(contrasts)

  #Check `quantiles`
  if(!is.null(quantiles)){
    stopifnot(is.numeric(quantiles))
    if(any(quantiles > 1 | quantiles < 0)) stop('All elements of `quantiles` must be between 0 and 1.')
  }

  do_excur <- !is.null(excursion_type) && (!identical(excursion_type, "none"))
  if (do_excur) {
    if(length(excursion_type) == 1) excursion_type <- rep(excursion_type, nC)
    if(length(gamma) == 1) gamma <- rep(gamma, nC)
    if(length(alpha) == 1) alpha <- rep(alpha, nC)
    if(length(gamma) != nC) stop('Length of gamma must match number of contrasts or be equal to one.')
    if(length(alpha) != nC) stop('Length of alpha must match number of contrasts or be equal to one.')
    if(length(excursion_type) != nC) stop('Length of excursion_type must match number of contrasts or be equal to one.')
  } else {
    excursion_type <- 'none'
  }

  out <- vector("list", nM)
  names(out) <- model_names

  for (mm in seq(nM)) {

    if (nM>1) { if (verbose>0) cat(model_names[mm], " ~~~~~~~~~~~\n") }
    results_mm <- lapply(results, function(x){ x$BayesGLM_results[[mm]] })

    # `Mask`
    Mask <- lapply(results_mm, function(x){ x$mask })
    if (length(unique(vapply(Mask, length, 0))) != 1) { 
      stop("Unequal mask lengths--check that the input files are in the same resolution.") 
    }
    Mask <- do.call(rbind, Mask)
    Mask_sums <- colSums(Mask)
    need_Mask <- !all(Mask_sums %in% c(0, nrow(Mask)))
    Mask <- apply(Mask, 2, all)

    # `mesh`, `spde`, `Amat`
    mesh <- results_mm[[1]]$mesh
    if (need_Mask) {
      mesh <- retro_mask_mesh(mesh, Mask[results_mm[[1]]$mask])
    }
    if (use_INLA) {
      spde <- INLA::inla.spde2.matern(mesh)
      Amat <- INLA::inla.spde.make.A(mesh) #Psi_{km} (for one task and subject, a VxN matrix, V=num_vox, N=num_mesh)
      Amat <- Amat[mesh$idx$loc,]
    } else {
      stop()
    }
    Amat.tot <- bdiag(rep(list(Amat), nK)) #Psi_m from paper (VKxNK)

    # Collecting theta posteriors from subject models
    Qmu_theta <- Q_theta <- 0
    # Collecting X and y cross-products from subject models (for posterior distribution of beta)
    Xcros.all <- Xycros.all <- vector("list", nN)
    for (nn in seq(nN)) {
      # [TO DO] test this
      if (need_Mask) {
        results_mm[[nn]] <- retro_mask_BGLM(
          results_mm[[nn]], Mask[results_mm[[nn]]$mask]
        )
      }

      # Check that mesh has same neighborhood structure
      if (!all.equal(results_mm[[nn]]$mesh$faces, mesh$faces, check.attribute=FALSE)) {
        stop(paste0(
          'Subject ', nn,
          ' does not have the same mesh neighborhood structure as subject 1.',
          ' Check meshes for discrepancies.'
        ))
      }

      #Collect posterior mean and precision of hyperparameters
      mu_theta_mm <- results_mm[[nn]]$INLA_model_obj$misc$theta.mode
      Q_theta_mm <- solve(results_mm[[nn]]$INLA_model_obj$misc$cov.intern)
      #iteratively compute Q_theta and mu_theta (mean and precision of q(theta|y))
      Qmu_theta <- Qmu_theta + as.vector(Q_theta_mm%*%mu_theta_mm)
      Q_theta <- Q_theta + Q_theta_mm
      rm(mu_theta_mm, Q_theta_mm)

      # compute Xcros = Psi'X'XPsi and Xycros = Psi'X'y
      # (all these matrices for a specific subject mm)
      y_vec <- results_mm[[nn]]$y
      X_list <- results_mm[[nn]]$X
      if (length(X_list) > 1) {
        n_sess <- length(X_list)
        X_list <- Matrix::bdiag(X_list)
        Amat.final <- Matrix::bdiag(rep(list(Amat.tot),n_sess))
      } else {
        X_list <- X_list[[1]]
        Amat.final <- Amat.tot
      }
      Xmat <- X_list %*% Amat.final
      Xcros.all[[nn]] <- Matrix::crossprod(Xmat)
      Xycros.all[[nn]] <- Matrix::crossprod(Xmat, y_vec)
    }
    rm(results_mm, y_vec, X_list, Xmat) # save memory

    if (use_INLA) {
      mu_theta <- solve(Q_theta, Qmu_theta) #mu_theta = poterior mean of q(theta|y) (Normal approximation) from paper, Q_theta = posterior precision
      #### DRAW SAMPLES FROM q(theta|y)
      #theta.tmp <- mvrnorm(nsamp_theta, mu_theta, solve(Q_theta))
      if (verbose>0) cat(paste0('Sampling ',nsamp_theta,' posterior samples of thetas \n'))
      theta.samp <- INLA::inla.qsample(n=nsamp_theta, Q = Q_theta, mu = mu_theta)
      #### COMPUTE WEIGHT OF EACH SAMPLES FROM q(theta|y) BASED ON PRIOR
      if (verbose>0) cat('Computing weights for each theta sample \n')
      logwt <- rep(NA, nsamp_theta)
      for (tt in seq(nsamp_theta)) {
        logwt[tt] <- F.logwt(theta.samp[,tt], spde, mu_theta, Q_theta, nN)
      }
      #weights to apply to each posterior sample of theta
      wt.tmp <- exp(logwt - max(logwt))
      wt <- wt.tmp/(sum(wt.tmp))
    } else {
      # theta.samp <- qsample(n=nsamp_theta, Q = Q_theta, mu = mu_theta) # ?
      mu_theta <- mu_theta / nN
      theta.samp <- as.matrix(mu_theta)
      wt <- 1
    }

    #get posterior quantities of beta, conditional on a value of theta
    if (verbose>0) cat(paste0('Sampling ',nsamp_beta,' betas for each value of theta \n'))
    if (is.null(num_cores)) {
      #6 minutes in simuation
      beta.posteriors <- apply(
        theta.samp,
        MARGIN = 2,
        FUN = beta.posterior.thetasamp,
        spde = spde,
        Xcros = Xcros.all,
        Xycros = Xycros.all,
        contrasts = contrasts,
        quantiles = quantiles,
        excursion_type = excursion_type,
        gamma = gamma,
        alpha = alpha,
        nsamp_beta = nsamp_beta
      )
    } else {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop(
          "`BayesGLM2` requires the `parallel` package. Please install it.",
          call. = FALSE
        )
      }

      #2 minutes in simulation (4 cores)
      max_num_cores <- min(parallel::detectCores() - 1, 25)
      num_cores <- min(max_num_cores, num_cores)
      cl <- parallel::makeCluster(num_cores)
      beta.posteriors <- parallel::parApply(
        cl, theta.samp,
        MARGIN=2,
        FUN=beta.posterior.thetasamp,
        spde=spde,
        Xcros = Xcros.all,
        Xycros = Xycros.all,
        contrasts=contrasts,
        quantiles=quantiles,
        excursion_type=excursion_type,
        gamma=gamma,
        alpha=alpha,
        nsamp_beta=nsamp_beta
      )
      parallel::stopCluster(cl)
    }

    ## Sum over samples using weights

    if (verbose>0) cat('Computing weighted summaries over beta samples \n')

    ## Posterior mean of each contrast
    betas.all <- lapply(beta.posteriors, function(x) return(x$mu))
    betas.wt <- mapply(
      function(x, a){return(x*a)},
      betas.all, wt, SIMPLIFY=FALSE
    ) #apply weight to each element of betas.all (one for each theta sample)
    betas.summ <- apply(abind::abind(betas.wt, along=3), MARGIN = c(1,2), sum)  #N x L (# of contrasts)
    dimnames(betas.summ) <- NULL

    ## Posterior quantiles of each contrast
    num_quantiles <- length(quantiles)
    if(num_quantiles > 0){
      quantiles.summ <- vector('list', num_quantiles)
      names(quantiles.summ) <- quantiles
      for(iq in 1:num_quantiles){
        quantiles.all_iq <- lapply(beta.posteriors, function(x) return(x$quantiles[[iq]]))
        betas.wt_iq <- mapply(function(x, a){return(x*a)}, quantiles.all_iq, wt, SIMPLIFY=FALSE) #apply weight to each element of quantiles.all_iq (one for each theta sample)
        quantiles.summ[[iq]] <- apply(abind::abind(betas.wt_iq, along=3), MARGIN = c(1,2), sum)  #N x L (# of contrasts)
        dimnames(quantiles.summ[[iq]]) <- NULL
      }
    } else {
      quantiles.summ <- NULL
    }

    ## Posterior probabilities and activations
    if(do_excur){
      ppm.all <- lapply(beta.posteriors, function(x) return(x$F))
      ppm.wt <- mapply(function(x, a){return(x*a)}, ppm.all, wt, SIMPLIFY=FALSE) #apply weight to each element of ppm.all (one for each theta sample)
      ppm.summ <- apply(abind::abind(ppm.wt, along=3), MARGIN = c(1,2), sum) #N x L (# of contrasts)
      dimnames(ppm.summ) <- NULL
      active <- array(0, dim=dim(ppm.summ))
      for (cc in seq(nC)) { active[ppm.summ[,cc] > (1-alpha[cc]),cc] <- 1 }
    } else {
      ppm.summ <- active <- NULL
    }

    ### Save results
    out[[mm]] <- list(
      estimates = betas.summ,
      quantiles = quantiles.summ,
      ppm = ppm.summ,
      active = active,
      mask = Mask,
      Amat = Amat # not Amat.final?
    )

    if (nM>1) { cat("~~~~~~~~~~~~~~~~~~~~\n\n") }
  }

  out <- list(
    model_results = out,
    contrasts = contrasts,
    excursion_type = excursion_type,
    task_names=task_names,
    session_names=session_names,
    gamma = gamma,
    alpha = alpha,
    nsamp_theta = nsamp_theta,
    nsamp_beta = nsamp_beta
  )
  class(out) <- "BayesGLM2"

  if (is_cifti) {
    out <- list(
      contrast_estimates_xii = as.xifti(
        out$model_results$cortex_left$estimates,
        out$model_results$cortex_left$mask,
        out$model_results$cortex_right$estimates,
        out$model_results$cortex_right$mask
      ),
      activations_xii = NULL,
      BayesGLM2_results = out
    )
    out$contrast_estimates_xii$meta$cifti$names <- names(contrasts)
    if (do_excur) {
      act_xii <- as.xifti(
        out$BayesGLM2_results$model_results$cortex_left$active,
        out$BayesGLM2_results$model_results$cortex_left$mask,
        out$BayesGLM2_results$model_results$cortex_right$active,
        out$BayesGLM2_results$model_results$cortex_right$mask
      )
      out$activations_xii <- convert_xifti(act_xii, "dlabel", colors='red')
      names(out$activations_xii$meta$cifti$labels) <- names(contrasts)
    }
    class(out) <- "BayesGLM2_cifti"
  }

  out
}

#' @rdname BayesGLM2
#' @export
BayesGLM_group <- function(
  results,
  contrasts = NULL,
  quantiles = NULL,
  excursion_type=NULL,
  gamma = 0,
  alpha = 0.05,
  nsamp_theta = 50,
  nsamp_beta = 100,
  num_cores = NULL,
  verbose = 1){

  BayesGLM2(
    results=results,
    contrasts=contrasts,
    quantiles=quantiles,
    excursion_type=excursion_type,
    gamma=gamma, alpha=alpha,
    nsamp_theta=nsamp_theta, nsamp_beta=nsamp_beta,
    num_cores=num_cores, verbose=verbose
  )
}
