#' Group-level Bayesian GLM
#'
#' Performs group-level Bayesian GLM estimation and inference using the joint
#'  approach described in Mejia et al. (2020).
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param results Either (1) a length \eqn{N} list of \code{"BGLM"} objects,
#'  or (2) a length \eqn{N} character vector of files storing \code{"BGLM"}
#'  objects saved with \code{\link{saveRDS}}. \code{"fit_bglm"} objects
#'  also are accepted.
#' @param contrasts (Optional) A list of contrast vectors that specify the
#'  group-level summaries of interest. If \code{NULL} (DEFAULT), use contrasts that
#'  compute the average of each field (field HRF) across all subjects/sessions.
#'
#'  Each contrast vector is length \eqn{KSN} specifying a group-level summary of
#'  interest, where \eqn{K} is the number of fields in the first-level design
#'  matrices, \eqn{S} is the number of sessions, and \eqn{N} is the number of
#'  subjects. The vector is grouped by fields, then sessions, then subjects.
#'
#'  For a single session/subject, the contrast vector for the first field would be:
#'
#'  \code{c0 <- c(1, rep(0, K-1)) #indexes the first field for a single session}
#'
#'  so the full contrast vector for the group *average over all sessions/subjects
#'  for the first field* would be:
#'
#'  \code{contrasts = rep(c0, S*N) /(S*N)}.
#'
#'  To obtain the group average for the first field, for *just the first session*,
#'  input zeros for the remaining sessions:
#'
#'  \code{c2 <- c(c0, rep(0, K*(S-1)))}
#'  \code{contrasts = rep(c2, N) /N}.
#'
#'  To obtain the group mean *difference between two sessions* (\eqn{S=2}) for the first field:
#'
#'  \code{c3 <- c(c0, -c0)}
#'  \code{contrasts = rep(c3, N) / N}.
#'
#'  To obtain the *mean over sessions* of the first field, just for the first subject:
#'
#'  \code{c4 <- rep(c0, S)}
#'  \code{c(c4, rep(0, K*S*(N-1))) / S}.
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
#' @importFrom ciftiTools as.xifti
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
    stop("`results` must be a list of all `'BGLM'` or all `'fit_bglm'` objects, or a character vector of files with `'BGLM(0)'` results.")
  }
  is_BGLM <- all(vapply(results, inherits, FALSE, "fit_bglm"))
  is_cifti <- all(vapply(results, inherits, FALSE, "BGLM"))
  if (!is_BGLM && !is_cifti) {
    stop("`results` must be a list of all `'BGLM'` or all `'fit_bglm'` objects, or a character vector of files with `'BGLM(0)'` results.")
  }
  rm(is_BGLM) # use `is_cifti`

  model_names <- if (is_cifti) {
    names(results[[1]]$BGLMs)[!vapply(results[[1]]$BGLMs, is.null, FALSE)]
  } else {
    "BayesGLM"
  }

  nM <- length(model_names)                 # models
  nN <- length(results)                     # subjects
  nS <- length(results[[1]]$session_names)  # sessions
  nK <- length(results[[1]]$field_names)     # fields

  session_names <- results[[1]]$session_names
  field_names <- results[[1]]$field_names

  # Check that every subject has the same models, sessions, fields
  for (nn in seq(nN)) {
    sub_nn <- results[[nn]]
    if (is_cifti) {
      stopifnot(identical(
        model_names,
        names(sub_nn$BGLMs)[!vapply(sub_nn$BGLMs, is.null, FALSE)]
      ))
    }
    stopifnot(identical(session_names, sub_nn$session_names))
    stopifnot(identical(field_names, sub_nn$field_names))
  }

  # Check `contrasts`.
  # `contrasts` should be fields * sessions * subjects
  if(!is.null(contrasts) & !is.list(contrasts)) contrasts <- list(contrasts)
  if(is.null(contrasts)) {
    if (verbose>0) cat('Using a contrast that computes the average across subjects for each field. If other contrasts are desired, provide `contrasts`.\n')
    contrasts <- vector('list', length=nK)
    names(contrasts) <- paste0(field_names, '_avg')
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

  # Get intersection mask.
  Masks <- intersect_mask(results)

  # If `BGLM` object, we will only be using the `BGLMs` list entry.
  # Delete everything else from here on for clarity.
  if (is_cifti) { results <- lapply(results, '[', "BGLMs") }

  spatial_sub <- NULL # only used for subcortex model

  # Do the group model
  for (mm in seq(nM)) {
    Mask <- Masks$Mdat[[mm]]

    if (nM>1) { if (verbose>0) cat(model_names[mm], " ~~~~~~~~~~~\n") }
    results_mm <- if (is_cifti) {
      lapply(results, function(x){ x$BGLMs[[mm]] })
    } else {
      results
    }

    # We know model names match, but still check `spatial_type` match.
    spatial_type <- vapply(results_mm, function(x){x$spatial$spatial_type}, "")
    spatial_type <- unique(spatial_type)
    if (length(spatial_type) != 1) {
      stop("`spatial_type` is not unique across subjects for model ", model_names[mm], ".")
    }

    # Get new `spatial`, `spde`, and `Amat`.
    mesh <- NULL # only used for vertex model
    if (spatial_type == "vertex") {
      # `spatial`.
      spatial <- results_mm[[1]]$spatial
      spatial$maskMdat <- Mask
      spatial$Mmap <- which(Mask)

      # `mesh` and `spde`.
      mesh <- results_mm[[1]]$spde$mesh
      spde <- results_mm[[1]]$spde

      # `Amat`
      Amat <- INLA::inla.spde.make.A(mesh) #Psi_{km} (for one field and subject, a VxN matrix, V=num_vox, N=num_mesh)
      Amat <- Amat[mesh$idx$loc,]

    } else if (spatial_type == "voxel") {
      # Check features of spatial that are expected to match for all subjects.
      spatials_expEq <- unique(lapply(results_mm, function(x){ x$spatial[c(
        "spatial_type", "labels",
        "trans_mat", "trans_units",
        "nbhd_order", "buffer")] }))

      if (length(unique(spatials_expEq)) != 1) {
        stop("`spatial`s for voxel model are expected to match in labels, trans_mat, trans_units, nbhd_order, and buffer.")
      }
      rm(spatials_expEq)

      # Get.
      spatial <- results_mm[[1]]$spatial

      # Update.
      toKeep <- Mask[results_mm[[1]]$spatial$maskMdat]
      spatial$labsMdat <- spatial$labsMdat[toKeep]
      spatial$maskMdat[spatial$maskMdat] <- toKeep
      spatial$Mmap <- spatial$Mmap[toKeep]
      spatial_sub <- spatial # for making the output xifti.
      if (length(unique(lapply(results_mm, function(q){q$logkappa_vec}))) > 1) {
        warning("`logkappa_vec` is no the same across subjects. Using the first subject's.")
      }
      if (length(unique(lapply(results_mm, function(q){q$logkappa_vec}))) > 1) {
        warning("`logtau_vec` is no the same across subjects. Using the first subject's.")
      }
      x <- SPDE_from_voxel(
        spatial,
        logkappa = results_mm[[1]]$logkappa_vec,
        logtau = results_mm[[1]]$logtau_vec
      )
      spde <- x$spde
      spatial <- x$spatial

      Amat <- make_A_mat(results_mm[[1]]$spatial)
    }

    Amat.tot <- bdiag(rep(list(Amat), nK)) #Psi_m from paper (VKxNK)

    # Update the results with the new SPDE.
    for (nn in seq(nN)) {
      cat(paste0("Checking data mask for subject ", nn, ".\n")) # [TO DO] option to hide?
      # [NOTE] for subcortex, we need to see the old `spatial` in order to
      # update `X`. So update `spatial` after `retro_mask_fit_bglm`, not before.
      results_mm[[nn]] <- retro_mask_fit_bglm(results_mm[[nn]], Mask)
      results_mm[[nn]]$spde <- spde
      results_mm[[nn]]$spatial <- spatial
    }

    # Collecting theta posteriors from subject models
    Qmu_theta <- Q_theta <- 0
    # Collecting X and y cross-products from subject models (for posterior distribution of beta)
    Xcros.all <- Xycros.all <- vector("list", nN)
    for (nn in seq(nN)) {
      # Check that mesh has same neighborhood structure
      if (!is.null(mesh)) {
        if (!all.equal(results_mm[[nn]]$spatial$mesh$faces, mesh$faces, check.attribute=FALSE)) {
          stop(paste0(
            'Subject ', nn,
            ' does not have the same mesh neighborhood structure as subject 1.',
            ' Check meshes for discrepancies.'
          ))
        }
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
        X_list <- Matrix::bdiag(X_list) #block-diagonialize over sessions
        Amat.final <- Matrix::bdiag(rep(list(Amat.tot),n_sess))
      } else {
        X_list <- X_list[[1]] #single-session case
        Amat.final <- Amat.tot
      }
      Xmat <- X_list #%*% Amat.final #already done within BayesGLM
      Xcros.all[[nn]] <- Matrix::crossprod(Xmat)
      Xycros.all[[nn]] <- Matrix::crossprod(Xmat, y_vec)
    }
    #rm(results_mm, y_vec, X_list, Xmat) # save memory

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

    # # Above, but trying to not use INLA.
    # # theta.samp <- qsample(n=nsamp_theta, Q = Q_theta, mu = mu_theta) # ?
    # mu_theta <- mu_theta / nN
    # theta.samp <- as.matrix(mu_theta)
    # wt <- 1

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

      if (verbose>0) cat(paste0('\t ... running in parallel with ',num_cores,' cores \n'))

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
      estimates = betas.summ, #includes boundary locations
      quantiles = quantiles.summ,
      ppm = ppm.summ,
      active = active,
      mask = lapply(Masks, '[[', mm),
      Amat = Amat # not Amat.final?
    )

    if (nM>1) { cat("\n") }
  }

  out <- list(
    model_results = out,
    contrasts = contrasts,
    excursion_type = excursion_type,
    field_names=field_names,
    session_names=session_names,
    gamma = gamma,
    alpha = alpha,
    nsamp_theta = nsamp_theta,
    nsamp_beta = nsamp_beta
  )
  class(out) <- "fit_bglm2"

  if (is_cifti) {

    # Set values in maskIn but not maskMdat to `NA`.
    # Mask with maskIn.
    result_oomSetNA <- out$model_results
    for (mm in seq(nM)) {
      results_mm <- lapply(results, function(x){ x$BGLMs[[mm]] })
      spatial_type <- unique(
        vapply(results_mm, function(x){x$spatial$spatial_type}, "")
      )

      if (spatial_type == "vertex") {
        result_oomSetNA[[mm]]$estimates[Masks$In[[mm]] & (!Masks$Mdat[[mm]]),] <- NA
        result_oomSetNA[[mm]]$estimates <- result_oomSetNA[[mm]]$estimates[Masks$In[[mm]],,drop=FALSE]
      } else {
        result_oomSetNA[[mm]]$estimates <- unmask_Mdat2In(
          result_oomSetNA[[mm]]$estimates[spatial_sub$Mmap,,drop=FALSE],
          spatial_sub$maskMdat[],
          spatial_sub$maskIn[]
        )
      }
    }

    out <- list(
      contrast_estimate_xii = as.xifti(
        cortexL = result_oomSetNA$cortexL$estimates,
        cortexL_mwall = Masks$In$cortexL,
        cortexR = result_oomSetNA$cortexR$estimates,
        cortexR_mwall = Masks$In$cortexR,
        c(NA, NaN),
        subcortVol = result_oomSetNA$subcort$estimates,
        subcortLabs = spatial_sub$labels,
        subcortMask = spatial_sub$maskIn
      ),
      activations_xii = NULL,
      masks = Masks,
      BayesGLM2_results = out
    )
    out$contrast_estimate_xii$meta$cifti$names <- names(contrasts)

    if (do_excur) {

      # Set values in maskIn but not maskMdat to `NA`.
      # Mask with maskIn.
      result_oomSetNA <- out$BayesGLM2_results$model_results
      for (mm in seq(nM)) {
        results_mm <- lapply(results, function(x){ x$BGLMs[[mm]] })
        spatial_type <- unique(
          vapply(results_mm, function(x){x$spatial$spatial_type}, "")
        )

        if (spatial_type == "vertex") {
          result_oomSetNA[[mm]]$estimates[Masks$In[[mm]] & (!Masks$Mdat[[mm]]),] <- NA
          result_oomSetNA[[mm]]$estimates <- result_oomSetNA[[mm]]$estimates[Masks$In[[mm]],,drop=FALSE]
        } else {
          result_oomSetNA[[mm]]$estimates <- unmask_Mdat2In(
            result_oomSetNA[[mm]]$estimates[spatial_sub$Mmap,,drop=FALSE],
            spatial_sub$maskMdat[],
            spatial_sub$maskIn[]
          )
        }
      }

      act_xii <- as.xifti(
        cortexL = result_oomSetNA$cortexL$active,
        cortexL_mwall = Masks$In$cortexL,
        cortexR = result_oomSetNA$cortexR$active,
        cortexR_mwall = Masks$In$cortexR,
        c(NA, NaN),
        subcortVol = result_oomSetNA$subcort$active,
        subcortLabs = spatial_sub$labsMdat,
        subcortMask = spatial_sub$maskMdat
      )
      out$activations_xii <- convert_xifti(act_xii, "dlabel", colors='red')
      out$activations_xii$meta$cifti$names <- names(contrasts)
      names(out$activations_xii$meta$cifti$labels) <- names(contrasts)
    }
    class(out) <- "BGLM2"
  }
  out
}

