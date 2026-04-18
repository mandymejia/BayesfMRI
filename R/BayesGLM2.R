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
#' @param design_matrix Design matrix for group-level summaries of interest. The
#' number of rows must be equal to \eqn{K\times S\times N}, where \eqn{K} is the
#' number of fields (task regressors) in the first-level design matrices, \eqn{S}
#' is the number of sessions per subject, and \eqn{N} is the number of subjects.
#' The rows are grouped by fields, then sessions, then subjects. The number of
#' columns, \eqn{p}, depends on the group-level design. If no intercept is included,
#' it will be added as long as it is not linearly dependent with the provided design.
#' @param contrast_matrix Matrix of group-level contrasts of interests, with \eqn{p}
#' columns and the number of rows corresponding to the number of contrasts to
#' estimate.
#' @param contrast_list (Only used if design_matrix and contrast_matrix not provided)
#'  A list of contrast vectors that specify the group-level summary or summaries of
#'  interest. If \code{NULL} (DEFAULT), a contrast will be automatically constructed
#'  to compute the average of each task regressor across all subjects/sessions.
#'
#'  Contrast vectors must be of length \eqn{K\times S\times N} specifying a group
#'  -level summary of interest, where \eqn{K} is the number of fields (task
#'  regressors) in the first-level design matrices, \eqn{S} is the number of
#'  sessions per subject, and \eqn{N} is the number of subjects. The vector is
#'  grouped by fields, then sessions, then subjects.
#'
#'  See Details for examples of contrast vectors for different group level summaries.
#'
#' @param contrasts Deprecated. Use `contrast_list` instead.
#' @param quantiles (Optional) Vector of posterior quantiles to return in
#'  addition to the posterior mean.
#' @param excursion_type (For inference only) The type of excursion function for
#'  the contrast (">", "<", "!="), or a vector thereof (each element
#'  corresponding to one contrast).  If \code{NULL}, no inference performed.
#' @param contrast_names (Optional) Names of contrasts.
#' @param gamma (For inference only) Activation threshold for the excursion set,
#'  or a vector thereof (each element corresponding to one contrast). Default:
#'  \code{0}.
#' @param alpha (For inference only, single-threshold mode only) Significance
#'  level for activation for the excursion set, or a vector thereof (each
#'  element corresponding to one contrast). Default: \code{.05}. This argument
#'  is ignored when \code{alpha_nested} is provided.
#' @param alpha_nested (For inference only, nested-threshold mode) Optional
#'  nested significance thresholds used to construct multiple activation maps
#'  and a nested activation map. If \code{NULL} (default), the function uses
#'  \code{alpha} and returns a single activation map for each contrast.
#'  Default: \code{NULL}.If provided, \code{alpha} is ignored.
#'
#'  \code{alpha_nested} can be specified in either of the following forms:
#'  \itemize{
#'    \item A numeric vector, in which case the same set of alpha thresholds is
#'    used for all contrasts.
#'    \item A list of numeric vectors of length equal to the number of contrasts,
#'    in which case each contrast uses its own set of alpha thresholds.
#'  }
#'
#'  Each alpha threshold must be between 0 and 1. Thresholds should be ordered
#'  from less strict to more strict (for example, \code{c(0.1, 0.05, 0.01)}),
#'  although the function may internally reorder and deduplicate them.
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
#' @details The easiest way to compute group-level contrasts is to specify the
#' arguments `design_matrix` and `contrast_matrix`.  If you wish to manually specify
#' contrasts instead, you can specify the `contrast_list` argument. Some examples
#' of contrast vectors for specific group-level summaries:
#'
#' **Example 1:** For a single session/subject, the contrast vector for the first field would be:
#'
#'  \code{c0 <- c(1, rep(0, K-1)) #indexes the first field for a single session}
#'
#'  so the full contrast vector for the group *average over all sessions/subjects
#'  for the first field* would be:
#'
#'  \code{contrast_list = rep(c0, S*N) /(S*N)}.
#'
#'  **Example 2:** To obtain the group average for the first field, for *just the first session*,
#'  input zeros for the remaining sessions:
#'
#'  \code{c2 <- c(c0, rep(0, K*(S-1)));}
#'  \code{contrast_list = rep(c2, N) /N}.
#'
#'  **Example 3:** To obtain the group mean *difference between two sessions* (\eqn{S=2}) for the first field:
#'
#'  \code{c3 <- c(c0, -c0);}
#'  \code{contrast_list = rep(c3, N) / N}.
#'
#'  **Example 4:** To obtain the *mean over sessions* of the first field, just for the first subject:
#'
#'  \code{c4 <- rep(c0, S);}
#'  \code{contrast_list = c(c4, rep(0, K*S*(N-1))) / S}.
#'
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag crossprod
#' @importFrom ciftiTools as.xifti
#'
#' @export
BayesGLM2 <- function(
    results,
    design_matrix = NULL,
    contrast_matrix = NULL,
    contrast_list = NULL,
    contrasts = NULL,
    quantiles = NULL,
    excursion_type = ">",
    contrast_names = NULL,
    gamma = 0,
    alpha = 0.05,
    alpha_nested = NULL,
    nsamp_theta = 50,
    nsamp_beta = 100,
    num_cores = NULL,
    verbose = 1){

  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("`BayesGLM2` requires the `abind` package. Please install it.", call. = FALSE)
  }

  # Check `results`, but do not read all files into memory at once.
  results_from_files <- is.character(results)

  if (results_from_files) {
    if (!all(file.exists(results))) {
      stop("`results` is a character vector, but not all elements are existing files.")
    }
  } else if (!is.list(results)) {
    stop("`results` must be a list of all `'BGLM'` or all `'fit_bglm'` objects, or a character vector of files with `'BGLM(0)'` results.")
  }

  .get_result <- function(nn) {
    if (results_from_files) {
      readRDS(results[[nn]])
    } else {
      results[[nn]]
    }
  }

  .get_model <- function(res, mm, is_cifti) {
    if (is_cifti) {
      res$BGLMs[[mm]]
    } else {
      res
    }
  }

  # Check the class of the first result to determine the object type and whether it's CIFTI. We will check that all other results match this type in the loop below.
  res1 <- .get_result(1)

  if (inherits(res1, "BGLM")) {
    object_type <- "BGLM"
    is_cifti <- TRUE
  } else if (inherits(res1, "fit_bglm")) {
    object_type <- "fit_bglm"
    is_cifti <- FALSE
  } else {
    stop("`results` must be a list of all `'BGLM'` or all `'fit_bglm'` objects, or a character vector of files with `'BGLM(0)'` results.")
  }

  model_names <- if (is_cifti) {
    names(res1$BGLMs)[!vapply(res1$BGLMs, is.null, FALSE)]
  } else {
    "BayesGLM"
  }

  nM <- length(model_names)                # models (brain structures, "cortexL", "cortexR", "subcort")
  nN <- length(results)                    # subjects
  nS <- length(res1$session_names)         # sessions
  nK <- length(res1$field_names)           # fields

  session_names <- res1$session_names
  field_names <- res1$field_names

  # Check that every subject has the same object type, models, sessions, and fields
  for (nn in seq_len(nN)) {
    sub_nn <- .get_result(nn)

    if (object_type == "BGLM") {
      if (!inherits(sub_nn, "BGLM")) {
        stop("Subject ", nn, " is not a 'BGLM' object.")
      }
      stopifnot(identical(
        model_names,
        names(sub_nn$BGLMs)[!vapply(sub_nn$BGLMs, is.null, FALSE)]
      ))
    } else if (object_type == "fit_bglm") {
      if (!inherits(sub_nn, "fit_bglm")) {
        stop("Subject ", nn, " is not a 'fit_bglm' object.")
      }
    }

    stopifnot(identical(session_names, sub_nn$session_names))
    stopifnot(identical(field_names, sub_nn$field_names))

    rm(sub_nn)
    gc(FALSE)
  }

  # Yunong added on Mar 02, 2026
  # Check `design_matrix` and `contrast_matrix`
  if (!is.null(design_matrix) || !is.null(contrast_matrix)) {

    if (xor(!is.null(design_matrix), !is.null(contrast_matrix))) {
      stop(
        "You must provide both `design_matrix` and `contrast_matrix`, or provide neither.\n",
        "If you want to specify contrasts manually, provide `contrast_list` instead.",
        call. = FALSE
      )
    }

    # contrast_list = (C %*% (X'X)^{-1} %*% X')
    # Here X = design_matrix, C = contrast_matrix
    if (!is.null(contrast_list)) {
      warning(
        "`design_matrix` and `contrast_matrix` were provided, so `contrast_list` will be ignored.",
        call. = FALSE
      )
      contrast_list <- NULL
    }

    # Validate numeric matrices
    design_matrix <- as.matrix(design_matrix)
    contrast_matrix <- as.matrix(contrast_matrix)

    if (!is.numeric(design_matrix)) stop("`design_matrix` must be numeric.", call. = FALSE)
    if (!is.numeric(contrast_matrix)) stop("`contrast_matrix` must be numeric.", call. = FALSE)

    J <- nK * nS * nN
    if (nrow(design_matrix) != J) {
      stop(
        sprintf(
          "`design_matrix` must have %d rows (= K*S*N = %d*%d*%d), but has %d.",
          J, nK, nS, nN, nrow(design_matrix)
        ),
        call. = FALSE
      )
    }

    # Check dimension match before adding intercept.
    if (ncol(contrast_matrix) != ncol(design_matrix)) {
      stop(
        sprintf(
          "`contrast_matrix` must have %d columns to match `design_matrix`, but has %d.",
          ncol(design_matrix), ncol(contrast_matrix)
        ),
        call. = FALSE
      )
    }

    # Try adding an intercept to `design_matrix`.
    added_intercept <- FALSE
    has_intercept <- any(apply(
      design_matrix, 2,
      function(col) isTRUE(all(abs(col - 1) < .Machine$double.eps^0.5))
    ))

    if (!has_intercept) {
      X0 <- cbind(`(Intercept)` = 1, design_matrix)
      # Check if adding the intercept increases the rank of `design_matrix`.
      if (qr(X0)$rank > qr(design_matrix)$rank) {
        design_matrix <- X0
        added_intercept <- TRUE
      }
    }

    # If we added an intercept to `design_matrix`, add a zero intercept column to `contrast_matrix`.
    if (added_intercept) {
      contrast_matrix <- cbind(`(Intercept)` = 0, contrast_matrix)
    }

    # Construct `contrast_list` from `design_matrix` and `contrast_matrix`.
    XtX <- crossprod(design_matrix)
    XtX_inv <- tryCatch(
      solve(XtX),
      error = function(e) stop("`design_matrix` is rank-deficient: cannot form (X'X)^{-1}.", call. = FALSE)
    )
    A <- contrast_matrix %*% XtX_inv %*% t(design_matrix)

    # Each row of A is a contrast vector of length J (=K*S*N)
    contrast_list <- lapply(seq_len(nrow(A)), function(i) as.numeric(A[i, ]))
    if (!is.null(rownames(contrast_matrix))) {
      names(contrast_list) <- rownames(contrast_matrix)
    } else {
      names(contrast_list) <- paste0("contrast_", seq_len(nrow(A)))
    }
  }


  # Check `contrast_list`.
  # Check `contrasts`.
  if (!is.null(contrasts)) {
    .Deprecated(msg = paste0(
      "`contrasts` is deprecated in BayesGLM2().\n",
      "Use `contrast_list` instead."
    ))
    if (is.null(contrast_list)) {
      contrast_list <- contrasts
    } else {
      stop("`contrasts` and `contrast_list` were both provided. Use just `contrast_list` instead.")
    }
  }
  # `contrast_list` should be fields * sessions * subjects
  if(!is.null(contrast_list) & !is.list(contrast_list)) contrast_list <- list(contrast_list)
  if(is.null(contrast_list)) {
    if (verbose>0) cat('Computing the average across subjects for each field. If other contrasts are desired, please provide `design_matrix` and `contrast_matrix`, or `contrast_list`.\n')
    contrast_list <- vector('list', length=nK)
    names(contrast_list) <- paste0(field_names, '_avg')
    for (kk in 1:nK) {
      # (1/J, 0, 0, ..., 0) for k=1,
      # (0, 1/J, 0, ..., 0) for k=2,
      # ...,
      # (0, 0, ..., 0, 1/J) for k=K
      # for each session, for each subject
      # where J == S * N
      contrast_1 <- c(rep(0, kk-1), 1/(nS*nN), rep(0, nK-kk)) # length nK
      contrast_list[[kk]] <- rep(rep(contrast_1, nS), nN)         # length nK*nS*nN
    }
  } else {
    #Check that each contrast vector is numeric and length J*K
    if(any(sapply(contrast_list, length) != nK*nS*nN)) {
      stop('Each contrast vector must be of length K*S*N (fields times sessions times subjects).')
    }
    if(any(!sapply(contrast_list, is.numeric))) {
      stop('Each contrast vector must be numeric, but at least one is not.')
    }
    if (is.null(names(contrast_list))) {
      names(contrast_list) <- paste0("contrast_", seq(length(contrast_list)))
    }
  }
  # Override `names(contrast_list)` with `contrast_names` if provided.
  if (!is.null(contrast_names)) {
    stopifnot(length(contrast_names) == length(contrast_list))
    names(contrast_list) <- contrast_names
  }
  nC <- length(contrast_list)

  # Check `quantiles`
  if(!is.null(quantiles)){
    stopifnot(is.numeric(quantiles))
    if(any(quantiles > 1 | quantiles < 0)) stop('All elements of `quantiles` must be between 0 and 1.')
  }

  # Check alpha mode
  use_nested_alpha <- !is.null(alpha_nested)
  if (use_nested_alpha && verbose > 0) {
    message("`alpha_nested` provided: ignoring `alpha` and using nested activation mode.")
  }

  # Validate and normalize `alpha_nested` if provided, and override `alpha` if `alpha_nested` is used.a
  .normalize_alpha_nested <- function(alpha_nested, nC, contrast_names = NULL) {
    if (is.numeric(alpha_nested)) {
      alpha_grid <- replicate(
        nC,
        sort(unique(alpha_nested), decreasing = TRUE),
        simplify = FALSE
      )
    } else if (is.list(alpha_nested)) {
      if (length(alpha_nested) != nC) {
        stop(
          "`alpha_nested` must be either a numeric vector or a list of length equal to the number of contrasts.",
          call. = FALSE
        )
      }
      alpha_grid <- lapply(alpha_nested, function(x) {
        if (!is.numeric(x)) {
          stop("Each element of `alpha_nested` must be numeric.", call. = FALSE)
        }
        sort(unique(x), decreasing = TRUE)
      })
    } else {
      stop(
        "`alpha_nested` must be either a numeric vector or a list of numeric vectors.",
        call. = FALSE
      )
    }
    bad_len <- vapply(alpha_grid, length, integer(1)) == 0
    if (any(bad_len)) {
      stop("Each contrast must have at least one alpha threshold in `alpha_nested`.", call. = FALSE)
    }
    bad_val <- vapply(
      alpha_grid,
      function(x) any(!is.finite(x) | x <= 0 | x >= 1),
      logical(1)
    )
    if (any(bad_val)) {
      stop("All values in `alpha_nested` must be finite and strictly between 0 and 1.", call. = FALSE)
    }
    if (!is.null(contrast_names)) {
      names(alpha_grid) <- contrast_names
    }
    alpha_grid
  }

  do_excur <- !is.null(excursion_type) && (!identical(excursion_type, "none"))

  alpha_grid <- NULL

  if (do_excur) {
    if (length(excursion_type) == 1) excursion_type <- rep(excursion_type, nC)
    if (length(gamma) == 1) gamma <- rep(gamma, nC)

    if (length(gamma) != nC) {
      stop("Length of `gamma` must match number of contrasts or be equal to one.", call. = FALSE)
    }
    if (length(excursion_type) != nC) {
      stop("Length of `excursion_type` must match number of contrasts or be equal to one.", call. = FALSE)
    }

    if (use_nested_alpha) {
      alpha_grid <- .normalize_alpha_nested(
        alpha_nested = alpha_nested,
        nC = nC,
        contrast_names = names(contrast_list)
      )
      alpha <- NULL
    } else {
      if (length(alpha) == 1) alpha <- rep(alpha, nC)
      if (length(alpha) != nC) {
        stop("Length of `alpha` must match number of contrasts or be equal to one.", call. = FALSE)
      }
      if (any(!is.finite(alpha) | alpha <= 0 | alpha >= 1)) {
        stop("All values in `alpha` must be finite and strictly between 0 and 1.", call. = FALSE)
      }
    }
  } else {
    excursion_type <- "none"
    alpha <- NULL
    alpha_grid <- NULL
  }

  out <- vector("list", nM)
  names(out) <- model_names

  # Get intersection mask without `intersect_mask` function, to avoid reading all data into memory at once.
  Masks <- list(
    In = vector("list", nM),
    Mdat = vector("list", nM)
  )
  names(Masks$In) <- names(Masks$Mdat) <- model_names

  for (mm in seq_len(nM)) {
    mod1 <- .get_model(res1, mm, is_cifti)
    Masks$In[[mm]] <- as.logical(mod1$spatial$maskIn)
    Masks$Mdat[[mm]] <- as.logical(mod1$spatial$maskMdat)
  }

  if (nN >= 2) {
    for (nn in 2:nN) {
      sub_nn <- .get_result(nn)
      for (mm in seq_len(nM)) {
        mod_nn <- .get_model(sub_nn, mm, is_cifti)
        Masks$In[[mm]] <- Masks$In[[mm]] & as.logical(mod_nn$spatial$maskIn)
        Masks$Mdat[[mm]] <- Masks$Mdat[[mm]] & as.logical(mod_nn$spatial$maskMdat)
      }
      rm(sub_nn, mod_nn)
      gc(FALSE)
    }
  }

  spatial_sub <- NULL # only used for subcortex model
  spatial_sub_by_model <- vector("list", nM)
  names(spatial_sub_by_model) <- model_names
  spatial_type_by_model <- character(nM)
  names(spatial_type_by_model) <- model_names

  # Do the group model (looping over models, which are brain structures in the cifti case).
  for (mm in seq(nM)) {

    Mask <- Masks$Mdat[[mm]]

    if (nM > 1) { if (verbose > 0) cat(model_names[mm], " ~~~~~~~~~~~\n") }

    res1_mm <- .get_model(res1, mm, is_cifti)

    # We know model names match, but still check `spatial_type` match.
    spatial_type <- res1_mm$spatial$spatial_type
    spatial_type_by_model[mm] <- spatial_type

    for (nn in seq_len(nN)) {
      sub_nn <- .get_result(nn)
      sub_nn_mm <- .get_model(sub_nn, mm, is_cifti)
      if (!identical(sub_nn_mm$spatial$spatial_type, spatial_type)) {
        stop("`spatial_type` is not unique across subjects for model ", model_names[mm], ".")
      }
      rm(sub_nn, sub_nn_mm)
      gc(FALSE)
    }

    # Get new `spatial`, `spde`, and `Amat`.
    mesh <- NULL # only used for vertex model

    if (spatial_type == "vertex") {
      # `spatial`
      spatial <- res1_mm$spatial
      spatial$maskMdat <- Mask
      spatial$Mmap <- which(Mask)

      # `mesh` and `spde`
      mesh <- res1_mm$spde$mesh
      spde <- res1_mm$spde

      # `Amat`
      Amat <- INLA::inla.spde.make.A(mesh)
      Amat <- Amat[mesh$idx$loc,]

    } else if (spatial_type == "voxel") {
      # Check voxel spatial features that are expected to match across subjects.
      spatial_ref <- res1_mm$spatial[c(
        "spatial_type", "labels",
        "trans_mat", "trans_units",
        "nbhd_order", "buffer"
      )]

      logkappa_ref <- res1_mm$logkappa_vec
      logtau_ref <- res1_mm$logtau_vec

      for (nn in seq_len(nN)) {
        sub_nn <- .get_result(nn)
        sub_nn_mm <- .get_model(sub_nn, mm, is_cifti)

        spatial_nn <- sub_nn_mm$spatial[c(
          "spatial_type", "labels",
          "trans_mat", "trans_units",
          "nbhd_order", "buffer"
        )]

        if (!identical(spatial_ref, spatial_nn)) {
          stop("`spatial`s for voxel model are expected to match in labels, trans_mat, trans_units, nbhd_order, and buffer.")
        }

        if (!identical(sub_nn_mm$logkappa_vec, logkappa_ref)) {
          warning("`logkappa_vec` is not the same across subjects. Using the first subject's.")
          break
        }

        rm(sub_nn, sub_nn_mm, spatial_nn)
        gc(FALSE)
      }

      for (nn in seq_len(nN)) {
        sub_nn <- .get_result(nn)
        sub_nn_mm <- .get_model(sub_nn, mm, is_cifti)

        if (!identical(sub_nn_mm$logtau_vec, logtau_ref)) {
          warning("`logtau_vec` is not the same across subjects. Using the first subject's.")
          break
        }

        rm(sub_nn, sub_nn_mm)
        gc(FALSE)
      }

      # Get
      spatial <- res1_mm$spatial

      # Update
      toKeep <- Mask[res1_mm$spatial$maskMdat]
      spatial$labsMdat <- spatial$labsMdat[toKeep]
      spatial$maskMdat[spatial$maskMdat] <- toKeep
      spatial$Mmap <- spatial$Mmap[toKeep]
      spatial_sub_by_model[[mm]] <- spatial

      x <- SPDE_from_voxel(
        spatial,
        logkappa = logkappa_ref,
        logtau = logtau_ref
      )
      spde <- x$spde
      spatial <- x$spatial

      Amat <- make_A_mat(res1_mm$spatial)
    }

    Amat.tot <- bdiag(rep(list(Amat), nK)) # Psi_m from paper (VKxNK)

    # Collect theta posteriors and X/y cross-products one subject at a time.
    Qmu_theta <- Q_theta <- 0
    Xcros.all <- Xycros.all <- vector("list", nN)

    for (nn in seq_len(nN)) {
      cat(paste0("Checking data mask for subject ", nn, ".\n"))

      sub_nn <- .get_result(nn)
      res_nn_mm <- .get_model(sub_nn, mm, is_cifti)

      # [NOTE] for subcortex, we need to see the old `spatial` in order to
      # update `X`. So update `spatial` after `retro_mask_fit_bglm`, not before.
      res_nn_mm <- retro_mask_fit_bglm(res_nn_mm, Mask)
      res_nn_mm$spde <- spde
      res_nn_mm$spatial <- spatial

      # Check that mesh has same neighborhood structure
      if (!is.null(mesh)) {
        if (!all.equal(res_nn_mm$spatial$mesh$faces, mesh$faces, check.attribute = FALSE)) {
          stop(paste0(
            'Subject ', nn,
            ' does not have the same mesh neighborhood structure as subject 1.',
            ' Check meshes for discrepancies.'
          ))
        }
      }

      # Collect posterior mean and precision of hyperparameters
      mu_theta_mm <- res_nn_mm$INLA_model_obj$misc$theta.mode
      Q_theta_mm <- solve(res_nn_mm$INLA_model_obj$misc$cov.intern)

      # Iteratively compute Q_theta and mu_theta
      Qmu_theta <- Qmu_theta + as.vector(Q_theta_mm %*% mu_theta_mm)
      Q_theta <- Q_theta + Q_theta_mm

      # Compute Xcros = Psi'X'XPsi and Xycros = Psi'X'y
      y_vec <- res_nn_mm$y
      X_list <- res_nn_mm$X

      if (length(X_list) > 1) {
        n_sess <- length(X_list)
        X_list <- Matrix::bdiag(X_list) # block-diagonalize over sessions
        Amat.final <- Matrix::bdiag(rep(list(Amat.tot), n_sess))
      } else {
        X_list <- X_list[[1]] # single-session case
        Amat.final <- Amat.tot
      }

      Xmat <- X_list #%*% Amat.final # already done within BayesGLM
      Xcros.all[[nn]] <- Matrix::crossprod(Xmat)
      Xycros.all[[nn]] <- Matrix::crossprod(Xmat, y_vec)

      rm(sub_nn, res_nn_mm, mu_theta_mm, Q_theta_mm, y_vec, X_list, Xmat, Amat.final)
      gc(FALSE)
    }

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
        contrasts = contrast_list,
        quantiles = quantiles,
        excursion_type = excursion_type,
        gamma = gamma,
        alpha = alpha,
        alpha_grid = alpha_grid,
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
        contrasts=contrast_list,
        quantiles=quantiles,
        excursion_type=excursion_type,
        gamma=gamma,
        alpha=alpha,
        alpha_grid = alpha_grid,
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
    if (do_excur) {

      if (use_nested_alpha) {

        ## Nested-alpha mode:
        ## beta.posterior.thetasamp() returns F_levels instead of F.
        ## For each contrast cc, F_levels[[cc]] is an n.mesh x n_alpha_cc matrix
        ## for a single theta draw. We now average these over theta using wt.

        ppm.summ <- NULL
        active <- NULL

        ppm.levels <- vector("list", nC)
        active_levels <- vector("list", nC)
        names(ppm.levels) <- names(contrast_list)
        names(active_levels) <- names(contrast_list)

        nested_code <- matrix(
          0L,
          nrow = nrow(betas.summ),
          ncol = nC
        )
        colnames(nested_code) <- names(contrast_list)

        for (cc in seq_len(nC)) {

          ## Collect theta-specific F_levels for this contrast
          ## Each element is an n.mesh x n_alpha_cc matrix
          F.levels.cc <- lapply(beta.posteriors, function(x) x$F_levels[[cc]])

          ## Apply theta weights
          F.levels.cc.wt <- mapply(
            function(x, a) x * a,
            F.levels.cc, wt,
            SIMPLIFY = FALSE
          )

          ## Weighted average over theta
          ppm.cc <- apply(
            abind::abind(F.levels.cc.wt, along = 3),
            MARGIN = c(1, 2),
            sum
          )

          ## Make sure ppm.cc stays a matrix even if there is only one alpha
          if (is.null(dim(ppm.cc))) {
            ppm.cc <- matrix(ppm.cc, ncol = 1)
          }

          colnames(ppm.cc) <- paste0("alpha_", alpha_grid[[cc]])
          ppm.levels[[cc]] <- ppm.cc

          ## For each alpha, threshold at 1 - alpha
          aa <- alpha_grid[[cc]]
          active.cc <- vapply(
            seq_along(aa),
            function(j) as.integer(ppm.cc[, j] > (1 - aa[j])),
            integer(nrow(ppm.cc))
          )

          ## Keep matrix shape if only one alpha
          if (is.null(dim(active.cc))) {
            active.cc <- matrix(active.cc, ncol = 1)
          }

          colnames(active.cc) <- paste0("alpha_", aa)
          active_levels[[cc]] <- active.cc

          ## Nested code:
          ## Because thresholds are nested, rowSums(active.cc) gives the level.
          nested_code[, cc] <- rowSums(active.cc)
        }

      } else {

        ## Single-alpha mode: original behavior
        ppm.all <- lapply(beta.posteriors, function(x) x$F)
        ppm.wt <- mapply(
          function(x, a) x * a,
          ppm.all, wt,
          SIMPLIFY = FALSE
        )
        ppm.summ <- apply(
          abind::abind(ppm.wt, along = 3),
          MARGIN = c(1, 2),
          sum
        )
        dimnames(ppm.summ) <- NULL

        active <- array(0L, dim = dim(ppm.summ))
        for (cc in seq_len(nC)) {
          active[ppm.summ[, cc] > (1 - alpha[cc]), cc] <- 1L
        }

        ppm.levels <- NULL
        active_levels <- NULL
        nested_code <- NULL
      }

    } else {
      ppm.summ <- NULL
      ppm.levels <- NULL
      active <- NULL
      active_levels <- NULL
      nested_code <- NULL
    }

    ### Save results
    out[[mm]] <- list(
      estimates = betas.summ,            # includes boundary locations
      quantiles = quantiles.summ,

      # Single-threshold mode:
      #   ppm is an n.mesh x nC matrix
      # Nested-threshold mode:
      #   ppm is NULL, and ppm_levels stores one matrix per contrast
      ppm = ppm.summ,
      ppm_levels = ppm.levels,

      # Single-threshold mode:
      #   active is an n.mesh x nC binary matrix
      # Nested-threshold mode:
      #   active is NULL, and active_levels stores one binary matrix per contrast
      active = active,
      active_levels = active_levels,

      # Nested-threshold mode only:
      #   integer code giving the nested activation level at each location
      nested_code = nested_code,

      mask = lapply(Masks, '[[', mm),
      Amat = Amat                      # not Amat.final?
    )

    if (nM>1) { cat("\n") }
  }

  out <- list(
    model_results = out,
    contrasts = contrast_list,
    excursion_type = excursion_type,
    field_names = field_names,
    session_names = session_names,
    gamma = gamma,

    # Single-threshold mode: alpha is kept, alpha_nested / alpha_grid are NULL
    # Nested-threshold mode: alpha is NULL, alpha_nested / alpha_grid are kept
    alpha = alpha,
    alpha_nested = alpha_nested,
    alpha_grid = alpha_grid,
    activation_mode = if (use_nested_alpha) "nested" else "single",

    nsamp_theta = nsamp_theta,
    nsamp_beta = nsamp_beta
  )
  class(out) <- "fit_bglm2"

  if (is_cifti) {

    # Set values in maskIn but not maskMdat to `NA`.
    # Mask with maskIn.
    result_oomSetNA <- out$model_results
    for (mm in seq(nM)) {
      spatial_type <- spatial_type_by_model[mm]
      spatial_sub <- spatial_sub_by_model[[mm]]

      if (spatial_type == "vertex") {
        result_oomSetNA[[mm]]$estimates[Masks$In[[mm]] & (!Masks$Mdat[[mm]]),] <- NA
        result_oomSetNA[[mm]]$estimates <- result_oomSetNA[[mm]]$estimates[Masks$In[[mm]],,drop=FALSE]

        if (!is.null(result_oomSetNA[[mm]]$ppm)) {
          result_oomSetNA[[mm]]$ppm[Masks$In[[mm]] & (!Masks$Mdat[[mm]]),] <- NA
          result_oomSetNA[[mm]]$ppm <- result_oomSetNA[[mm]]$ppm[Masks$In[[mm]],,drop=FALSE]
        }
        if (!is.null(result_oomSetNA[[mm]]$active)) {
          result_oomSetNA[[mm]]$active[Masks$In[[mm]] & (!Masks$Mdat[[mm]]),] <- NA
          result_oomSetNA[[mm]]$active <- result_oomSetNA[[mm]]$active[Masks$In[[mm]],,drop=FALSE]
        }

      } else {
        result_oomSetNA[[mm]]$estimates <- unmask_Mdat2In(
          result_oomSetNA[[mm]]$estimates[spatial_sub$Mmap,,drop=FALSE],
          spatial_sub$maskIn[],
          spatial_sub$maskMdat[]
        )

        if (!is.null(result_oomSetNA[[mm]]$ppm)) {
          result_oomSetNA[[mm]]$ppm <- unmask_Mdat2In(
            result_oomSetNA[[mm]]$ppm[spatial_sub$Mmap,,drop=FALSE],
            spatial_sub$maskIn[],
            spatial_sub$maskMdat[]
          )
        }
        if (!is.null(result_oomSetNA[[mm]]$active)) {
          result_oomSetNA[[mm]]$active <- unmask_Mdat2In(
            result_oomSetNA[[mm]]$active[spatial_sub$Mmap,,drop=FALSE],
            spatial_sub$maskIn[],
            spatial_sub$maskMdat[]
          )
        }
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
      nested_activations_xii = NULL,
      activation_levels_xii = NULL,
      masks = Masks,
      BayesGLM2_results = out
    )
    out$contrast_estimate_xii$meta$cifti$names <- names(contrast_list)

    if (do_excur) {

      result_oomSetNA <- out$BayesGLM2_results$model_results

      for (mm in seq(nM)) {
        spatial_type <- spatial_type_by_model[mm]
        spatial_sub <- spatial_sub_by_model[[mm]]

        if (!use_nested_alpha) {

          ## Single-threshold mode: unwrap `active`
          if (spatial_type == "vertex") {
            if (!is.null(result_oomSetNA[[mm]]$active)) {
              result_oomSetNA[[mm]]$active[Masks$In[[mm]] & (!Masks$Mdat[[mm]]), ] <- NA
              result_oomSetNA[[mm]]$active <- result_oomSetNA[[mm]]$active[Masks$In[[mm]], , drop = FALSE]
            }

          } else {

            if (!is.null(result_oomSetNA[[mm]]$active)) {
              result_oomSetNA[[mm]]$active <- unmask_Mdat2In(
                result_oomSetNA[[mm]]$active[spatial_sub$Mmap, , drop = FALSE],
                spatial_sub$maskIn[],
                spatial_sub$maskMdat[]
              )
            }
          }

        } else {

          ## Nested-threshold mode: unwrap `nested_code`
          if (spatial_type == "vertex") {
            if (!is.null(result_oomSetNA[[mm]]$nested_code)) {
              result_oomSetNA[[mm]]$nested_code[Masks$In[[mm]] & (!Masks$Mdat[[mm]]), ] <- NA
              result_oomSetNA[[mm]]$nested_code <- result_oomSetNA[[mm]]$nested_code[Masks$In[[mm]], , drop = FALSE]
            }

          } else {

            if (!is.null(result_oomSetNA[[mm]]$nested_code)) {
              result_oomSetNA[[mm]]$nested_code <- unmask_Mdat2In(
                result_oomSetNA[[mm]]$nested_code[spatial_sub$Mmap, , drop = FALSE],
                spatial_sub$maskIn[],
                spatial_sub$maskMdat[]
              )
            }
          }
        }
      }

      if (!use_nested_alpha) {

        ## Build the original binary activation dlabel
        act_xii <- as.xifti(
          cortexL = result_oomSetNA$cortexL$active,
          cortexL_mwall = Masks$In$cortexL,
          cortexR = result_oomSetNA$cortexR$active,
          cortexR_mwall = Masks$In$cortexR,
          c(NA, NaN),
          subcortVol = result_oomSetNA$subcort$active,
          subcortLabs = spatial_sub$labels,
          subcortMask = spatial_sub$maskIn
        )

        out$activations_xii <- convert_xifti(act_xii, "dlabel", colors = "red")
        out$activations_xii$meta$cifti$names <- names(contrast_list)
        names(out$activations_xii$meta$cifti$labels) <- names(contrast_list)

      } else {

        ## Build nested activation dlabel
        nested_xii <- as.xifti(
          cortexL = result_oomSetNA$cortexL$nested_code,
          cortexL_mwall = Masks$In$cortexL,
          cortexR = result_oomSetNA$cortexR$nested_code,
          cortexR_mwall = Masks$In$cortexR,
          c(NA, NaN),
          subcortVol = result_oomSetNA$subcort$nested_code,
          subcortLabs = spatial_sub$labels,
          subcortMask = spatial_sub$maskIn
        )

        ## Basic color palette:
        ## level 0 = not active
        ## higher levels = more stringent nested significance
        max_nested_level <- max(vapply(alpha_grid, length, integer(1)))

        nested_colors <- grDevices::colorRampPalette(
          c("gold", "darkorange", "red3", "darkred")
        )(max_nested_level)

        out$nested_activations_xii <- convert_xifti(
          nested_xii,
          "dlabel",
          colors = nested_colors
        )

        out$nested_activations_xii$meta$cifti$names <- names(contrast_list)
        names(out$nested_activations_xii$meta$cifti$labels) <- names(contrast_list)

        ## Build per-contrast exact activation-level xifti objects
        out$activation_levels_xii <- vector("list", nC)
        names(out$activation_levels_xii) <- names(contrast_list)

        for (cc in seq_len(nC)) {

          result_level <- out$BayesGLM2_results$model_results
          nlev <- length(alpha_grid[[cc]])

          for (mm in seq(nM)) {
            spatial_type <- spatial_type_by_model[mm]
            spatial_sub <- spatial_sub_by_model[[mm]]

            ## Use nested_code to build exact level maps:
            ## level j = locations with nested_code exactly equal to j
            code_vec <- result_level[[mm]]$nested_code[, cc]

            lev_mat <- vapply(
              seq_len(nlev),
              function(j) as.integer(code_vec == j),
              integer(length(code_vec))
            )

            if (is.null(dim(lev_mat))) {
              lev_mat <- matrix(lev_mat, ncol = 1)
            }

            if (spatial_type == "vertex") {
              lev_mat[Masks$In[[mm]] & (!Masks$Mdat[[mm]]), ] <- NA
              lev_mat <- lev_mat[Masks$In[[mm]], , drop = FALSE]
            } else {
              lev_mat <- unmask_Mdat2In(
                lev_mat[spatial_sub$Mmap, , drop = FALSE],
                spatial_sub$maskIn[],
                spatial_sub$maskMdat[]
              )
            }

            result_level[[mm]]$activation_levels_exact <- lev_mat
          }

          level_xii <- as.xifti(
            cortexL = result_level$cortexL$activation_levels_exact,
            cortexL_mwall = Masks$In$cortexL,
            cortexR = result_level$cortexR$activation_levels_exact,
            cortexR_mwall = Masks$In$cortexR,
            c(NA, NaN),
            subcortVol = result_level$subcort$activation_levels_exact,
            subcortLabs = spatial_sub$labels,
            subcortMask = spatial_sub$maskIn
          )

          level_xii <- convert_xifti(level_xii, "dlabel", colors = "red")
          level_xii$meta$cifti$names <- paste0(
            names(contrast_list)[cc],
            "_level_",
            seq_len(nlev)
          )
          names(level_xii$meta$cifti$labels) <- level_xii$meta$cifti$names

          out$activation_levels_xii[[cc]] <- level_xii
        }
      }
    }
    class(out) <- "BGLM2"
  }
  out
}

