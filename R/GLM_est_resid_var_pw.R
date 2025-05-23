#' Standardize data variance, and prewhiten if applicable
#'
#' Standardize data variance and prewhiten if applicable, for the GLM.
#' @param BOLD,design,spatial See \code{fit_bayesglm}.
#' @param session_names,field_names,design_type See \code{fit_bayesglm}.
#' @param valid_cols,nT,do_pw See \code{fit_bayesglm}.
#' @return List of results
#' @keywords internal
GLM_est_resid_var_pw <- function(
  BOLD, design, spatial,
  session_names, field_names, design_type,
  valid_cols, nT,
  ar_order, ar_smooth, aic, n_threads,
  do_pw
){

  nS <- length(session_names)
  nK <- length(field_names)
  nV_D <- get_nV(spatial)$mdata

  AR_coefs_avg <- var_avg <- max_AIC <- NULL

  var_resid <- array(dim = c(nV_D, nS))
  if (do_pw) {
    AR_coefs <- array(dim = c(nV_D, ar_order, nS))
    AR_AIC <- if (aic) { array(dim = c(nV_D, nS)) } else { NULL }
  }

  # Estimate parameters for each session.
  for (ss in seq(nS)) {
    vcols_ss <- valid_cols[ss,]

    if (design_type == "regular") {
      resid_ss <- fMRItools::nuisance_regression(
        BOLD[[ss]], design[[ss]][,vcols_ss]
      )
    } else if (design_type == "per_location") {
      resid_ss <- BOLD[[ss]]*0
      for (dd in seq(nV_D)) {
        resid_ss[,dd] <- fMRItools::nuisance_regression(
          BOLD[[ss]][,dd,drop=FALSE],
          matrix(design[[ss]][,vcols_ss,dd,drop=FALSE], nrow=dim(design[[ss]])[1])
        )
      }
    } else { stop() }

    if (do_pw) {
      pw_est_ss <- pw_estimate(resid_ss, ar_order, aic=aic)
      var_resid[,ss] <- pw_est_ss$sigma_sq
      AR_coefs[,,ss] <- pw_est_ss$phi
      if (aic) { AR_AIC[,ss] <- pw_est_ss$aic }
    } else {
      var_resid[,ss] <- matrixStats::colVars(resid_ss)
    }
  }

  # Average across sessions.
  var_avg <- rowMeans(as.matrix(var_resid))
  if (!do_pw) { var_avg <- var_avg/mean(var_avg, na.rm=TRUE) }
  if (do_pw) {
    AR_coefs_avg <- apply(AR_coefs, seq(2), mean)
    if (aic) { max_AIC <- apply(AR_AIC, 1, max) }
  }

  # Smooth prewhitening parameters.
  if (do_pw && ar_smooth > 0) {
    x <- pw_smooth(
      spatial=spatial,
      AR=AR_coefs_avg, var=var_avg,
      FWHM=ar_smooth
    )
    AR_coefs_avg <- x$AR
    var_avg <- x$var
    rm(x)
  }

  sqrtInv_all <- lapply(nT, function(q){ make_sqrtInv_all(q,
    nV_D, do_pw, n_threads, ar_order, AR_coefs_avg, var_avg
  )})

  list(
    var_resid=var_resid, sqrtInv_all=sqrtInv_all,
    AR_coefs_avg=AR_coefs_avg, var_avg=var_avg, max_AIC=max_AIC
  )
}

#' Make \code{sqrtInv_all}
#'
#' Make \code{sqrtInv_all} for prewhitening
#' @param nT,nV,do_pw,n_threads,ar_order,AR_coefs_avg,var_avg See \code{\link{GLM_est_resid_var_pw}}.
#' @return \code{sqrtInv_all}
#' @keywords internal
make_sqrtInv_all <- function(
  nT, nV, do_pw, n_threads, ar_order, AR_coefs_avg, var_avg){

  if (do_pw) {
    # Case 1A: Prewhitening; not parallel.
    if (is.null(n_threads) | n_threads < 2) {
      # Initialize the block diagonal covariance matrix
      template_pw <- Matrix::bandSparse(
        n = nT, k = seq(0, ar_order+1), symmetric = TRUE
      )
      template_pw_list <- rep(list(template_pw), nV)
      for (vv in seq(nV)) {
        #if (vv %% 100 == 0) if (verbose>1) { cat("\tLocation",vv,"of",nV,"\n") }
        template_pw_list[[vv]] <- .getSqrtInvCpp(
          AR_coefs = AR_coefs_avg[vv,],
          nTime = nT,
          avg_var = var_avg[vv]
        )
      }

    # Case 1B: Prewhitening; not parallel.
    } else {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Prewhitening in parallel requires the `parallel` package. Please install it.", call. = FALSE)
      }
      num_threads <- min(max(parallel::detectCores(), 25), n_threads)
      cl <- parallel::makeCluster(num_threads)
      template_pw_list <- parallel::clusterMap(
        cl,
        .getSqrtInvCpp,
        AR_coefs = split(AR_coefs_avg, row(AR_coefs_avg)),
        nTime = nT,
        avg_var = var_avg,
        SIMPLIFY = FALSE
      )
      parallel::stopCluster(cl)
    }

    #consider using a variant of bdiag_m if this is very slow.  See help(Matrix::bdiag)
    sqrtInv_all <- Matrix::bdiag(template_pw_list)

  # Case 2: Prewhitening; not parallel.
  } else if (!do_pw) {
    diag_values <- rep(1/sqrt(var_avg), each = nT)
    sqrtInv_all <- Diagonal(x = rep(1/sqrt(var_avg), each = nT))
  } else { stop() }

  sqrtInv_all
}
