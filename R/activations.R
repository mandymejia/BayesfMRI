#' Identify field activations
#'
#' Identify areas of activation for each field from the result of \code{BayesGLM}
#'  or \code{fit_bayesglm}.
#'
#' @param x Result of \code{BayesGLM} or \code{fit_bayesglm} model
#'  call, of class \code{"BGLM"} or \code{"fit_bglm"}.
#' @param Bayes Use spatial Bayesian modeling to identify activations based on
#' the joint posterior distribution? Default: \code{TRUE}. If \code{FALSE},
#' activations will be based on classical (massive univariate) GLM model, with
#' multiple comparisons correction (see \code{correction}). Note that \code{TRUE}
#' is only applicable if \code{x} includes Bayesian results (i.e.
#' \code{x <- BayesGLM(..., Bayes = TRUE)} was run.)
#' @param gamma Activation threshold, for example \code{1} for 1 percent
#'  signal change if \code{scale_BOLD=="mean"} during model estimation. Setting
#'  a \code{gamma} is required for the Bayesian method; \code{NULL}
#'  (default) will use a \code{gamma} of zero for the classical method.
#' @param alpha Significance level for inference. Default: \code{0.05}.
#' @param correction For the classical method only: Type of multiple comparisons
#'  correction: \code{"FWER"} (Bonferroni correction, the default), \code{"FDR"}
#'  (Benjamini Hochberg), or \code{"none"}.
#' @param fields The field(s) to identify activations for. Give either the name(s)
#'  as a character vector, or the numerical indices. If \code{NULL} (default),
#'  analyze all fields.
#' @param sessions The session(s) to identify activations for. Give either the
#'  name(s) as a character vector, or the numerical indices. If \code{NULL}
#' (default), analyze the first session.
#' @inheritParams verbose_Param
#'
# @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @importFrom ciftiTools convert_xifti
#' @importFrom fMRItools is_posNum is_1
#' @importFrom utils packageVersion
#' @importFrom viridisLite plasma
#'
#' @return An \code{"act_BGLM"} or \code{"act_fit_bglm"} object, a
#'  list which indicates the activated locations along with related information.
#' @export
#'
activations <- function(
  x,
  Bayes=TRUE,
  gamma=NULL,
  alpha=0.05,
  correction = c("FWER", "FDR", "none"),
  fields=NULL,
  sessions=NULL,
  verbose = 1){

  # If 'BGLM', we will loop over the brain structures.
  is_cifti <- inherits(x, "BGLM")
  if (is_cifti) {
    cifti_obj <- x
    x <- cifti_obj$BGLMs
  } else {
    if (!inherits(x, "fit_bglm")) {
      stop("`x` is not a `'BGLM'` or 'fit_bglm' object.")
    }
    x <- list(bglm=x)
  }
  spatial <- lapply(x, '[[', "spatial" )
  models <- names(x) #which brainstructures are in x
  idx1 <- min(which(!vapply(x, is.null, FALSE))) #first brainstructure in x

  # Argument checks.
  # Get `fields` and `sessions`.
  stopifnot(is.null(fields) || is.character(fields) || is.numeric(fields))
  if (is.character(fields)) { stopifnot(all(fields %in% x[[idx1]]$field_names)) }
  if (is.numeric(fields)) {
    stopifnot(all.equal(fields, round(fields))==TRUE)
    stopifnot(min(fields) >= 1)
    stopifnot(max(fields) <= length(x[[idx1]]$field_names))
    stopifnot(length(fields) == length(unique(fields)))
    fields <- x[[idx1]]$field_names[fields]
  }
  if (is.null(fields)) { fields <- x[[idx1]]$field_names }
  stopifnot(is.null(sessions) || is.character(sessions) || is.numeric(sessions))
  if (is.character(sessions)) { stopifnot(all(sessions %in% x[[idx1]]$session_names)) }
  if (is.numeric(sessions)) {
    stopifnot(all.equal(sessions, round(sessions))==TRUE)
    stopifnot(min(sessions) >= 1)
    stopifnot(max(sessions) <= length(x[[idx1]]$session_names))
    stopifnot(length(sessions) == length(unique(sessions)))
    sessions <- x[[idx1]]$session_names[sessions]
  }
  if (is.null(sessions)) { sessions <- x[[idx1]]$session_names[1] }
  if(!(Bayes %in% c(TRUE, FALSE))) stop("`Bayes` must be TRUE or FALSE.")
  stopifnot(is_posNum(alpha) && alpha < 1)
  stopifnot(is.null(gamma) || (is.numeric(gamma) && all(gamma >= 0)))
  correction <- match.arg(correction, c("FWER", "FDR", "none"))
  stopifnot(is_posNum(verbose, zero_ok=TRUE))

  # Check that Bayesian results are available, if requested.
  if (Bayes && is.null(x[[idx1]]$INLA_model_obj)) {
    stop("`Bayes = TRUE` but only classical model results are available in `x`.
         Set Bayes to FALSE or re-run `BayesGLM` with `Bayes = TRUE`")
  }

  # Check gamma was set, if computing Bayesian activations.
  if (is.null(gamma)) {
    if (Bayes) stop("Must specify `gamma` when `Bayes` is `TRUE`.")
    if (!Bayes) gamma <- 0
  } else {
    gamma <- sort(gamma)
  }

  # Get activation function arguments that won't change.
  if(Bayes) actFUN <- activations.posterior else actFUN <- activations.classical
  actArgs <- list(fields=fields, alpha=alpha)
  if (!Bayes) {
    actArgs <- c(actArgs, list(correction=correction))
  } else {
    correction <- "not applicable"
  }

  nS <- length(sessions)
  nB <- length(models) #number of brainstructures
  nG <- length(gamma) #number of gamma levels
  nK <- length(fields)

  # Loop over models, sessions, and gamma to compute activations.
  activations <- setNames(rep(list(setNames(vector("list", nS), sessions)), nB), models)

  # Loop over models (brain structures, if `is_cifti`).
  for (bb in 1:nB) {
    if (is.null(x[[bb]])) { next }
    if (Bayes && identical(attr(x[[bb]]$INLA_model_obj, "format"), "minimal")) {
      stop("Bayesian activations are not available because `return_INLA` was set to `'minimal'` in the `BayesGLM()` call. Request the classical activations or re-run `BayesGLM()`.")
    }
    if (Bayes && models[bb]!="bglm") {
      if (verbose>0) cat(paste0("Identifying Bayesian GLM activations in ",models[bb],'\n'))
    }

    # Loop over sessions.
    for (ss in 1:nS) {

      # Loop over gamma.
      q <- paste0("gamma=", gamma)
      q <- setNames(vector("list", nG), q)
      for (gg in seq(nG)) {
        q[[gg]] <- do.call(actFUN,
          c(actArgs, list(x=x[[bb]], session=sessions[ss], gamma=gamma[gg]))
        )
        # Apply mask to get data locs, removing boundary locs.
        if (Bayes & !is.null(spatial[[bb]]$maskMdat)) {
          q[[gg]]$active <- q[[gg]]$active[spatial[[bb]]$Mmap,,drop=FALSE]
        }
        # Both `actFUN`, posterior and classical, work on the `Mdat`
        #   locations. So here, we need to apply the mask to get the data
        #   locations, removing boundary locations.
        if (!Bayes) {
          q[[gg]]$p_values <- unmask_Mdat2In(q[[gg]]$p_values, spatial[[bb]]$maskIn, spatial[[bb]]$maskMdat)
          q[[gg]]$p_values_adj <- unmask_Mdat2In(q[[gg]]$p_values_adj, spatial[[bb]]$maskIn, spatial[[bb]]$maskMdat)
        }
        q[[gg]]$active <- unmask_Mdat2In(q[[gg]]$active, spatial[[bb]]$maskIn, spatial[[bb]]$maskMdat)
      }
      activations[[bb]][[ss]] <- q
    }
  }

  result <- list(
    activations = activations,
    Bayes=Bayes,
    alpha=alpha,
    gamma=gamma,
    correction=correction,
    spatial=spatial,
    field_names = fields,
    session_names = sessions
  )

  # If class(x) = fit_bglm
  if (!is_cifti) {
    result$activations <- result$activations[[1]]
    class(result) <- "act_fit_bglm"
    return(result)
  }

  act_xii <- vector("list", length(sessions))
  names(act_xii) <- sessions
  for (ss in 1:nS) {
    sess <- sessions[ss]
    the_xii <- cifti_obj$estimate_xii$classical[[sess]]
    act_xii_ss <- 0*select_xifti(the_xii, match(fields, the_xii$meta$cifti$names))
    for (bs in names(the_xii$data)) {
      bs2 <- switch(bs,
        cortex_left="cortexL",
        cortex_right="cortexR",
        subcort="subcort"
      )
      if (!is.null(the_xii$data[[bs]])) {
        dat <- Reduce("+", lapply(activations[[bs2]][[sess]], function(q){1*q$active}))
        act_xii_ss$data[[bs]] <- dat
      }
    }

    act_labs <- paste0("Active, gamma=", gamma)

    stopifnot(is.xifti(act_xii_ss))

    viridis_params <- list(
      n = nG,
      begin = c(1/2, 1/3, 0)[min(nG, 3)],
      end = c(1/2, 2/3, 1)[min(nG, 3)],
      direction=-1
    )

    if (utils::packageVersion("ciftiTools") < "0.13.1") { # approximate package version
      act_xii_ss <- convert_xifti(
        act_xii_ss, "dlabel",
        values=setNames(seq(0, nG), c("Inactive", act_labs)),
        colors=do.call(viridisLite::plasma, viridis_params)
      )
    } else {
      act_xii_ss <- convert_to_dlabel(
        act_xii_ss,
        levels_old=seq(0, nG),
        labels=c("Inactive", act_labs),
        colors=do.call(viridisLite::plasma, viridis_params)
      )
    }
    act_xii[[sess]] <- act_xii_ss
  }

  result <- c(list(activations_xii=act_xii), result)
  class(result) <- "act_BGLM"
  result
}

#' Identify activations using joint posterior probabilities
#'
#' Identifies areas of activation given an activation threshold and significance
#'  level using joint posterior probabilities
#'
#' For a given latent field, identifies locations that exceed a certain activation
#'  threshold (e.g. 1 percent signal change) at a given significance level, based on the joint
#'  posterior distribution of the latent field.
#'
#' @param x Result of \code{BayesGLM}, of class \code{"BGLM"}.
#' @param fields,session,alpha,gamma See \code{\link{activations}}.
#' @return A list with two elements: \code{active}, which gives a matrix of zeros
#' and ones of the same dimension as \code{x$field_estimates${session}},
#' and \code{excur_result}, an object of class \code{"excurobj"} (see \code{\link{excursions.inla}} for
#'  more information).
#'
#' @importFrom excursions excursions.inla
#'
#' @keywords internal
activations.posterior <- function(
  x,
  fields, session,
  alpha=0.05, gamma){

  stopifnot(inherits(x, "fit_bglm"))

  sess_ind <- which(x$session_names == session)
  n_vox <- x$spde$n.spde
  #indices of beta vector corresponding to specified session
  inds <- (1:n_vox) + (sess_ind-1)*n_vox

	#loop over latent fields
	excur <- vector('list', length=length(fields))
	act <- matrix(NA, nrow=n_vox, ncol=length(fields))
	colnames(act) <- fields
	for (field in fields) {

		#if(is.null(area.limit)){
		res.exc <- excursions.inla(
        x$INLA_model_obj,
        name=field, ind=inds, u=gamma, type='>', alpha=alpha, method="EB",
        verbose=FALSE
      )
		#} else {
		#	res.exc <- excursions.inla.no.spurious(x$INLA_model_obj, mesh=mesh, name=f, ind=inds, u=gamma, type='>', method=excur_method, alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
		#}
	  which_f <- which(fields==field)
		act[,which_f] <- res.exc$E[inds] == 1
    excur[[which_f]] <- res.exc
	}
	result <- list(active=act, excursions_result=excur)

  # [TO DO]: get the mesh
  # #compute size of activations
  # areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
  # areas_act <- apply(act, 2, function(x) sum(areas_all[x==1]))

  # result$areas_all <- areas_all
  # result$areas_act <- areas_act

	result
}

#' Identification of areas of activation in a General Linear Model using classical methods
#'
#' @param x A \code{BayesGLM} object
#' @param gamma,alpha,correction,fields,session See \code{\link{activations}}.
#' @param mesh (Optional) An \code{"inla.mesh"} object (see \code{\link{make_mesh}} for
#'  surface data). Only necessary for computing surface areas of identified activations.
#'
#' @return A matrix corresponding to the
#'   0-1 activation status for the model coefficients.
#'
#' @importFrom stats sd pt p.adjust
#' @importFrom matrixStats colVars
#'
#' @keywords internal
activations.classical <- function(x,
                                  gamma = 0,
                                  alpha = 0.05,
                                  correction = c("FWER", "FDR", "none"),
                                  fields,
                                  session,
                                  mesh = NULL) {

  # Argument checks ------------------------------------------------------------
  if (!inherits(x, "fit_bglm")) {
    stop(paste0(
      "The model object is of class ",
      paste0(class(x), collapse=", "),
      " but should be of class 'fit_bglm'."
    ))
  }
  # fields, session, alpha, gamma checked in `activations`
  correction <- match.arg(correction, c("FWER","FDR","none"))

  fields_idx <- x$field_names %in% fields
  stopifnot(any(fields_idx))
  mask_In2Mdat <- x$spatial$maskMdat[x$spatial$maskIn]
  beta_est <- x$result_classical[[session]]$estimates[mask_In2Mdat,fields_idx,drop=FALSE]
  se_beta <- x$result_classical[[session]]$SE_estimates[mask_In2Mdat,fields_idx,drop=FALSE]
  DOF <- x$result_classical[[session]]$DOF

  nV_input <- nrow(beta_est)
  #if(any(!(fields %in% 1:K))) stop(paste0('fields must be between 1 and the number of fields, ',K))
  nK <- length(fields)

  #Compute t-statistics and p-values
  t_star <- (beta_est - gamma) / se_beta
  if(!is.matrix(t_star)) t_star <- matrix(t_star, nrow=nV_input)
  #perform multiple comparisons correction
  p_values <- p_values_adj <- active <- matrix(NA, nV_input, nK)

  for (kk in seq(nK)) {
    p_values_k <- sapply(t_star[,kk], pt, df = DOF, lower.tail = F)
    p_vals_adj_k <- switch(correction,
      FWER = p.adjust(p_values_k, method='bonferroni'),
      FDR = p.adjust(p_values_k, method='BH'),
      none = p_values_k
    )
    p_values[,kk] <- p_values_k
    p_values_adj[,kk] <- p_vals_adj_k
    active[,kk] <- (p_vals_adj_k < alpha)
  }

  #compute size of activations
  if(!is.null(mesh)){
    mask <- x$result_classical[[session]]$mask
    if(sum(mask) != nrow(mesh$loc)) stop('Supplied mesh is not consistent with mask in x.')
    areas_all <- diag(INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 0, output = list("c0"))$c0) #area of each vertex
    areas_act <- apply(active[mask==1,,drop=FALSE], 2, function(x) sum(areas_all[x==1]))
  } else {
    areas_all <- areas_act <- NULL
  }

  na_pvalues <- which(is.na(p_values[,1]))
  if (length(na_pvalues) > 0) {
    p_values <- p_values[-na_pvalues,,drop=FALSE]
    p_values_adj <- p_values_adj[-na_pvalues,,drop=FALSE]
    active <- active[-na_pvalues,,drop=FALSE]
  }

  result <- list(
    p_values = p_values,
    p_values_adj = p_values_adj,
    active = active,
    correction = correction,
    alpha = alpha,
    gamma = gamma,
    areas_all = areas_all,
    areas_act = areas_act
  )
  result
}

#' @rdname activations
#' @export
id_activations <- function(
  x,
  Bayes=TRUE,
  gamma=NULL,
  alpha=0.05,
  correction = c("FWER", "FDR", "none"),
  fields=NULL,
  sessions=NULL,
  verbose = 1){

  activations(
    x=x,
    Bayes=Bayes,
    gamma=gamma,
    alpha=alpha,
    correction = correction,
    fields=fields,
    sessions=sessions,
    verbose = verbose
  )
}
