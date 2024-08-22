#' Activations prevalence.
#'
#' @param act_list List of activations from \code{\link{activations}}. All
#'  should have the same sessions, fields, and brainstructures.
#' @param gamma_idx If activations at multiple thresholds were computed, which
#'  threshold should be used for prevalence? Default: the first (lowest).
#' @param p_test For inference: the expected baseline rate of activation for all
#'  data locations, under the null hypothesis. Default: \code{NULL} (do not
#'  perform hypothesis testing).
#' @param alpha Significance level for inference. Default: \code{.05}.
#' @param correction For the classical method only: Type of multiple comparisons
#'  correction: \code{"FWER"} (Bonferroni correction, the default), \code{"FDR"}
#'  (Benjamini Hochberg), or \code{"none"}.
#'
#' @return A list containing the prevalences of activation, as a proportion of
#'  the results from \code{act_list}.
#'
#' @importFrom stats setNames pbinom
#' @importFrom fMRItools is_1
#' @importFrom ciftiTools convert_xifti
#'
#' @export
prevalence <- function(
  act_list, gamma_idx=1,
  p_test=NULL, alpha=.05, correction=c("FWER", "FDR", "none")){

  # Determine if `act_BGLM` or `act_fit_bglm`.
  is_cifti <- all(vapply(act_list, function(q){ inherits(q, "act_BGLM") }, FALSE))
  if (!is_cifti) {
    if (!all(vapply(act_list, function(q){ inherits(q, "act_fit_bglm") }, FALSE))) {
      stop("All objects in `act_list` must be the same type of result from `activations`: either `act_BGLM` or `act_fit_bglm`.")
    }
  }

  # Get the number of results, sessions, and fields, and brainstructures (for CIFTI).
  # Ensure sessions, fields, and brainstructures match for all results.
  # [TO DO] could check that the number of locations is also the same.
  #   but maybe not here because that's more complicated for CIFTI.
  nN <- length(act_list)
  session_names <- act_list[[1]]$session_names
  nS <- length(session_names)
  field_names <- act_list[[1]]$field_names
  nK <- length(field_names)
  if (is_cifti) {
    bs_names <- names(act_list[[1]]$activations)
  } else {
    bs_names <- "activations"
  }
  nB <- length(bs_names)
  for (nn in seq(2, nN)) {
    if (length(act_list[[nn]]$session_names) != nS) {
      stop("Result ", nn, " has a different number of sessions than the first result.")
    }
    if (!all(act_list[[nn]]$session_names == session_names)) {
      warning("Result ", nn, " has different session names than the first result.")
    }
    if (length(act_list[[nn]]$field_names) != nK) {
      stop("Result ", nn, " has a different number of fields than the first result.")
    }
    if (!all(act_list[[nn]]$field_names == field_names)) {
      warning("Result ", nn, " has different field names than the first result.")
    }
    if (is_cifti) {
      if (length(act_list[[nn]]$activations) != nB) {
        stop("Result ", nn, " has a different number of brain structures than the first result.")
      }
      if (!all(names(act_list[[nn]]$activations) == bs_names)) {
        warning("Result ", nn, " has different brain structure names than the first result.")
      }
    }
  }

  # Get and apply intersection mask
  Masks <- intersect_mask(act_list)
  for (nn in seq(nN)) {
    cat(paste0("Checking data mask for subject ", nn, ".\n")) # [TO DO] option to hide?
    act_list[[nn]] <- retro_mask_act(act_list[[nn]], Masks)
  }

  # Compute prevalence, for every session and every field.
  prev <- setNames(rep(list(setNames(vector("list", nS), session_names)), nB), bs_names)
  if (!is.null(p_test)) { prev_pval <- prev_test <- prev }
  for (bb in seq(nB)) {
    for (ss in seq(nS)) {

      # Prevalence
      x <- lapply(act_list, function(y){
        y <- if (is_cifti) { y$activations[[bb]] } else { y$activations }
        y[[ss]][[gamma_idx]]$active
      })
      n_act <- Reduce("+", x)
      prev[[bb]][[ss]] <- n_act/nN

      # Inference
      if (!is.null(p_test)) {
        # Arg checks
        stopifnot(fMRItools::is_1(p_test, "numeric"))
        stopifnot(p_test >= 0)
        stopifnot(p_test <= 1)
        correction <- match.arg(correction, c("FWER", "FDR", "none"))

        # Compute pvals
        n_act <- as.matrix(n_act)
        pvals <- pbinom(n_act[], nN, p_test, lower.tail = FALSE, log.p = FALSE)

        # Multiple testing correction
        pvals_adj <- pvals * 0
        for (kk in seq(ncol(pvals))) {
          pvals_adj[,kk] <- switch(correction,
            FWER = p.adjust(pvals[,kk], method='bonferroni'),
            FDR = p.adjust(pvals[,kk], method='BH'),
            none = pvals[,kk]
          )
        }

        prev_pval[[bb]][[ss]] <- pvals_adj
        prev_test[[bb]][[ss]] <- (pvals_adj < alpha) * 1
      }
    }
  }

  if (!is_cifti) {
    prev <- prev[[1]]
    if (!is.null(p_test)) {
      prev_pval <- prev_pval[[1]]
      prev_test <- prev_test[[1]]
    }
  }

  result <- list(
    prevalence = prev,
    n_results = nN,
    field_names = field_names,
    session_names = session_names
  )
  if (!is.null(p_test)) {
    result$prev_pval <- prev_pval
    result$prev_test <- prev_test
  }

  # If fit_bglm, return.
  if (!is_cifti) {
    class(result) <- "prev_fit_bglm"
    return(result)
  }

  # If BGLM, create 'xifti' with activations.
  prev_xii <- vector("list", nS)
  names(prev_xii) <- session_names
  if (!is.null(p_test)) { prev_test_xii <- prev_xii }
  for (session in session_names) {
    prev_xii_ss <- 0*convert_xifti(act_list[[1]]$activations_xii[[session]], "dscalar")
    prev_xii_ss$meta$cifti$names <- field_names
    if (!is.null(p_test)) { prev_test_xii_ss <- prev_xii_ss }
    for (bs in names(prev_xii_ss$data)) {
      bs2 <- switch(bs,
        cortex_left="cortexL",
        cortex_right="cortexR",
        subcort="subcort"
      )
      if (!is.null(prev_xii_ss$data[[bs]])) {
        dat <- prev[[bs2]][[session]]
        colnames(dat) <- NULL
        prev_xii_ss$data[[bs]] <- dat
        if (!is.null(p_test)) {
          dat <- prev_test[[bs2]][[session]]
          colnames(dat) <- NULL
          prev_test_xii_ss$data[[bs]] <- dat
        }
      }
    }
    prev_xii[[session]] <- prev_xii_ss
    if (!is.null(p_test)) {
      prev_test_xii[[session]] <- convert_xifti(
        prev_test_xii_ss,
        to = "dlabel",
        levels=c(0, 1),
        labels=c("Not Significant", "Significant"),
        colors="red"
      )
    }
  }

  q <- if (is.null(p_test)) { NULL } else { list(prev_test_xii=prev_test_xii) }
  result <- c(list(prev_xii=prev_xii), q, result)
  class(result) <- "prev_BGLM"

  result
}
