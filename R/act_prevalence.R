#' Activations prevalence.
#' 
#' @param act_list List of activations from \code{\link{id_activations}}. All
#'  should have the same sessions, fields, and brainstructures.
#' @param gamma_idx If activtions at multiple thresholds were computed, which 
#'  threshold should be used for prevalence? Default: the first (lowest).
#'  
#' @return A list containing the prevalances of activation, as a proportion of
#'  the results from \code{act_list}.
#' 
#' @importFrom stats setNames
#' @importFrom ciftiTools convert_xifti
#' 
#' @export 
act_prevalance <- function(act_list, gamma_idx=1){

  # Determine if `act_BGLM` or `act_BGLM0`.
  is_cifti <- all(vapply(act_list, function(q){ inherits(q, "act_BGLM") }, FALSE))
  if (!is_cifti) {
    if (!all(vapply(act_list, function(q){ inherits(q, "act_BGLM0") }, FALSE))) {
      stop("All objects in `act_list` must be the same type of result from `id_activations`: either `act_BGLM` or `act_BGLM0`.")
    }
  }

  # Get the number of results, sessions, and fields, and brainstructures (for CIFTI).
  # Ensure sessions, fields, and brainstructures match for all results.
  # [TO DO] could check that the number of locations is also the same.
  #   but maybe not here because that's more complicated for CIFTI.
  nA <- length(act_list)
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
  for (aa in seq(2, nA)) {
    if (length(act_list[[aa]]$session_names) != nS) {
      stop("Result ", aa, " has a different number of sessions than the first result.")
    }
    if (!all(act_list[[aa]]$session_names == session_names)) {
      warning("Result ", aa, " has different session names than the first result.")
    }
    if (length(act_list[[aa]]$field_names) != nK) {
      stop("Result ", aa, " has a different number of fields than the first result.")
    }
    if (!all(act_list[[aa]]$field_names == field_names)) {
      warning("Result ", aa, " has different field names than the first result.")
    }
    if (is_cifti) {
      if (length(act_list[[aa]]$activations) != nB) {
        stop("Result ", aa, " has a different number of brain structures than the first result.")
      }
      if (!all(names(act_list[[aa]]$activations) == bs_names)) {
        warning("Result ", aa, " has different brain structure names than the first result.")
      }
    }
  }

  # Compute prevalance, for every session and every field.
  prev <- setNames(rep(list(setNames(vector("list", nS), session_names)), nB), bs_names)
  for (bb in seq(nB)) {
    for (ss in seq(nS)) {
      x <- lapply(act_list, function(y){
        y <- if (is_cifti) { y$activations[[bb]] } else { y$activations }
        y[[ss]][[gamma_idx]]$active
      })
      prev[[bb]][[ss]] <- Reduce("+", x)/nA
    }
  }

  if (!is_cifti) { prev <- prev[[1]] }

  result <- list(
    prevalence = prev,
    n_results = nA,
    field_names = field_names,
    session_names = session_names
  )

  # If BGLM0, return.
  if (!is_cifti) {
    class(result) <- "prev_BGLM0"
    return(result)
  }

  # If BGLM, create 'xifti' with activations.
  prev_xii <- vector("list", nS)
  names(prev_xii) <- session_names
  for (session in session_names) {
    prev_xii_ss <- 0*convert_xifti(act_list[[1]]$activations_xii[[session]], "dscalar")
    prev_xii_ss$meta$cifti$names <- field_names
    for (bs in names(prev_xii_ss$data)) {
      if (!is.null(prev_xii_ss$data[[bs]])) {
        dat <- prev[[bs]][[session]]
        colnames(dat) <- NULL
        prev_xii_ss$data[[bs]] <- dat
      }
    }
    prev_xii[[session]] <- prev_xii_ss
  }

  result <- c(list(prev_xii=prev_xii), result)
  class(result) <- "prev_BGLM"
  result
}
