#' Activations prevalence.
#' 
#' @param act_list List of activations from \code{\link{id_activations}}. All
#'  should have the same sessions, fields, and brainstructures.
#'  
#' @return A list containing the prevalances of activation, as a proportion of
#'  the results from \code{act_list}.
#' 
#' @importFrom stats setNames
#' @importFrom ciftiTools convert_xifti
#' 
#' @export 
act_prevalance <- function(act_list){

  # Determine if `act_BayesGLM_cifti` or `act_BayesGLM`.
  is_cifti <- all(vapply(act_list, function(q){ inherits(q, "act_BayesGLM_cifti") }, FALSE))
  if (!is_cifti) {
    if (!all(vapply(act_list, function(q){ inherits(q, "act_BayesGLM") }, FALSE))) {
      stop("All objects in `act_list` must be the same type of result from `id_activations`: either `act_BayesGLM_cifti` or `act_BayesGLM`.")
    }
  }

  # Get the number of results, sessions, and fields, and brainstructures (for CIFTI).
  # Ensure sessions, fields, and brainstructures match for all results.
  # [TO DO] could check that the number of locations is also the same.
  #   but maybe not here because that's more complicated for CIFTI.
  nA <- length(act_list)
  session_names <- act_list[[1]]$session_names
  nS <- length(session_names)
  task_names <- act_list[[1]]$task_names
  nK <- length(task_names)
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
    if (length(act_list[[aa]]$task_names) != nK) {
      stop("Result ", aa, " has a different number of tasks than the first result.")
    }
    if (!all(act_list[[aa]]$task_names == task_names)) {
      warning("Result ", aa, " has different task names than the first result.")
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

  # Compute prevalance, for every session and every task.
  prev <- setNames(rep(list(setNames(vector("list", nS), session_names)), nB), bs_names)
  for (bb in seq(nB)) {
    for (ss in seq(nS)) {
      x <- lapply(act_list, function(y){
        y <- if (is_cifti) { y$activations[[bb]] } else { y$activations }
        y[[ss]]$active
      })
      prev[[bb]][[ss]] <- Reduce("+", x)/nA
    }
  }

  if (!is_cifti) { prev <- prev[[1]] }

  result <- list(
    prevalence = prev,
    n_results = nA,
    task_names = task_names,
    session_names = session_names
  )

  # If BayesGLM, return.
  if (!is_cifti) {
    class(result) <- "prev_BayesGLM"
    return(result)
  }

  # If BayesGLM_cifti, create 'xifti' with activations.
  prev_xii <- vector("list", nS)
  names(prev_xii) <- session_names
  for (session in session_names) {
    prev_xii_ss <- 0*convert_xifti(act_list[[1]]$activations_xii[[session]], "dscalar")
    prev_xii_ss$meta$cifti$names <- task_names
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
  class(result) <- "prev_BayesGLM_cifti"
  result
}
