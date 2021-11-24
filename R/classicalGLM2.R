# This is a function to perform the group analysis for the classical GLM
#' Group-level Classical GLM
#'
#' @param results Either (1) a list of length M of objects of class
#'  BayesGLM_cifti with classical results, or (2) a character vector of length M
#'  of file names output from the BayesGLM_cifti function. M is the number of
#'  subjects.
#' @param brainstructure Character string indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface) or \code{"right"} (right
#'  cortical surface) or \code{"subcortical"}. Default: \code{c("left","right")}
#'  (entire cortical surface).
#' @param session_name (character) The name of the session that should be
#' examined. If \code{NULL} (default), the average across all sessions is used.
#' @param gamma (For inference only) List of vectors of activation thresholds
#'   for the excursion set (each element corresponding to one contrast).
#'   Remember that if a contrast is not specified, the average is found.
#' @param correction  Type of multiple comparisons correction, 'FDR'
#'   (Benjamini Hochberg) or 'FWER' (Bonferroni correction)
#' @param alpha Significance level (e.g. 0.05)
#'
#' @importFrom stats t.test p.adjust
#'
#' @return a list containing...
#' @export
classicalGLM2 <- function(results,
                          brainstructure = NULL,
                          session_name = NULL,
                          gamma = 0,
                          correction = c("FWER", "FDR"),
                          alpha = 0.05) {
  result_classes <- sapply(results, class)
  results_class <- unique(result_classes)
  if (length(results_class) > 1)
    stop("All elements in results should have class
                                     BayesGLM_cifti or character")
  correction <- match.arg(correction, c("FWER", "FDR"))
  p_method <- switch(correction, FWER = "bonferroni", FDR = "fdr")

  if (results_class == "character") {
    first_result <- readRDS(results[1])$GLMs_classical
  }
  if (results_class == "BayesGLM_cifti") {
    first_result <- results[[1]]$GLMs_classical
  }
  brainstructures <-
    names(first_result)[which(!sapply(first_result, is.null))]
  if (is.null(brainstructure)) {
    message(
      paste0(
        "No brainstructure specified. Using first non-null brainstructure in the object, which is ",
        brainstructures[1],
        ".\n"
      )
    )
    brainstructure <- brainstructures[1]
  }
  sess_names <- names(first_result[[brainstructure]])
  has_avg <- "avg" %in% sess_names
  if (is.null(session_name)) {
    if (has_avg) {
      message(
        "No session specified, but the session named avg, denoting an average across sessions, is present and will be used."
      )
      session_name <- "avg"
    } else {
      message(
        paste0(
          "No session specified, and avg is not present. The first available session, ",
          sess_names[1],
          ", will be used"
        )
      )
      session_name <- sess_names[1]
    }
  }
  if (results_class == "character") {
    results <-
      sapply(results, function(x)
        readRDS(x)$GLMs_classical[[brainstructure]][[session_name]]$estimates, simplify = "array")
  }
  if (results_class == "BayesGLM_cifti") {
    results <-
      sapply(results, function(x)
        x$GLMs_classical[[brainstructure]][[session_name]]$estimates, simplify = "array")
  }
  avg_estimate <- apply(results, 1:2, mean)

  gamma_result <- sapply(gamma, function(g) {
    p_vals <-
      apply(results, 1:2, function(x)
        stats::t.test(x, alternative = "greater", mu = g)$p.value)
    p_vals_adj <- apply(p_vals, 2, stats::p.adjust, method = p_method)
    active <- p_vals_adj < alpha
    return(list(
      p_vals = p_vals,
      p_vals_adj = p_vals_adj,
      active = active
    ))
  }, simplify = F)
  names(gamma_result) <- paste0("gamma = ", gamma)
  out <- list(avg_estimate = avg_estimate,
              active_result = gamma_result)
  return(out)
}
