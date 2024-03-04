#' BOLD
#'
#' @param BOLD fMRI timeseries data in CIFTI format ("*.dtseries.nii").
#'  For single-session analysis this can be a file path to a CIFTI file or a
#'  \code{"xifti"} object from the \code{ciftiTools} package. For multi-session
#'  analysis this can be a vector of file paths or a list of \code{"xifti"}
#'  objects.
#'
#'  If \code{BOLD} is a \code{"xifti"} object(s), the surfaces, if any, will be
#'  used for the spatial model. However, if \code{surfL} and \code{surfR} are 
#'  provided, they will override any surfaces in \code{BOLD}.
#'
#' @name BOLD_Param_BayesGLM_cifti
NULL

#' design
#' 
#' @param design A numeric matrix or \code{data.frame}, or a 
#'  \code{"BayesfMRI_design"} object from \code{\link{make_design}}.
#' 
#' @name design_Param_BayesGLM_cifti
NULL

#' TR
#' 
#' @param TR Temporal resolution of the data, in seconds.
#' 
#' @name TR_Param_BayesGLM_cifti 
NULL

#' brainstructures
#'
#' @param brainstructures Character vector indicating which brain structure(s)
#'  of \code{BOLD} to analyze: \code{"left"} cortex; \code{"right"} cortex;
#'  and/or \code{"subcortical"} structures. Default: \code{c("left","right")} 
#'  (cortex only).
#'
#' @name brainstructures_Param_BayesGLM_cifti
NULL

#' surfaces 
#' 
#' @param surfL,surfR For cortex spatial model. Left and right cortex surface 
#'  geometry in GIFTI format ("*.surf.gii"). These can be a file path to
#'  a GIFTI file or a \code{"surf"} object from \code{ciftiTools}. 
#' 
#'  Surfaces can alternatively be provided through the \code{$surf} metadata in
#'  \code{BOLD} if it is \code{"xifti"} data. If neither are provided, the
#'  HCP group-average inflated surface included in \code{ciftiTools} will be
#'  used for the cortex spatial model. 
#   [TO DO: does this work?]
#' 
#' @name surfaces_Param_BayesGLM_cifti
NULL

#' resamp_res
#' 
#' @param resamp_res For cortex spatial model. The number of vertices to which
#'  each cortical surface should be resampled, or \code{NULL} to not resample. 
#'  
#'  For computational feasibility, a value of \code{10000} (default) or lower is
#'  recommended for Bayesian spatial modeling. If \code{Bayes=FALSE}, 
#'  \code{resamp_res} can be set to \code{NULL} for full-resolution classical
#'  modeling. 
#' 
#' @name resamp_res_Param_BayesGLM_cifti
NULL

#' nuisance 
#' 
#' @param nuisance (Optional) A \eqn{T \times N} matrix of nuisance signals,
#'  where \eqn{T} is the number of timepoints and \eqn{N} is the number of
#'  nuisance signals, or a list of these for multi-session analysis. Nuisance
#'  signals are regressed from the fMRI data and design matrix prior to GLM
#'  computation. Nuisance signals can include motion regressors, HRF derivatives
#'  not being modeled as tasks, and other sources of noise.
#' 
#'  Detrending/high-pass filtering is accomplished by adding DCT bases to the
#'  nuisance matrix; see the parameters \code{hpf} and \code{DCT}.
#' 
#' @name nuisance_Param_BayesGLM_cifti 
NULL

#' hpf
#' 
#' @param hpf Add DCT bases to \code{nuisance} to apply a temporal high-pass 
#'  filter to the data, for detrending? \code{hpf} is the filter frequency.
#'  Use \code{NULL} to skip detrending. Detrending is strongly recommended for
#'  fMRI data, to help reduce the autocorrelation in the residuals.
#'
#'  Using at least two DCT bases is as sufficient for detrending as using linear
#'  and quadratic drift terms in the nuisance matrix. So if DCT detrending is 
#'  being used here, there is no need to add linear and quadratic drift terms to
#'  \code{nuisance}.
#' 
#' @name hpf_Param_BayesGLM_cifti 
NULL