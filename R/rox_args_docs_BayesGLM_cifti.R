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
#' @name BOLD_param_BayesGLM_cifti
NULL

#' brainstructures
#'
#' @param brainstructures Character vector indicating which brain structure(s)
#'  of \code{BOLD} to analyze: \code{"left"} cortex; \code{"right"} cortex;
#'  and/or \code{"subcortical"} structures. Default: \code{c("left","right")} 
#'  (cortex only).
#'
#' @name brainstructures_param_BayesGLM_cifti
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
#' @name nuisance_param_BayesGLM_cifti 
NULL

#' detrending
#' 
#' @param TR Temporal resolution of the data, in seconds.
#' @param hpf,DCT Add DCT bases to \code{nuisance} to apply a temporal high-pass 
#'  filter to the data, for detrending? Only one of these should be provided:
#'  the filter frequency \code{hpf}, or the number of DCT bases \code{DCT}.
#'  If \code{hpf} is provided, \code{TR} must be provided too and the 
#'  corresponding number of DCT bases will be calculated. 
#'
#'  Default: \code{DCT=4}. For typical \code{TR} and length of \code{BOLD}, four
#'  DCT bases amounts to a lower frequency cutoff than the approximately .01 Hz
#'  used in most studies. We selected this default to err on the side of 
#'  retaining more low-frequency information, but recommend setting these 
#'  arguments to values most appropriate for the analysis at hand.
#'
#'  Using at least two DCT bases is as sufficient for detrending as using linear
#'  and quadratic drift terms in the design matrix. So if DCT detrending is 
#'  being used used, there is no need to add linear and quadratic drift terms to
#'  \code{nuisance}.
#' 
#' @name detrending_param_BayesGLM_cifti 
NULL