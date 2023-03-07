#' INLA
#'
#' @section INLA Requirement:
#'  This function requires the \code{INLA} package, which is not a CRAN package.
#'  See \url{https://www.r-inla.org/download-install} for easy installation instructions.
#'
#' @name INLA_Description
NULL

#' aic
#'
#' @param aic Use the AIC to select AR model order between \code{0} and 
#'  \code{ar_order}? Default: \code{FALSE}.
#'
#' @name aic_Param
NULL

#' ar_order
#'
#' @param ar_order (numeric) Controls prewhitening. If greater than zero, this
#'  should be a number indicating the order of the autoregressive model to use
#'  for prewhitening. If zero, do not prewhiten. Default: \code{6}.
#'
#' @name ar_order_Param
NULL

#' ar_smooth
#'
#' @param ar_smooth (numeric) FWHM parameter for smoothing. Remember that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} or \code{NULL}
#'  to not do any smoothing. Default: \code{5}.
#'
#' @name ar_smooth_Param
NULL

#' avg_sessions
#'
#' @param avg_sessions Average estimates for betas over multiple
#'  sessions? Default: \code{TRUE}.
#'
#' @name avg_sessions_Param
NULL

#'  Bayes
#'
#' @param Bayes If \code{TRUE} (default), will fit a spatial Bayesian GLM in 
#'  addition to the classical GLM. Classical GLM results are always returned. 
#'
#' @name Bayes_Param
NULL

#' contrasts
#'
#' @param contrasts List of contrast vectors to be passed to \code{inla::inla}.
#'
#' @name contrasts_Param
NULL

#' EM
#'
#' @param EM (logical) Should the EM implementation of the Bayesian GLM be used?
#'  Default: \code{FALSE}. This method is still in development.
#'
#' @name EM_Param
NULL

#' emTol
#'
#' @param emTol The stopping tolerance for the EM algorithm. Default: 
#'  \code{1e-3}.
#'
#' @name emTol_Param
NULL

#'  faces
#'
#' @param faces An \eqn{F x 3} matrix, where each row contains the vertex
#'  indices for a given triangular face in the mesh. \eqn{F} is the number of
#'  faces in the mesh.
#'
#' @name faces_Param
NULL

#' mask: vertices
#'
#' @param mask  A length \eqn{V} logical vector indicating if each vertex is
#'  within the input mask.
#'
#' @name mask_Param_vertices
NULL

#' mesh: INLA only
#'
#' @param mesh An \code{"inla.mesh"} object (see \code{\link{make_mesh}} for
#'  surface data).
#'
#' @name mesh_Param_inla
NULL

#' mesh: either
#'
#' @param mesh An \code{"inla.mesh"} object (see \code{\link{make_mesh}} for
#'  surface data)
#  or \code{"BayesfMRI.spde"} object (see \code{\link{create_spde_vol3D}} for subcortical data).
#'
#' @name mesh_Param_either
NULL

#' max.threads
#'
#' @param max.threads The maximum number of threads to use in the inla-program
#'  for model estimation. \code{0} (default) will use the maximum number of
#'  threads allowed by the system.
#'
#' @name max.threads_Param
NULL

#' num.threads
#'
#' @param num.threads The maximum number of threads to use in the inla-program
#'  for model estimation. Default: \code{4}.
#'
#' @name num.threads_Param
NULL

#' outfile
#'
#' @param outfile File name where results will be written (for use by
#'  \code{BayesGLM2}).
#'
#' @name outfile_Param
NULL

#' return_INLA_result
#'
#' @param return_INLA_result Return the INLA model object? (It can be large.)
#'  Default: \code{TRUE}. Required for running \code{id_activations}
#'  after, but not for running BayesGLM_joint after to get
#'  posterior quantities of group means or contrasts.
#'
#' @name return_INLA_result_Param
NULL

#' scale_BOLD
#'
#' @param scale_BOLD Option for scaling the BOLD response.
#' 
#' 	If \code{"auto"} (default), will use mean scaling except if demeaned data
#' 	is detected, in which case sd scaling will be used instead.
#' 
#' 	\code{"mean"} scaling will scale the data to percent local signal change.
#' 
#' 	\code{"sd"} scaling will scale the data by local standard deviation.
#' 
#' 	\code{"none"} will only center the data, not scale it. 
#'
#' @name scale_BOLD_Param
NULL

#' scale_design
#'
#' @param scale_design Scale the design matrix by dividing each column by its
#'  maximum and then subtracting the mean? Default: \code{TRUE}.
#'
#' @name scale_design_Param
NULL

#' seed
#'
#' @param seed Random seed (optional). Default: \code{NULL}.
#'
#' @name seed_Param
NULL

#' session_names
#'
#' @param session_names (Optional, and only relevant for multi-session modeling)
#'  Names of each session. Default: \code{NULL}. In \code{\link{BayesGLM}} this
#'  argument will overwrite the names of the list entries in \code{data}, if
#'  both exist. 
#'
#' @name session_names_Param
NULL

#' task_names
#'
#' @param task_names (Optional) Names of tasks represented in design matrix.
#'
#' @name task_names_Param
NULL

#' trim_INLA
#'
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'  result be trimmed to only what is necessary to use \code{id_activations}? 
#'  Default: \code{TRUE}.
#'
#' @name trim_INLA_Param
NULL

#' verbose: INLA only
#'
#' @param verbose Should INLA be run in verbose mode? Default: \code{FALSE}.
#'
#' @name verbose_Param_inla
NULL

#' verbose: direct only, TRUE
#'
#' @param verbose Should occasional updates be printed? Default: \code{TRUE}.
#'
#' @name verbose_Param_direct_TRUE
NULL

#' verbose: direct only, FALSE
#'
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @name verbose_Param_direct_FALSE
NULL

#'  vertices
#'
#' @param vertices A \eqn{V x 3} matrix, where each row contains the Euclidean
#'  coordinates at which a given vertex in the mesh is located. \eqn{V} is the
#'  number of vertices in the mesh
#'
#' @name vertices_Param
NULL
