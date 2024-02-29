#' Connectome Workbench
#' 
#' @section Connectome Workbench Requirement:
#'  
#'  This function uses a system wrapper for the 'wb_command' executable. The
#'  user must first download and install the Connectome Workbench, available
#'  from https://www.humanconnectome.org/software/get-connectome-workbench .
#' 
#' @name Connectome_Workbench_Description
NULL

#' INLA
#'
#' @section INLA Requirement:
#'  This function requires the \code{INLA} package, which is not a CRAN package.
#'  See \url{https://www.r-inla.org/download-install} for easy installation instructions.
#'
#' @name INLA_Description
NULL

#' INLA Latent Fields
#' 
#' @section INLA Latent Fields Limit:
#'  INLA computation times increase greatly when the number of columns in the
#'  design matrix exceeds five: when there are more than five tasks, or more
#'  than three tasks each with a temporal derivative modeled as a field. In 
#'  cases like the latter, we recommend modeling the temporal derivatives as 
#'  nuisance signals using the option \code{dHRF_as="nuisance"}, rather than 
#'  modeling the temporal derivatives as fields.
#' 
#' @name INLA_Latent_Fields_Limit_Description
NULL

#' aic
#'
#' @param aic (For prewhitening) Use the Akaike information criterion (AIC) to
#'  select AR model orders between \code{0} and \code{ar_order}? Default: 
#'  \code{FALSE}.
#'
#' @name aic_Param
NULL

#' ar_order
#'
#' @param ar_order (For prewhitening) The order of the autoregressive (AR) model
#'  to use for prewhitening. If \code{0}, do not prewhiten. Default: \code{6}. 
#' 
#'  For multi-session modeling, note that a single AR model is used; its 
#'  coefficients will be the average estimate from each session.
#'
#' @name ar_order_Param
NULL

#' ar_smooth
#'
#' @param ar_smooth (For prewhitening) The FWHM parameter for spatially 
#'  smoothing the coefficient estimates for the AR model to use for 
#'  prewhitening. Recall that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} to not smooth
#'  the estimates. Default: \code{5}.
#'
# [TO DO] vol vs surf?
#' @name ar_smooth_Param
NULL

#'  Bayes
#'
#' @param Bayes Fir the spatial Bayesian GLM? Default: \code{TRUE}. If
#'  \code{FALSE}, only fit the classical GLM. (The classical GLM is always
#'  returned, whether \code{Bayes} is \code{TRUE} or \code{FALSE}.)
#'
#' @name Bayes_Param
NULL

#' buffer 
#' 
#' @param buffer For volumetric model. The number of extra voxel layers around 
#'  the bounding box. Set to \code{NULL} for no buffer. (We recommend not
#'  changing \code{buffer} unless you know what you're doing. Instead, to reduce 
#'  the number of boundary voxels, adjust \code{nbhd_order}).
#' 
#' @name buffer_Param 
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
#' @param faces An \eqn{F \times 3} matrix, where each row contains the vertex
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

#' mean and variance tolerance 
#' 
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for both.
#' 
#' @name mean_var_Tol_Param
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
#  or \code{"BayesfMRI.spde"} object (see \code{\link{make_spde_vol3D}} for subcortical data).
#'
#' @name mesh_Param_either
NULL

#' max_threads
#'
#' @param max_threads The maximum number of threads to use in the inla-program
#'  for model estimation. \code{0} (default) will use the maximum number of
#'  threads allowed by the system.
#'
#' @name max_threads_Param
NULL

#' nbhd_order 
#' 
#' @param nbhd_order For volumetric model. What order neighborhood around data
#' locations to keep? \code{0} for no neighbors, \code{1} for 1st-order 
#'  neighbors, \code{2} for 1st- and 2nd-order neighbors, etc. Smaller values 
#'  will provide greater computational efficiency at the cost of higher variance
#'  around the edge of the data.
#' 
#' @name nbhd_order_Param
#' 
NULL

#' n_threads
#'
#' @param n_threads The maximum number of threads to use for parallel
#'  computations: prewhitening parameter estimation, and the inla-program model
#'  estimation. Default: \code{4}. Note that parallel prewhitening requires the
#'  \code{parallel} package.
#'
#' @name n_threads_Param
NULL

#' return_INLA
#'
#' @param return_INLA Return the INLA model object? (It can be large.) Use 
#'  \code{"trimmed"} (default) returns the results sufficient for 
#'  \code{\link{id_activations}} and \code{\link{BayesGLM2}}; \code{"minimal"}
#'  returns enough for \code{\link{BayesGLM2}} but not 
#'  \code{\link{id_activations}}; \code{"full"} returns the full \code{inla} 
#'  output. 
#'
#' @name return_INLA_Param
NULL

#' scale_BOLD
#'
#' @param scale_BOLD Controls scaling the BOLD response at each location. 
#'  \describe{
#'    \item{"auto":}{   (default) Use \code{"mean"} scaling, except if 
#'      demeaned data is detected (any location's mean < 1), use \code{"sd"} 
#'      scaling.}
#'    \item{"mean":}{   Scale the data to percent local signal change.}
#'    \item{"sd":}{   Scale the data by local standard deviation.}
#'    \item{"none":}{   Center the data but do not scale it.}
#' }
#' 
#' @name scale_BOLD_Param
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
#'  Names of each session. Default: \code{NULL}. Will overwrite the names of the
#'  \code{BOLD} data, if both are provided. 
#'
#' @name session_names_Param
NULL

#' field_names
#'
#' @param field_names (Optional) Names of fields represented in design matrix.
#'
#' @name field_names_Param
NULL

#' trim_INLA
#'
#' @param trim_INLA (logical) should the \code{INLA_model_obj} within the
#'  result be trimmed to only what is necessary to use \code{id_activations}? 
#'  Default: \code{TRUE}.
#'
#' @name trim_INLA_Param
NULL

#' verbose
#'
#' @param verbose \code{1} (default) to print occasional updates during model 
#'  computation; \code{2} for occasional updates as well as running INLA in
#'  verbose mode (if \code{Bayes}), or \code{0} for no printed updates.
#'
#' @name verbose_Param
NULL

#' vertices
#'
#' @param vertices A \eqn{V \times 3} matrix, where each row contains the Euclidean
#'  coordinates at which a given vertex in the mesh is located. \eqn{V} is the
#'  number of vertices in the mesh
#'
#' @name vertices_Param
NULL
