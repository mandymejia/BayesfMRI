#' INLA
#'
#' @section INLA Requirement:
#'  This function requires the \code{INLA} package, which is not a CRAN package.
#'  See \url{https://www.r-inla.org/download-install} for easy installation instructions.
#'
#' @name INLA_Description
NULL

#' avg_sessions
#'
#' @param avg_sessions Average estimates for betas over multiple
#'  sessions? Default: \code{FALSE}.
#'
#' @name avg_sessions_Param
NULL

#' contrasts: INLA
#'
#' @param contrasts List of contrast vectors to be passed to \code{inla::inla}.
#'
#' @name contrasts_Param_inla
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

#' return_INLA_result: TRUE
#'
#' @param return_INLA_result Return the INLA model object? (It can be large.)
#'  Default: \code{TRUE}. Required for running \code{id_activations}
#'  after, but not for running BayesGLM_joint after to get
#'  posterior quantities of group means or contrasts.
#'
#' @name return_INLA_result_Param_TRUE
NULL

#' return_INLA_result: FALSE
#'
#' @param return_INLA_result Return the INLA model object? (It can be large.)
#'  Default: \code{FALSE}. Required for running \code{id_activations}
#'  after, but not for running BayesGLM_joint after to get
#'  posterior quantities of group means or contrasts.
#'
#' @name return_INLA_result_Param_FALSE
NULL


#' scale_BOLD
#'
#' @param scale_BOLD Scale timeseries data so estimates represent percent signal
#'  change? Default: \code{TRUE}. If \code{FALSE}, the data will just be
#'  centered.
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

#'  Bayes
#'
#' @param Bayes If TRUE, will fit a spatial Bayesian GLM in addition to the
#' classical GLM. Classical GLM results are always returned.
#'
#' @name Bayes_Param
NULL

