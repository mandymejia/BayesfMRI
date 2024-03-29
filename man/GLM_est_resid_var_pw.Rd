% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLM_est_resid_var_pw.R
\name{GLM_est_resid_var_pw}
\alias{GLM_est_resid_var_pw}
\title{Standardize data variance, and prewhiten if applicable}
\usage{
GLM_est_resid_var_pw(
  BOLD,
  design,
  spatial,
  spatial_type,
  session_names,
  field_names,
  design_type,
  valid_cols,
  nT,
  nD,
  ar_order,
  ar_smooth,
  aic,
  n_threads,
  do_pw
)
}
\arguments{
\item{BOLD, design, spatial, spatial_type}{See \code{BayesGLM}.}

\item{session_names, field_names, design_type}{See \code{BayesGLM}.}

\item{valid_cols, nT, nD, do_pw}{See \code{BayesGLM}.}
}
\value{
List of results
}
\description{
Standardize data variance and prewhiten if applicable, for the GLM.
}
\keyword{internal}
