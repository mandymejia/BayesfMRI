% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_bayesglm_utils.R
\name{extract_estimates}
\alias{extract_estimates}
\title{Extract Estimates of Activation}
\usage{
extract_estimates(INLA_model_obj, session_names, spatial, spde, stat = "mean")
}
\arguments{
\item{INLA_model_obj}{An object of class \code{"inla"}, a result of a call to
\code{inla}.}

\item{session_names}{Vector of fMRI session names}

\item{stat}{A string representing the posterior summary statistic to be returned}
}
\value{
Estimates from inla model
}
\description{
Obtains the posterior mean or other summary statistic for each latent field
}
\keyword{internal}
