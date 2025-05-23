% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesGLM2_Bayes_utils.R
\name{qsample}
\alias{qsample}
\title{Sample from a multivariate normal with mean and precision}
\usage{
qsample(n, mu, Q)
}
\arguments{
\item{n}{number of samples}

\item{mu}{mean vector (length = p)}

\item{Q}{sparse p x p positive definite precision matrix (class = dgCMatrix)}
}
\value{
An n x p matrix of samples
}
\description{
Sample from a multivariate normal with mean and precision
}
\keyword{internal}
