% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesGLM2_Bayes_utils.R
\name{cholQsample}
\alias{cholQsample}
\title{Sample from the multivariate normal distribution with Cholesky(Q)}
\usage{
cholQsample(n, mu, cholQ)
}
\arguments{
\item{n}{number of samples}

\item{mu}{mean vector}

\item{cholQ}{Cholesky decomposition of the precision (found via \code{Matrix::Cholesky(Q)})}
}
\value{
An \eqn{n \times p} matrix of samples from the MVN distribution,
where \eqn{p} is the length of \code{mu}.
}
\description{
Sample from the multivariate normal distribution with Cholesky(Q)
}
\keyword{internal}
