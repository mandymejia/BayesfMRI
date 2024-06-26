% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prevalence.R
\name{prevalence}
\alias{prevalence}
\title{Activations prevalence.}
\usage{
prevalence(act_list, gamma_idx = 1)
}
\arguments{
\item{act_list}{List of activations from \code{\link{activations}}. All
should have the same sessions, fields, and brainstructures.}

\item{gamma_idx}{If activations at multiple thresholds were computed, which
threshold should be used for prevalence? Default: the first (lowest).}
}
\value{
A list containing the prevalences of activation, as a proportion of
the results from \code{act_list}.
}
\description{
Activations prevalence.
}
