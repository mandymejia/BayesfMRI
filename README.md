# BayesfMRI

BayesfMRI R package. Includes the main function `BayesGLM`, which implements the spatial Bayesian GLM for task fMRI on the cortical surface proposed by Mejia et al. 2019a [https://doi.org/10.1080/01621459.2019.1611582].  Also contains two wrapper functions:

* `BayesGLM_cifti` – implements `BayesGLM` on CIFTI-format (greyordinates) fMRI data
* `BayesGLM_slice` - implements `BayesGLM` on slice-wise volumetric fMRI data


## Important Note on Dependencies:

The INLA package is required, which, due to a CRAN policy, will not be installed automatically. You can obtain the two packages by running `install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=FALSE)`.  The default R-INLA binaries are built on Ubuntu1604. Instructions on how to obtain binaries for other Linux builds are available at http://www.r-inla.org/events/alternativelinuxbuilds.

  
