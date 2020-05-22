# BayesfMRI
BayesfMRI R package 

Two main functions:
* BayesGLM – Implements the spatial Bayesian GLM for task fMRI on the cortical surface proposed by Mejia et al. 2019a [https://doi.org/10.1080/01621459.2019.1611582]
* templateICA – Implements the template ICA model for estimating functional brain networks in individual subjects proposed by Mejia et al. 2019b [https://doi.org/10.1080/01621459.2019.1611582]


## Important Note on Dependencies:

The INLA package is required, which, due to a CRAN policy, will not be installed automatically. You can obtain the two packages by running `install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`.  

The default R-INLA binaries are built on Ubuntu1604. Instructions on how to obtain binaries for other Linux builds are available at http://www.r-inla.org/events/alternativelinuxbuilds.

  
