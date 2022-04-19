# BayesfMRI

`BayesfMRI` R package. Includes the main function `BayesGLM`, which implements the spatial Bayesian GLM for task fMRI on the cortical surface proposed by Mejia et al. 2019a [https://doi.org/10.1080/01621459.2019.1611582].  Also contains two wrapper functions:

* `BayesGLM_cifti` – implements `BayesGLM` on CIFTI-format (greyordinates) fMRI data
* `BayesGLM_slice` - implements `BayesGLM` on slice-wise volumetric fMRI data

## Important Note on Dependencies:

`BayesfMRI` depends on the `ciftiTools` package, which requires an installation of Connectome Workbench.  It can be installed from the [HCP website](https://www.humanconnectome.org/software/get-connectome-workbench).

The INLA package is required, which, due to a CRAN policy, will not be installed automatically. You can obtain it by running `install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=FALSE)`.  The default R-INLA binaries are built on Ubuntu1604. Instructions on how to obtain binaries for other Linux builds are available at https://www.r-inla.org/events/alternativelinuxbuilds.  **Note: INLA must be installed before installing `BayesfMRI`.**

An INLA-PARDISO license is also required for computational efficiency when using the INLA implementation of the Bayesian GLM. **As of March 21, 2022, PARDISO is no longer offering free academic licenses. Our team is currently at work to implement the Bayesian GLM using an EM algorithm that does not depend on INLA or PARDISO. Please watch this space for updates concerning this new implementation. For now, a beta version that does not require PARDISO can be found on the 1.8.EM branch.**  To obtain an INLA-PARDISO license, run `inla.pardiso()` in R after running `library(INLA)`. Once you obtain a license, point to it using `INLA::inla.setOption(pardiso.license = "pardiso.lic")` followed by `INLA::inla.pardiso.check()` to ensure that PARDISO is successfully installed and running.
