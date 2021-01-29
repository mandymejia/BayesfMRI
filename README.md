# BayesfMRI

`BayesfMRI` R package. Includes the main function `BayesGLM`, which implements the spatial Bayesian GLM for task fMRI on the cortical surface proposed by Mejia et al. 2019a [https://doi.org/10.1080/01621459.2019.1611582].  Also contains two wrapper functions:

* `BayesGLM_cifti` – implements `BayesGLM` on CIFTI-format (greyordinates) fMRI data
* `BayesGLM_slice` - implements `BayesGLM` on slice-wise volumetric fMRI data

## Installation

The `INLA` package must be installed before `BayesfMRI`. Due to a CRAN policy it not be installed automatically. You can obtain `INLA` by running `install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=FALSE)`. The default R-INLA binaries are built on Ubuntu1604. Instructions on how to obtain binaries for other Linux builds are available at http://www.r-inla.org/events/alternativelinuxbuilds.

After installing `INLA`, you can install the development version of `BayesfMRI` from Github with:

``` r
devtools::install_github("mandymejia/BayesfMRI")
```

Additionally, the CIFTI-related functions use the `ciftiTools` package, which requires the Connectome Workbench. It can be installed from the [HCP website](https://www.humanconnectome.org/software/get-connectome-workbench). You will need to point `ciftiTools` to the Workbench location:

``` r
library(ciftiTools)
ciftiTools.setOption("wb_path", "path/to/workbench")
library(BayesfMRI)
```