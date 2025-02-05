
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesfMRI <img src="man/figures/logo.png" align="right" height="139" alt="BayesfMRI sticker" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/BayesfMRI)](https://cran.r-project.org/package=BayesfMRI)
[![R-CMD-check](https://github.com/mandymejia/BayesfMRI/workflows/R-CMD-check/badge.svg)](https://github.com/mandymejia/BayesfMRI/actions)
[![Codecov test
coverage](https://codecov.io/gh/mandymejia/BayesfMRI/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mandymejia/BayesfMRI?branch=master)
<!-- badges: end -->

The `BayesfMRI` R package includes the main function `BayesGLM`, which
implements a spatial Bayesian GLM for task fMRI. It also contains a
wrapper function `BayesGLM_cifti`, for CIFTI cortical surface fMRI data.

<!-- * `BayesGLM_vol` - implements `BayesGLM` on NIFTI subcortical voxel fMRI data -->

## Citation

If you use `BayesfMRI` please cite the following papers:

| Name                                                                                   | APA Citation                                                                                                                                                                                                                       |
|----------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Spatial Bayesian GLM](https://doi.org/10.1080/01621459.2019.1611582)                  | Mejia, A. F., Yue, Y., Bolin, D., Lindgren, F., & Lindquist, M. A. (2020). A Bayesian general linear modeling approach to cortical surface fMRI data analysis. Journal of the American Statistical Association, 115(530), 501-520. |
| [Multi-session Spatial Bayesian GLM](https://doi.org/10.1016/j.neuroimage.2022.118908) | Spencer, D., Yue, Y. R., Bolin, D., Ryan, S., & Mejia, A. F. (2022). Spatial Bayesian GLM on the cortical surface produces reliable task activations in individuals and groups. NeuroImage, 249, 118908.                           |

You can also obtain citation information from within R like so:

``` r
citation("BayesfMRI")
```

## Installation

You can install `BayesfMRI` from [CRAN](https://cran.r-project.org/)
with:

``` r
install.packages("BayesfMRI")
```

## Vignette

See [this
link](https://htmlpreview.github.io/?https://github.com/mandymejia/BayesfMRI/blob/10.0/vignettes/BayesfMRI_vignette.html)
to view the tutorial vignette.

## Important Note on Dependencies:

`BayesfMRI` depends on the `ciftiTools` package, which requires an
installation of Connectome Workbench. It can be installed from the [HCP
website](https://www.humanconnectome.org/software/get-connectome-workbench).

<!-- By default, the spatial Bayesian model in `BayesGLM` is implemented using an expectation-maximization algorithm written in C++. To instead use the INLA package, set `EM=FALSE`. The INLA package will be required, as well as an INLA-PARDISO license for computational efficiency. -->

The INLA package is required, which, due to a CRAN policy, will not be
installed automatically. You can obtain it by running
`install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=FALSE)`.
For more information, see [the INLA
website](https://www.r-inla.org/download-install). **Note: INLA must be
installed before installing `BayesfMRI`.**

On Mac platforms, an installation of
[Xcode](https://mac.r-project.org/tools/) is necessary to build the C++
code included in `BayesfMRI`.

<!--An INLA-PARDISO license is also required for computational efficiency.  To obtain an INLA-PARDISO license, run `inla.pardiso()` in R after running `library(INLA)`. Once you obtain a license, point to it using `INLA::inla.setOption(pardiso.license = "pardiso.lic")` followed by `INLA::inla.pardiso.check()` to ensure that PARDISO is successfully installed and running. 
-->
