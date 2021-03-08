# Build --> Install and Restart
library(testthat)
library(INLA)
inla.setOption(pardiso.license = '~/pardiso.lic')
inla.pardiso.check()
library(ciftiTools)
ciftiTools.setOption('wb_path', '~/workbench/')
source("tests/testthat/test-cifti.R")
