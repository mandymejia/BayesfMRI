# Build --> Install and Restart
library(testthat)
library(INLA)
# PARDISO license
inla.setOption(pardiso.license = '../INLA/pardiso.lic')
inla.pardiso.check()
library(ciftiTools)
# Connectome Workbench
ciftiTools.setOption('wb_path', '../workbench/')
source("tests/testthat/test-cifti.R")
