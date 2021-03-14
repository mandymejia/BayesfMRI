# Build --> Install and Restart
library(testthat)
library(INLA)
# [Edit this] path to your PARDISO license
inla.setOption(pardiso.license = '../INLA/pardiso.lic')
inla.pardiso.check()
library(ciftiTools)
# [Edit this] path to your Connectome Workbench
ciftiTools.setOption('wb_path', '../workbench/')
# [Edit this] path to test data
dir_data <- "../BayesfMRI_testData"
# [Edit this] path to results from tests
dir_results <- file.path(dir_data, "results")
source("tests/testthat/test-classical.R")
source("tests/testthat/test-Bayes.R")