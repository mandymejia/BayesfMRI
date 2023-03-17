# Build --> Install and Restart

# [Edit these]
# ## path to your PARDISO license
# my_pardiso <- "../INLA/pardiso.lic"
## path to your Connectome Workbench
my_wb <- "~/Desktop/workbench"
## path to test data
dir_data <- "tests/data"
## path to results from tests
dir_results <- file.path(dir_data, "results")

library(testthat)
# library(INLA)
# # [Edit this] path to your PARDISO license
# if (interactive()) {
#   inla.setOption(pardiso.license = my_pardiso)
#   inla.pardiso.check()
# }
library(brainSim)
library(ciftiTools)
if (interactive()) { ciftiTools.setOption('wb_path', my_wb) }
library(BayesfMRI)

tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }
source(file.path(tests_dir, "test-classical.R"))
source(file.path(tests_dir, "test-misc.R"))
