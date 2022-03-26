# Build --> Install and Restart
# [Edit these] paths to Workbench and result output.
my_wb <- "../workbench"
save_path <- "tests/output" # in gitignore

library(testthat)
library(BayesfMRI)
library(ciftiTools)
if (interactive()) { ciftiTools.setOption("wb_path", my_wb) }
library(brainSim)

tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }
# source(file.path(tests_dir, "test-Bayes.R"))
# source(file.path(tests_dir, "test-classical.R"))
# source(file.path(tests_dir, "test-misc.R"))
source(file.path(tests_dir, "test-single.R")) # won't require interactive
source(file.path(tests_dir, "test-new.R")) # will require interactive
