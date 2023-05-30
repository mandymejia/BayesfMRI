# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 9000
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_results <- "tests/results_notInPackage"
thisResultName <- gsub(
  ".", "_",
  as.character(packageVersion("BayesfMRI")[[1]]), fixed=TRUE
)
dir_resultsThis <- file.path(dir_results, thisResultName, "sim")
if (!overwriteResults && dir.exists(dir_resultsThis)) { stop("Results exist already.") }
if (!dir.exists(dir_resultsThis)) { dir.create(dir_resultsThis) }

library(testthat)
if (doINLA) { library(INLA) }
library(fMRItools)
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
library(BayesfMRI)

# Simulate data ----------------------------------------------------------------
bsim <- readRDS(file.path(dir_results, "simData/bsim.rds"))

# BayesGLM ---------------------------------------------------------------------
bfmri_1sess <- lapply(bsim$simulated_cifti[seq(3)], BayesGLM_cifti,
  brainstructures = "both",
  design = bsim$design,
  Bayes = TRUE,
  ar_order = 6,
  ar_smooth = 5,
  resamp_res = 3000,
  verbose = FALSE,
  return_INLA = "trimmed"
)
saveRDS(bfmri_1sess, file.path(dir_resultsThis, "bfmri_1sess.rds"))

bfmri_2sess <- lapply(
  list(
    bsim$simulated_cifti[c(1,4)],
    bsim$simulated_cifti[c(2,5)],
    bsim$simulated_cifti[c(3,6)]
  ),
  BayesGLM_cifti,
  brainstructures = "left",
  design = list(bsim$design, bsim$design),
  Bayes = TRUE,
  ar_order = 6,
  ar_smooth = 5,
  resamp_res = 3000,
  verbose = FALSE,
  combine_sessions = FALSE,
  return_INLA = "trimmed"
)
saveRDS(bfmri_2sess, file.path(dir_resultsThis, "bfmri_2sess.rds"))

# BayesGLM2 --------------------------------------------------------------------
bfmri_1sess <- readRDS(file.path(dir_resultsThis, "bfmri_1sess.rds"))
b2 <- BayesGLM2(bfmri_1sess, excursion_type = ">", quantiles=.7)
saveRDS(b2, "~/Downloads/temp.rds")
