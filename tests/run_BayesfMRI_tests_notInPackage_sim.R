# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 6000
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_results <- "tests/results_notInPackage"
thisResultName <- gsub(
  ".", "_",
  as.character(packageVersion("BayesfMRI")[[1]]), fixed=TRUE
)
dir_resultThis <- file.path(dir_results, thisResultName, "sim")
if (!overwriteResults && dir.exists(dir_resultThis)) { stop("Results exist already.") }
if (!dir.exists(dir_resultThis)) { dir.create(dir_resultThis) }

library(testthat)
if (doINLA) { library(INLA) }
library(fMRItools)
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
library(BayesfMRI)

# Simulate data ----------------------------------------------------------------
cat("Simulating data\n~~~~~~~~~~~~~~~~\n")
source("tests/simulate_cifti.R")
set.seed(0)
bsim <- simulate_cifti_multiple(
  wb_path=my_wb,
  brainstructures="both",
  n_subjects=1,
  n_sessions=2,
  n_runs=1,
  ntasks=3,
  ntime=400,
  session_var=3,
  run_var=2,
  resamp_res=resamp_res
)

# BayesGLM ---------------------------------------------------------------------
cat("BayesGLM\n~~~~~~~~~~~~~~~~\n")

# combinatorial parameters
Bayes <- c(FALSE, TRUE)
bs <- c("right", "both")
# parameters w/ rising power
params_pwr <- c("sess1-avgNo-pwYes", "sess2-avgNo-pwYes", "sess2-avgYes-pwYes", "sess2-avgYes-pwNo")
params <- expand.grid(params_pwr=params_pwr, bs=bs, Bayes=Bayes, stringsAsFactors=FALSE)
# add  no smoothing + no PW + avg for classical
params$smooth <- 5
params <- rbind(
  data.frame(bs="right", Bayes=FALSE, params_pwr="sess2-avgYes-pwNo", smooth=0),
  params
)
# convert `params_pwr` column to individual params
params$sess <- as.numeric(gsub("sess", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 1)))
params$avg <- "Yes" == gsub("avg", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 2))
params$pw <- "Yes" == gsub("pw", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 3))
params$params_pwr <- NULL

# Test each combination.
for (ii in seq(nrow(params))) {
  print(params[ii,])

  xii_ii <- bsim$simulated_cifti[seq(params$sess[ii])]
  # if (params$bs[ii] != "both") {
  #   xii_ii <- lapply(xii_ii, remove_xifti, "cortex_left")
  # }

  # Do BayesGLM_cifti
  exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
    cifti_fname = bsim$simulated_cifti[seq(params$sess[ii])],
    surfL_fname=ciftiTools.files()$surf["left"],
    surfR_fname=ciftiTools.files()$surf["right"],
    brainstructures = params$bs[ii],
    design = switch(params$sess[ii], bsim$design, list(bsim$design, bsim$design)),
    Bayes = params$Bayes[ii],
    ar_order = ifelse(params$pw[ii], 6, 0),
    ar_smooth = params$smooth[ii],
    resamp_res = ifelse(params$Bayes[ii], resamp_res/2, resamp_res) / ifelse(params$bs[ii]=="both", 2, 1),
    verbose = FALSE,
    return_INLA = "trimmed",
    avg_sessions = params$avg[ii]
  ))
  print(exec_time)

  if (saveResults) {
    saveRDS(
      bfmri_ii,
      file.path(dir_resultThis, paste0("bfmri", ii, "_HCP.rds"))
    )
  }
}

# # Activations and plot ---------------------------------------------------------
# for (ii in seq(nrow(params))) {
#   if (FALSE) {
#     bgroup_fake <- list(bfmri_ii, bfmri_ii)
#     z <- BayesGLM2(bgroup_fake)
#   }
#
#   # Plot GLM results.
#   if (saveResults) {
#     plot(
#       bfmri_ii, idx=1,
#       title=paste0("Win, params ", ii),
#       fname=file.path(dir_resultThis, paste0("bglm_", ii))
#     )
#   }
#
#   # [Manual visual comparison if `resamp_res` was high enough]
#   # stop()
#   # contii <- read_xifti(file.path(dir_data, "derivatives surface_pipeline sub-MSC01 task_contrasts_cifti motor sub-MSC01-motor_contrasts_32k_fsLR.dscalar.nii"))
#   # contii$meta$cifti$names
#   # plot(...)
#
#   # Identify the activations.
#   act_ii <- id_activations_cifti(bfmri_ii, threshold=.01, method=ifelse(params$Bayes[ii], "Bayesian", "classical"), alpha=0.05)
#   # plot(act_ii$activations_xifti)
#   if (ii == 1) {
#     # Test the other arguments too.
#     act_ii_temp <- id_activations_cifti(
#       bfmri_ii, threshold=.1, method='classical', alpha=0.1, correction='FWER'
#       #, excur_method='QC'
#     )
#     rm(act_ii_temp)
#   }
#   if (saveResults) {
#     plot(
#       act_ii$activations_xifti, idx=1,
#       title=paste0("Win activations, params ", ii),
#       fname=file.path(dir_resultThis, paste0("act_", ii))
#     )
#   }
#
#   if (saveResults) {
#     saveRDS(
#       list(bfmri=bfmri_ii, act=act_ii, exec_time=exec_time),
#       file.path(dir_resultThis, paste0("params", ii, "_HCP.rds"))
#     )
#   }
# }
#
# # BayesGLM2?
