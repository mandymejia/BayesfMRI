# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 6000
my_pardiso <- "~/Documents/pardiso.lic" # INLA PARDISO license
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_data <- "tests/data_notInPackage"
dir_results <- "tests/results_notInPackage"
thisResultName <- gsub(
  ".", "_",
  as.character(packageVersion("BayesfMRI")[[1]]), fixed=TRUE
)
dir_resultThis <- file.path(dir_results, thisResultName)
if (!overwriteResults && dir.exists(dir_resultThis)) { stop("Results exist already.") }
if (!dir.exists(dir_resultThis)) { dir.create(dir_resultThis) }

library(testthat)
library(stats) # setNames
if (doINLA) {
  library(INLA)
  inla.setOption(pardiso.license = my_pardiso)
  inla.pardiso.check()
}
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
library(BayesfMRI)

# Get file names.
fnames <- list(
  #tmask = "derivatives surface_pipeline sub-MSC01 processed_task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_bold_32k_fsLR_tmask.txt"
  cifti_1 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-01_bold.dtseries.nii",
  events_1 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-01_events.tsv",
  cifti_2 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-02_bold.dtseries.nii",
  events_2 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-02_events.tsv"
)
fnames <- lapply(fnames, function(x){file.path(dir_data, x)})

events <- lapply(fnames[c("events_1", "events_2")], read.table, header=TRUE)
task_names <- c("RHand", "LHand", "RFoot", "LFoot", "Tongue")
events <- lapply(events, function(x){
  setNames(lapply(task_names, function(y){
    cbind(
      onset = x$onset[x$trial_type==y],
      duration = x$duration[x$trial_type==y]
    )
  }), task_names)
})

# BayesGLM ---------------------------------------------------------------------

n_sess <- 2
# Do BayesGLM_cifti
exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
  cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(n_sess)],
  surfL_fname=ciftiTools.files()$surf["left"],
  surfR_fname=ciftiTools.files()$surf["right"],
  brainstructures = "left",
  onsets = switch(n_sess, events[[1]], events[seq(2)]),
  TR = 2.2,
  Bayes = TRUE,
  hpf = 0,
  ar_order = 6,
  ar_smooth = 5,
  resamp_res = resamp_res,
  verbose = FALSE,
  return_INLA = TRUE,
  combine_sessions = TRUE
))
print(exec_time)

act_ii <- id_activations(bfmri_ii, gamma=.2, alpha=0.01)
plot(act_ii, idx=seq(5), together="idx", together_ncol=1, fname=file.path(dir_resultThis, "MSC/bfmri_act.png"))

## [TO DO] and if I add more noise, what happens?

ctast <- read_xifti(
  file.path(dir_data, "derivatives surface_pipeline sub-MSC01 task_contrasts_cifti motor sub-MSC01-motor_contrasts_32k_fsLR.dscalar.nii")
)
plot(ctast, idx=seq(5), zlim=c(-5, 5), together="idx", together_ncol=1, fname=file.path(dir_resultThis, "MSC/motorContrasts.png"))
