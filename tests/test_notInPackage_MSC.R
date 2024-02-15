# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 7000
# my_pardiso <- "~/Documents/pardiso.lic" # INLA PARDISO license
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
  # inla.setOption(pardiso.license = my_pardiso)
  # inla.pardiso.check()
}
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
library(BayesfMRI)

# Get file names.
fnames <- list(
  tmask = "derivatives surface_pipeline sub-MSC01 processed_task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_bold_32k_fsLR_tmask.txt",
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
### Classical vs Bayes; Single- vs Multi-session -----
BayesGLM_cifti_args <- function(n_sess, resamp_factor=1){
  list(
    cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(n_sess)],
    surfL_fname=ciftiTools.files()$surf["left"],
    surfR_fname=ciftiTools.files()$surf["right"],
    brainstructures = "both",
    onsets = switch(n_sess, events[[1]], events[seq(2)]),
    TR = 2.2,
    #Bayes = TRUE,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = resamp_res * resamp_factor,
    verbose = TRUE,
    return_INLA = "trim"
  )
}

##### First pass to detect errors
bglm_c1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1, resamp_factor=.1)))
bglm_b1 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(1, resamp_factor=.1)))
bglm_c2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2, resamp_factor=.1)))
bglm_b2 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(2, resamp_factor=.1)))

##### Second pass to get results of decent resolution
bglm_c1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1)))
bglm_b1 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(1)))
bglm_c2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2)))
bglm_b2 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(2)))

act_c1 <- id_activations(bglm_c1, gamma=.01, sessions=1)
act_b1 <- id_activations(bglm_b1, gamma=.01, sessions=1)
act_c2 <- id_activations(bglm_c2, gamma=.01, sessions=seq(2))
act_b2 <- id_activations(bglm_b2, gamma=.01, sessions=seq(2))

### Misc. cases; not checking these results, but checking for errors
### Last updated: 5.1
bglm_m1 <- BayesGLM_cifti(
  cifti_fname = fnames$cifti_1,
  brainstructures = "right",
  onsets = events[[1]],
  TR = 2.2,
  dHRF_as="field",
  Bayes = TRUE,
  resamp_res=resamp_res,
  ar_order = 0,
  verbose = 0
)
act_m1 <- id_activations(bglm_m1, alpha=.1, gamma=.05, fields=1)

bglm_m2 <- BayesGLM_cifti(
  cifti_fname = c(fnames$cifti_1, fnames$cifti_2),
  brainstructures = "left",
  onsets = events[seq(2)],
  TR = 2.2,
  dHRF=0,
  Bayes = FALSE,
  resamp_res=resamp_res*2,
  nuisance=nuis,
  ar_order = 2,
  ar_smooth = 0,
  verbose = 2
)
act_m2 <- id_activations(bglm_m2)

# Save ---
save(list=ls(), file=file.path(dir_resultThis, "test_notInPackage_MSC.rda"))

# Plot ---
plot(bglm_c1, fname=file.path(dir_resultThis, "MSC_bglm_c1"), together="idx")
plot(bglm_b1, fname=file.path(dir_resultThis, "MSC_bglm_b1"), together="idx")
plot(bglm_c2, fname=file.path(dir_resultThis, "MSC_bglm_c2"), together="idx")
plot(bglm_b2, fname=file.path(dir_resultThis, "MSC_bglm_b2"), together="idx")

plot(act_c1, fname=file.path(dir_resultThis, "MSC_act_c1"), together="idx")
plot(act_b1, fname=file.path(dir_resultThis, "MSC_act_b1"), together="idx")
plot(act_c2, fname=file.path(dir_resultThis, "MSC_act_c2"), together="idx")
plot(act_b2, fname=file.path(dir_resultThis, "MSC_act_b2"), together="idx")

plot(bglm_m1, fname=file.path(dir_resultThis, "MSC_bglm_m1"), together="idx")
plot(bglm_b2, fname=file.path(dir_resultThis, "MSC_bglm_m2"), together="idx")

plot(act_m1, fname=file.path(dir_resultThis, "MSC_act_m1"), together="idx")
plot(act_m2, fname=file.path(dir_resultThis, "MSC_act_m2"), together="idx")
