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

# data.frame where each row describes a combination of arguments to test.
sess_avg <- c("1_TRUE", "2_TRUE", "2_FALSE")
prewhiten <- c(TRUE, FALSE)
params <- expand.grid(sess_avg=sess_avg, prewhiten=prewhiten, stringsAsFactors=FALSE)
params$sess <- as.numeric(gsub("_.*", "", sess_avg))
params$avg <- as.logical(gsub(".*_", "", sess_avg))
params$bs <- "left"
params$smooth <- 5
params$Bayes <- FALSE
# Add non-combinatorial test: both hemispheres; no smoothing
params <- rbind(
  params,
  data.frame(
    sess_avg="1_TRUE",
    prewhiten=TRUE,
    sess=1,
    avg=TRUE,
    bs=c("both", "right"),
    smooth=c(5, 0),
    Bayes=FALSE
  )
)
if (doINLA) {
  # Add Bayesian modeling tests: just three
  params <- rbind(
    params,
    data.frame(
      sess_avg=c("1_TRUE", "2_TRUE", "2_FALSE"),
      prewhiten=c(FALSE, TRUE, FALSE),
      sess=c(1, 2, 2),
      avg=c(TRUE, TRUE, FALSE),
      bs=c("right", "both", "right"),
      smooth=c(5, 5, 0),
      Bayes=TRUE
    )
  )
}

# Test each combination.
for (ii in seq(nrow(params))) {
  # Print a description of the combination to test.
  cat(
    params$sess[ii],
    ifelse(
      params$sess[ii]>1,
      ifelse(params$avg[ii], "sessions (w/ averaging),", "sessions (not avg),"),
      "session,"
    ), "and",
    ifelse(params$prewhiten[ii], "with prewhitening", "without prewhitening"),
    ifelse(params$Bayes[ii], ", Bayesian model", ", classical model" ), "\n\n"
  )

  # Do BayesGLM_cifti
  exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
    cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(params$sess[ii])],
    surfL_fname=ciftiTools.files()$surf["left"],
    surfR_fname=ciftiTools.files()$surf["right"],
    brainstructures = params$bs[ii],
    onsets = switch(params$sess[ii], events[[1]], events[seq(2)]),
    TR = 2.2,
    Bayes = params$Bayes[ii],
    ar_order = ifelse(params$prewhiten[ii], 6, 0),
    ar_smooth = params$smooth[ii],
    resamp_res = ifelse(params$Bayes[ii], resamp_res/2, resamp_res) / ifelse(params$bs[ii]=="both", 2, 1),
    verbose = FALSE,
    return_INLA = TRUE,
    combine_sessions = params$avg[ii]
  ))
  print(exec_time)

  if (FALSE) {
    bgroup_fake <- list(bfmri_ii, bfmri_ii)
    z <- BayesGLM2(bgroup_fake)
  }

  # Plot GLM results.
  if (saveResults) {
    plot(
      bfmri_ii, idx=5,
      title=paste0("Tongue, params ", ii),
      fname=file.path(dir_resultThis, paste0("bglm_", ii))
    )
  }

  # [Manual visual comparison if `resamp_res` was high enough]
  # stop()
  # contii <- read_xifti(file.path(dir_data, "derivatives surface_pipeline sub-MSC01 task_contrasts_cifti motor sub-MSC01-motor_contrasts_32k_fsLR.dscalar.nii"))
  # contii$meta$cifti$names
  # plot(...)

  # Identify the activations.
  act_ii <- id_activations_cifti(bfmri_ii, threshold=.01, method=ifelse(params$Bayes[ii], "Bayesian", "classical"), alpha=0.05)
  # plot(act_ii$activations_xifti)
  if (ii == 1) {
    # Test the other arguments too.
    act_ii_temp <- id_activations_cifti(
      bfmri_ii, threshold=.1, method='classical', alpha=0.1, correction='FWER'
      #, excur_method='QC'
    )
    rm(act_ii_temp)
  }
  if (saveResults) {
    plot(
      act_ii$activations_xifti, idx=5,
      title=paste0("Tongue activations, params ", ii),
      fname=file.path(dir_resultThis, paste0("act_", ii))
    )
  }

  if (saveResults) {
    saveRDS(
      list(bfmri=bfmri_ii, act=act_ii, exec_time=exec_time),
      file.path(dir_resultThis, paste0("params", ii, "_MSC.rds"))
    )
  }
}

file.remove(file.path(dir_results, "bfmri_out_left.rds"))
file.remove(file.path(dir_results, "bfmri_out_right.rds"))

# BayesGLM2?

# Export summary graphics ------------------------------------------------------
params$sess_avg <- NULL
library(gridExtra)
png(
  file.path(dir_resultThis, "params.png"),
  height = 50*nrow(params), width = 200*ncol(params)
)
grid.table(params)
dev.off()

ctast <- read_xifti(file.path(dir_data, "derivatives surface_pipeline sub-MSC01 task_contrasts_cifti motor sub-MSC01-motor_contrasts_32k_fsLR.dscalar.nii"))
plot(
  ctast, idx=5, zlim=c(-5, 5),
  title=paste0("MSC tongue contrast"),
  fname=file.path(dir_resultThis, paste0("MSC_tongue_contrast"))
)

ciftiTools:::view_comp(
  file.path(dir_resultThis, c(paste0("bglm_", seq(1, 11), ".png"), "MSC_tongue_contrast.png")),
  fname=file.path(dir_resultThis, "bglm_comp.png")
)

ciftiTools:::view_comp(
  file.path(dir_resultThis, c(paste0("act_", seq(1, 11), ".png"), "MSC_tongue_contrast.png")),
  fname=file.path(dir_resultThis, "act_comp.png")
)
