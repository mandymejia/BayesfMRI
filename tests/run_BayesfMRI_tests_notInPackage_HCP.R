# [Build --> Install and Restart]

# TO DO: diff plot names for HCP vs MSC

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 4000
my_pardiso <- "~/Documents/pardiso.lic" # INLA PARDISO license
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_data <- "/Users/ddpham/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/O365-BL-STAT-StatMIND-Projects - General/Data/bfMRI"
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
  cifti_1 = "151526/tfMRI_GAMBLING_LR_Atlas.dtseries.nii",
  events_1w = "151526/win_event.txt",
  events_1l = "151526/loss_event.txt",
  events_1n = "151526/neut_event.txt",
  rp_1 = "151526/Movement_Regressors.txt",
  cifti_2 = "151526_retest/tfMRI_GAMBLING_LR_Atlas.dtseries.nii",
  events_2w = "151526_retest/win_event.txt",
  events_2l = "151526_retest/loss_event.txt",
  events_2n = "151526_retest/neut_event.txt",
  rp_2 = "151526/Movement_Regressors.txt"
)
fnames <- lapply(fnames, function(x){file.path(dir_data, x)})

events <- lapply(fnames[grepl("events", names(fnames))], read.table, header=FALSE)
events <- lapply(events, function(x){
  x <- x[,seq(2)]
  colnames(x) <- c("onset", "duration")
  x
})
names(events) <- rep(c("win", "loss", "neut"), 2)

nuis <- lapply(fnames[grepl("rp", names(fnames))], read.table, header=FALSE)
nuis <- lapply(nuis, as.matrix)

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
    bs=c("both", "left"),
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
      bs=c("left", "both", "left"),
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

  if (!params$Bayes[ii]) { next }

  # Do BayesGLM_cifti
  exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
    cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(params$sess[ii])],
    surfL_fname=ciftiTools.files()$surf["left"],
    surfR_fname=ciftiTools.files()$surf["right"],
    brainstructures = params$bs[ii],
    onsets = switch(params$sess[ii], events[seq(3)], list(events[seq(3)], events[seq(4,6)])),
    TR = 2.2,
    dHRF=2,
    nuisance=switch(params$sess[ii], nuis$rp_1, nuis),
    Bayes = params$Bayes[ii],
    ar_order = ifelse(params$prewhiten[ii], 6, 0),
    ar_smooth = params$smooth[ii],
    resamp_res = ifelse(params$Bayes[ii], resamp_res/2, resamp_res) / ifelse(params$bs[ii]=="both", 2, 1),
    verbose = FALSE,
    return_INLA = TRUE,
    outfile = file.path(dir_results, "bfmri_out"),
    avg_sessions = params$avg[ii]
  ))
  print(exec_time)
  stop()

  if (FALSE) {
    bgroup_fake <- list(bfmri_ii, bfmri_ii)
    z <- BayesGLM2(bgroup_fake)
  }

  # Plot GLM results.
  if (saveResults) {
    plot(
      bfmri_ii, idx=1,
      title=paste0("Win, params ", ii),
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
      act_ii$activations_xifti, idx=1,
      title=paste0("Win activations, params ", ii),
      fname=file.path(dir_resultThis, paste0("act_", ii))
    )
  }

  if (saveResults) {
    saveRDS(
      list(bfmri=bfmri_ii, act=act_ii, exec_time=exec_time),
      file.path(dir_resultThis, paste0("params", ii, "_HCP.rds"))
    )
  }
}

file.remove(file.path(dir_results, "bfmri_out_left.rds"))
file.remove(file.path(dir_results, "bfmri_out_right.rds"))

# BayesGLM2?

