# [Build --> Install and Restart]

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
n_sess <- 1
BayesGLM_cifti_args <- list(
  cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(n_sess)],
  surfL_fname=ciftiTools.files()$surf["left"],
  surfR_fname=ciftiTools.files()$surf["right"],
  brainstructures = "both",
  onsets = switch(n_sess, events[seq(3)], list(events[seq(3)], events[seq(4,6)])),
  TR = 2.2,
  dHRF=2,
  nuisance=switch(n_sess, nuis$rp_1, nuis),
  Bayes = TRUE,
  ar_order = 0,
  ar_smooth = 3,
  resamp_res = 800,
  verbose = FALSE,
  return_INLA = "trim",
  combine_sessions = FALSE
)
bglm <- do.call(BayesGLM_cifti, BayesGLM_cifti_args)
bglm2 <- do.call(BayesGLM_cifti, c(list(EM=TRUE), BayesGLM_cifti_args))

# stop here

# Change `brainstructures`, `n_sess`, and `meanTol`.
# Expect same classical model results, where they exist.
n_sess <- 1
BayesGLM_cifti_args <- list(
  cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(n_sess)],
  surfL_fname=ciftiTools.files()$surf["left"],
  surfR_fname=ciftiTools.files()$surf["right"],
  brainstructures = "left",
  onsets = switch(n_sess, events[rev(seq(3))], list(events[rev(seq(3))], events[rev(seq(4,6))])),
  TR = 2.2,
  dHRF=2,
  nuisance=switch(n_sess, nuis$rp_1, nuis),
  Bayes = TRUE,
  ar_order = 0,
  ar_smooth = 3,
  resamp_res = 800,
  verbose = TRUE,
  return_INLA = "full",
  meanTol=9999,
  combine_sessions = FALSE
)
bglm2 <- do.call(BayesGLM_cifti, BayesGLM_cifti_args)

z <- bglm$BayesGLM_results$cortex_left$result_classical
y <- bglm2$BayesGLM_results$cortex_left$result_classical
q <- y$single_session$estimates - z$session1$estimates[,c(3,2,1)]
