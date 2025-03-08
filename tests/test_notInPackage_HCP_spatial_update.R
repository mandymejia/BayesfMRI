# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 2200
# my_pardiso <- "~/Documents/pardiso.lic" # INLA PARDISO license
my_wb <- "~/Applications/workbench" # path to your Connectome Workbench

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
library(stats)
if (doINLA) {
  library(INLA)
  # inla.setOption(pardiso.license = my_pardiso)
  # inla.pardiso.check()
}
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
roxygen2::roxygenize("~/Documents/GitHub/BayesfMRI")
#library(BayesfMRI)

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

nTime <- 253#vapply(lapply(fnames[c("cifti_1", "cifti_2")], read_cifti), ncol, 0)

# From `onsets`.
des <- lapply(
  list(events[seq(3)], events[seq(4,6)]),
  function(q) { make_design(EVs=q, TR=.72, nTime=nTime, dHRF=1) }
)
for (ss in seq(2)) {
  nuis[[ss]] <- cbind(nuis[[ss]], des[[ss]]$design[,seq(4,6)])
  des[[ss]] <- des[[ss]]$design[,seq(3)]
}


# BayesGLM ---------------------------------------------------------------------
BOLD <- lapply(fnames[c("cifti_1", "cifti_2")], read_cifti)
BOLD[[2]] <- BOLD[[1]]
BOLD[[1]]$data$cortex_left[c(1,11,111),] <- NA
BOLD[[2]]$data$cortex_left[c(30),] <- NA
BOLD[[1]]$data$subcort[which(BOLD[[1]]$meta$subcort$labels=="Amygdala-L")[c(4,14)],] <- 0
BOLD[[2]]$data$subcort[which(BOLD[[2]]$meta$subcort$labels=="Thalamus-L")[c(4,14)],] <- 0

bglm1 <- BayesGLM(
  BOLD = BOLD[[1]],
  design = des[[1]],
  nuisance = nuis[[1]],
  scrub = c(rep(FALSE, 253-5), rep(TRUE, 5)),
  TR = 0.72,
  brainstructures = c("left", "sub"),
  subROI = c("Amygdala-L", "Caudate-R", "Thalamus-L"),
  surfL=ciftiTools.files()$surf["left"],
  surfR=ciftiTools.files()$surf["right"],
  resamp_res = resamp_res,
  hpf=.01,
  Bayes=TRUE,
  ar_order = 1,
  ar_smooth = 0,
  return_INLA = "trim",
  verbose = TRUE
)

bglm2 <- BayesGLM(
  BOLD = BOLD[[2]],
  design = des[[2]],
  nuisance = nuis[[2]],
  scrub = c(rep(FALSE, 253-5), rep(TRUE, 5)),
  TR = 0.72,
  brainstructures = c("left", "sub"),
  subROI = c("Amygdala-L", "Caudate-R", "Thalamus-L"),
  surfL=ciftiTools.files()$surf["left"],
  surfR=ciftiTools.files()$surf["right"],
  resamp_res = resamp_res,
  hpf=.01,
  Bayes=TRUE,
  ar_order = 1,
  ar_smooth = 0,
  return_INLA = "trim",
  verbose = TRUE
)

z <- BayesGLM2(list(bglm1,  bglm2), nsamp_theta=5, nsamp_beta=10)

bglm3 <- BayesGLM(
  BOLD = BOLD,
  design = des,
  nuisance = nuis,
  scrub = list(
    c(rep(FALSE, 253-5), rep(TRUE, 5)),
    c(rep(FALSE, 253-5), rep(TRUE, 5))
  ),
  TR = 0.72,
  brainstructures = c("left", "sub"),
  subROI = c("Amygdala-L", "Caudate-R", "Thalamus-L"),
  surfL=ciftiTools.files()$surf["left"],
  surfR=ciftiTools.files()$surf["right"],
  resamp_res = resamp_res,
  hpf=.01,
  Bayes=TRUE,
  ar_order = 1,
  ar_smooth = 0,
  return_INLA = "full",
  verbose = TRUE
)

y <- BayesGLM2(list(bglm3, bglm3, bglm3), nsamp_theta=5, nsamp_beta=10)
yy <- BayesGLM2(list(bglm3, bglm3, bglm3, bglm3, bglm3, bglm3), nsamp_theta=5, nsamp_beta=10)

