# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 7000
# my_pardiso <- "~/Documents/pardiso.lic" # INLA PARDISO license
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
library(stats)
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
des <- make_design(
  onsets = list(events[seq(3)], events[seq(4,6)]),
  TR = 0.72, nTime=253
)

### Classical vs Bayes; Single- vs Multi-session -----
BayesGLM_cifti_args <- function(n_sess, resamp_factor=1, dtype=c("single", "multi")){
  dtype <- match.arg(dtype, c("single", "multi"))
  if (dtype == "multi") {
    des$design <- lapply(des$design, function(ddd){
      abind::abind(
        ddd,
        ddd + rnorm(prod(dim(ddd)), sd=.1),
        ddd + rnorm(prod(dim(ddd)), sd=.5), along=3
      )
    })
  }
  des1 <- des; des1$design <- des1$design[1]
  list(
    cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(n_sess)],
    surfL_fname=ciftiTools.files()$surf["left"],
    surfR_fname=ciftiTools.files()$surf["right"],
    brainstructures = "both",
    design=switch(n_sess, des1, des),
    TR = 0.72,
    dHRF=1,
    nuisance=switch(n_sess, nuis$rp_1, nuis),
    #Bayes = TRUE,
    ar_order = 1,
    ar_smooth = 3,
    resamp_res = resamp_res * resamp_factor,
    verbose = TRUE,
    return_INLA = "trim"
  )
}

##### First pass to detect errors
bglm_c1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1, resamp_factor=.1)))
bglm_b1 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(1, resamp_factor=.1)))
bglm_m1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1, resamp_factor=.1, dtype="multi")))
bglm_c2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2, resamp_factor=.1)))
bglm_b2 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(2, resamp_factor=.1)))
bglm_m2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2, resamp_factor=.1, dtype="multi")))


##### Second pass to get results of decent resolution
bglm_c1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1)))
bglm_b1 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(1)))
bglm_m1 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(1, dtype="multi")))
bglm_c2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2)))
bglm_b2 <- do.call(BayesGLM_cifti, c(list(Bayes=TRUE), BayesGLM_cifti_args(2)))
bglm_m2 <- do.call(BayesGLM_cifti, c(list(Bayes=FALSE), BayesGLM_cifti_args(2, dtype="multi")))

act_c1 <- id_activations(bglm_c1, sessions=1)
act_b1 <- id_activations(bglm_b1, gamma=.01, sessions=1)
act_m1 <- id_activations(bglm_b1, gamma=.01, sessions=1)
act_c2 <- id_activations(bglm_c2, gamma=.01, sessions=seq(2))
act_b2 <- id_activations(bglm_b2, gamma=.01, sessions=seq(2))
act_m2 <- id_activations(bglm_b1, gamma=.01, sessions=1)

### Misc. cases; not checking these results, but checking for errors
### Last updated: 5.2
des1 <- make_design(onsets = list(events[seq(3)]), TR=.72, nTime=253)
bglm_x1 <- BayesGLM_cifti(
  cifti_fname = fnames$cifti_1,
  brainstructures = "right",
  design = des1,
  TR = 0.72,
  dHRF_as="field",
  Bayes = TRUE,
  resamp_res=resamp_res,
  ar_order = 0,
  verbose = 0
)
act_m1 <- id_activations(bglm_x1, alpha=.1, gamma=.05, fields=1)

bglm_x2 <- BayesGLM_cifti(
  cifti_fname = c(fnames$cifti_1, fnames$cifti_2),
  brainstructures = "left",
  design=des,
  TR = 0.72,
  dHRF=0,
  Bayes = FALSE,
  resamp_res=resamp_res*2,
  nuisance=nuis,
  ar_order = 2,
  ar_smooth = 0,
  verbose = 2
)
act_m2 <- id_activations(bglm_x2)

# Save ---
save(list=ls(), file=file.path(dir_resultThis, "test_notInPackage_HCP.rda"))

# Plot ---
plot(bglm_c1, fname=file.path(dir_resultThis, "HCP_bglm_c1"), together="idx")
plot(bglm_b1, fname=file.path(dir_resultThis, "HCP_bglm_b1"), together="idx")
plot(bglm_m1, fname=file.path(dir_resultThis, "HCP_bglm_m1"), together="idx")
plot(bglm_c2, fname=file.path(dir_resultThis, "HCP_bglm_c2"), together="idx")
plot(bglm_b2, fname=file.path(dir_resultThis, "HCP_bglm_b2"), together="idx")
plot(bglm_m2, fname=file.path(dir_resultThis, "HCP_bglm_m2"), together="idx")

plot(act_c1, fname=file.path(dir_resultThis, "HCP_act_c1"), together="idx")
plot(act_b1, fname=file.path(dir_resultThis, "HCP_act_b1"), together="idx")
plot(act_m1, fname=file.path(dir_resultThis, "HCP_act_m1"), together="idx")
plot(act_c2, fname=file.path(dir_resultThis, "HCP_act_c2"), together="idx")
plot(act_b2, fname=file.path(dir_resultThis, "HCP_act_b2"), together="idx")
plot(act_m2, fname=file.path(dir_resultThis, "HCP_act_m2"), together="idx")

plot(bglm_x1, fname=file.path(dir_resultThis, "HCP_bglm_x1"), together="idx")
plot(bglm_b2, fname=file.path(dir_resultThis, "HCP_bglm_x2"), together="idx")

plot(act_m1, fname=file.path(dir_resultThis, "HCP_act_m1"), together="idx")
plot(act_m2, fname=file.path(dir_resultThis, "HCP_act_m2"), together="idx")
