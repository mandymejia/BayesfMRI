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

# Make version of `events` with no stimuli for some tasks, for testing.
eventsB <- events
eventsB[c(2, 6)] <- NA

# Test `make_design` -----------------------------------------------------------

des <- make_design(events[1], TR=.72, nTime=nTime,dHRF=0)

# # quick little test w/ the other options
# win_len <- c(0, .1, .3, .5, 1, 3, 5, 10)
# pdf("~/Desktop/win_duration.pdf")
# for (ww in seq(length(win_len))) {
#   events$win[,2] <- win_len[ww]
#   des <- make_design(
#     events[seq(3)], TR=.72, nTime=nTime, scale_design = FALSE,
#     dHRF=0, onset="win", offset=c("win")
#   )
#   print(matplot(
#     des$design, type="l", col=c("blue", "pink", "#ffcc99", "lightblue", "darkblue"),
#     main=paste0("win duration: ", win_len[ww]), ylim=c(-100, 700),
#     lty=1, lwd=1.5
#   ))
# }
# dev.off()
# pdf("~/Desktop/win_duration_scaled.pdf")
# for (ww in seq(length(win_len))) {
#   events$win[,2] <- win_len[ww]
#   des <- make_design(
#     events[seq(3)], TR=.72, nTime=nTime,
#     dHRF=0, onset="win", offset=c("win")
#   )
#   print(matplot(
#     des$design, type="l", col=c("blue", "pink", "#ffcc99", "lightblue", "darkblue"),
#     main=paste0("win duration: ", win_len[ww]), ylim=c(-.1, .7)/.7,
#     lty=1, lwd=1.5
#   ))
# }
# dev.off()

des <- make_design(events[seq(3)], TR=.72, nTime=nTime, dHRF=2, onset="win", offset=c("loss", "neut"))
des <- make_design(events[seq(3)], TR=.72, nTime=nTime, dHRF=0)

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
### Classical vs Bayes; Single- vs Multi-session -----
BGLM_cii_args <- function(sess, resamp_factor=1){
  BOLD_ss <- c(fnames$cifti_1, fnames$cifti_2)[sess]
  design_ss <- des[sess]
  nuis_ss <- nuis[sess]
  if (length(sess)==1) {
    design_ss <- design_ss[[1]]
    nuis_ss <- nuis_ss[[1]]
  }

  list(
    BOLD = BOLD_ss,
    design = design_ss,
    nuisance = nuis_ss,
    TR = 0.72,
    brainstructures = c("left", "sub"),
    surfL=ciftiTools.files()$surf["left"],
    surfR=ciftiTools.files()$surf["right"],
    resamp_res = resamp_res * resamp_factor,
    hpf=.01,
    ar_order = 1,
    ar_smooth = 3,
    return_INLA = "trim",
    verbose = TRUE
  )
}

# ##### First pass to detect errors
#bglm_c1 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(1, resamp_factor=.1)))
#bglm_c2 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(2, resamp_factor=.1)))
bglm_b1 <- do.call(BayesGLM, c(list(Bayes=TRUE), BGLM_cii_args(1, resamp_factor=.1)))
bglm_b2 <- do.call(BayesGLM, c(list(Bayes=TRUE), BGLM_cii_args(2, resamp_factor=.1)))

bglm_x1 <- bglm_b1$BGLMs$cortexL
bglm_x2 <- bglm_b2$BGLMs$cortexL
bglm2a <- BayesGLM2(list(bglm_x1, bglm_x2))
bglm_x1 <- bglm_b1; bglm_x1$BGLMs$subcort <- NULL
bglm_x2 <- bglm_b2; bglm_x2$BGLMs$subcort <- NULL
bglm2b <- BayesGLM2(list(bglm_x1, bglm_x2))

bglm_x1 <- bglm_b1$BGLMs$subcort
bglm_x2 <- bglm_b2$BGLMs$subcort
bglm2c <- BayesGLM2(list(bglm_x1, bglm_x2))

# Misc
BOLD <- as.matrix(read_cifti(fnames$cifti_1))
design <- abind::abind(des[[1]], des[[2]], along=3)
nuisance <- nuis$rp_1

### multiGLM
multiGLM(
  BOLD = BOLD,
  design = design,
  nuisance=cbind(nuisance, dct_bases(253, 5)),
  verbose = TRUE
)

### Per-location design
bglmA <- BayesGLM(
  BOLD = read_cifti(fnames$cifti_1),
  design = cbind(des[[1]], 1),
  nuisance=cbind(nuis$rp_1, dct_bases(253, 5), 1), Bayes=FALSE,
  verbose=TRUE, hpf=.01, ar_order=0, resamp_res=100, TR=.72
)
bglmB <- BayesGLM(
  BOLD = read_cifti(fnames$cifti_1),
  design = des[[2]],
  nuisance=cbind(nuis$rp_1, dct_bases(253, 5)), Bayes=FALSE,
  verbose=TRUE, hpf=.01, ar_order=0, resamp_res=100, TR=.72
)
# alternate every voxel A & B
desX <- do.call(cbind, des)
desX <- array(
  c(rep(c(desX), floor(169/2)), desX[,seq(3)]),
  dim=c(253, 3, 169)
)
bglmX <- BayesGLM(
  BOLD = read_cifti(fnames$cifti_1),
  design = desX,
  #nuisance=cbind(nuis$rp_1, dct_bases(253, 5)), Bayes=FALSE,
  verbose=TRUE,# hpf=.01, ar_order=0,
  resamp_res=100, TR=.72
)
## looks good

### Subcortex
BOLD <- read_cifti(fnames$cifti_1, brainstructures="sub")
bglmA <- BayesGLM(
  BOLD = BOLD,
  design = cbind(des[[1]], 1),
  nuisance=cbind(nuis$rp_1, dct_bases(253, 5)), Bayes=FALSE,
  verbose=TRUE, hpf=.01, ar_order=0, TR=.72
)

##### Second pass to get results of decent resolution
bglm_c1 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(1)))
bglm_b1 <- do.call(BayesGLM, c(list(Bayes=TRUE), BGLM_cii_args(1)))
#bglm_m1 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(1, dtype="multi")))
bglm_c2 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(2)))
bglm_b2 <- do.call(BayesGLM, c(list(Bayes=TRUE), BGLM_cii_args(2)))
#bglm_m2 <- do.call(BayesGLM, c(list(Bayes=FALSE), BGLM_cii_args(2, dtype="multi")))

act_c1 <- activations(bglm_c1, sessions=1)
act_b1 <- activations(bglm_b1, gamma=.01, sessions=1)
#act_m1 <- activations(bglm_b1, gamma=.01, sessions=1)
act_c2 <- activations(bglm_c2, gamma=.01, sessions=seq(2))
act_b2 <- activations(bglm_b2, gamma=.01, sessions=seq(2))
#act_m2 <- activations(bglm_b1, gamma=.01, sessions=1)

# Save ---
save(list=ls(), file=file.path(dir_resultThis, "test_notInPackage_HCP.rda"))

# Plot ---
plot(bglm_c1, fname=file.path(dir_resultThis, "HCP_bglm_c1"), together="idx")
plot(bglm_b1, fname=file.path(dir_resultThis, "HCP_bglm_b1"), together="idx")
plot(bglm_c2, fname=file.path(dir_resultThis, "HCP_bglm_c2"), together="idx")
plot(bglm_b2, fname=file.path(dir_resultThis, "HCP_bglm_b2"), together="idx")

plot(act_c1, fname=file.path(dir_resultThis, "HCP_act_c1"), together="idx")
plot(act_b1, fname=file.path(dir_resultThis, "HCP_act_b1"), together="idx")
plot(act_c2, fname=file.path(dir_resultThis, "HCP_act_c2"), together="idx")
plot(act_b2, fname=file.path(dir_resultThis, "HCP_act_b2"), together="idx")
