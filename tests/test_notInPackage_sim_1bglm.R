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
  ntime=155,
  resamp_res=resamp_res
)
add_nze <- function(xii, sd=1){ xii + matrix(rnorm(prod(dim(xii)))*sd, nrow=nrow(xii)) }
bsim$simulated_cifti$subj1_sess1_run1 <- add_nze(bsim$simulated_cifti$subj1_sess1_run1, sd=40)
bsim$simulated_cifti$subj1_sess2_run1 <- add_nze(bsim$simulated_cifti$subj1_sess2_run1, sd=40)

# BayesGLM ---------------------------------------------------------------------
cat("BayesGLM\n~~~~~~~~~~~~~~~~\n")

# combinatorial parameters
Bayes <- c(FALSE, TRUE)
bs <- c("right", "both")
# parameters w/ rising power
params_pwr <- c("sess1-cmbNo-pwYes", "sess2-cmbNo-pwYes", "sess2-cmbYes-pwYes", "sess2-cmbYes-pwNo")
params <- expand.grid(params_pwr=params_pwr, bs=bs, Bayes=Bayes, stringsAsFactors=FALSE)
# add  no smoothing + no PW + cmb for classical
params$smooth <- 5
params <- rbind(
  data.frame(bs="right", Bayes=FALSE, params_pwr="sess2-cmbYes-pwYes", smooth=0),
  params
)
# convert `params_pwr` column to individual params
params$sess <- as.numeric(gsub("sess", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 1)))
params$cmb <- "Yes" == gsub("cmb", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 2))
params$pw <- "Yes" == gsub("pw", "", vapply(strsplit(params$params_pwr, "-"), '[[', "", 3))
params$params_pwr <- NULL

# Test each combination.
for (ii in seq(nrow(params))) {
  print(params[ii,])
  outf_ii <- file.path(dir_resultThis, paste0("bfmri", ii, ".rds"))
  #if (file.exists(outf_ii)) { next }

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
    resamp_res = resamp_res / ((params$bs[ii]=="both") + 1), # / (params$Bayes[ii] + 1)
    verbose = FALSE,
    return_INLA = "trimmed",
    combine_sessions = params$cmb[ii]
  ))

  act_ii <- id_activations(
    bfmri_ii,
    threshold=.1,
    method=ifelse(params$Bayes[ii], "Bayesian", "classical"),
    correction="FDR",
    alpha=0.05
  )

  if (saveResults) {
    saveRDS(
      list(bfmri=bfmri_ii, act=act_ii, exec_time=exec_time),
      outf_ii
    )
  }
}

for (ii in seq(nrow(params))) {
  x_ii <- readRDS(file.path(dir_resultThis, paste0("bfmri", ii, ".rds")))
  pfname_pre <- file.path(dir_resultThis, paste0("bfmri", ii))
  ptitle <- if (ii==1) {
    "sess2-cmbYes-pwYes-noSmooth"
  } else {
    params_pwr[((ii-2) %% 4) + 1]
  }
  ptitle <- gsub("-", ", ", ptitle, fixed=TRUE)
  plot(
    x_ii$bfmri,
    idx=seq(3), together="idx", together_ncol=1, zlim=c(-10, 10),
    together_title=ptitle, legend_fname=NULL, legend_embed=FALSE,
    fname=paste0(pfname_pre, "_est.png")
  )
  plot(
    x_ii$act,
    idx=seq(3), together="idx", together_ncol=1,
    together_title=ptitle, legend_fname=NULL, legend_embed=FALSE,
    fname=paste0(pfname_pre, "_act.png")
  )
}

pcomp <- params[seq(4)*4,seq(2)]
pcomp$Bayes <- ifelse(pcomp$Bayes, "Bayes", "classical")
pcomp_names <- apply(pcomp, 1, paste0, collapse="-")
for (jj in seq(4)) {
  idx <- seq(2,5) + (jj-1)*4
  view_comp(
    img = file.path(dir_resultThis, paste0("bfmri", idx, "_est.png")),
    nrow=1, title=gsub("-", ", ", pcomp_names[jj], fixed=TRUE),
    fname = file.path(dir_resultThis, paste0(pcomp_names[jj], "_est.png"))
  )
  view_comp(
    img = file.path(dir_resultThis, paste0("bfmri", idx, "_act.png")),
    nrow=1, title=gsub("-", ", ", pcomp_names[jj], fixed=TRUE),
    fname = file.path(dir_resultThis, paste0(pcomp_names[jj], "_act.png"))
  )
}
