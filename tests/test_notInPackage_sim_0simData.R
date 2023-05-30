# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
resamp_res <- 9000
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_results <- "tests/results_notInPackage"
dir_resultsThis <- file.path(dir_results, "simData")

library(testthat)
if (doINLA) { library(INLA) }
library(fMRItools)
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)

# Simulate data ----------------------------------------------------------------
cat("Simulating data\n~~~~~~~~~~~~~~~~\n")
source("tests/simulate_cifti.R")
set.seed(0)
bsim <- simulate_cifti_multiple(
  wb_path=my_wb,
  brainstructures="both",
  n_subjects=3,
  n_sessions=2,
  n_runs=1,
  ntasks=3,
  ntime=145,
  subject_var=0,
  session_var=0,
  run_var=0,
  ar_error=c(.1, .08, .03),
  resamp_res=resamp_res
)
for (ss in seq(length(bsim$simulated_cifti))) {
  bsim$simulated_cifti[[ss]] <- add_noise(bsim$simulated_cifti[[ss]], sd=50)
}
saveRDS(bsim, file.path(dir_resultsThis, "bsim.rds"))

# Plot data --------------------------------------------------------------------
# bsim <- readRDS(file.path(dir_resultsThis, "bsim.rds"))
for (ss in seq(length(bsim$simulated_cifti))) {
  bsim_ss_name <- names(bsim$simulated_cifti)[ss]
  plot(
    bsim$coef_cifti[[ss]],
    idx=seq(3), together="idx", together_ncol=1,
    together_title=bsim_ss_name, zlim=c(0, 1),
    legend_fname=NULL, legend_embed=FALSE,
    fname=file.path(dir_resultsThis, bsim_ss_name)
  )
}
view_comp(
  img = file.path(
    dir_resultsThis,
    paste0(names(bsim$simulated_cifti[seq(3)]), ".png")
  ),
  nrow=1, title="Comp",
  fname = file.path(dir_resultsThis, "comp.png")
)
