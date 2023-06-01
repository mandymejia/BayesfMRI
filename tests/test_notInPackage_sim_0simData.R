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
set.seed(1)
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
  ar_error=c(.2, .1, .03),
  resamp_res=resamp_res
)
for (ss in seq(length(bsim$simulated_cifti))) {
  bsim$simulated_cifti[[ss]] <- add_noise(bsim$simulated_cifti[[ss]], sd=80)
}
saveRDS(bsim, file.path(dir_resultsThis, "bsim.rds"))

bsim2 <- simulate_cifti_multiple(
  wb_path=my_wb,
  brainstructures="both",
  n_subjects=3,
  n_sessions=2,
  n_runs=1,
  ntasks=3,
  ntime=145,
  subject_var=350,
  session_var=0,
  ar_error=c(.2, .1, .03),
  resamp_res=resamp_res
)
for (ss in seq(length(bsim2$simulated_cifti))) {
  bsim2$simulated_cifti[[ss]] <- add_noise(bsim2$simulated_cifti[[ss]], sd=80)
}
saveRDS(bsim2, file.path(dir_resultsThis, "bsim_wSpatialVar.rds"))

# Plot data --------------------------------------------------------------------
bsim <- readRDS(file.path(dir_resultsThis, "bsim.rds"))
bsim2 <- readRDS(file.path(dir_resultsThis, "bsim_wSpatialVar.rds"))

for (bb in seq(2)) {
  bsim_bb <- list(bsim, bsim2)[[bb]]
  for (ss in seq(length(bsim_bb$simulated_cifti))) {
    bsim_ss_name <- paste0(names(bsim_bb$simulated_cifti)[ss], "_bsim", bb)
    plot(
      bsim_bb$coef_cifti[[ss]],
      idx=seq(3), together="idx", together_ncol=1,
      together_title=bsim_ss_name, zlim=c(0, 1),
      legend_fname=NULL, legend_embed=FALSE,
      fname=file.path(dir_resultsThis, bsim_ss_name)
    )
  }
  view_comp(
    img = file.path(
      dir_resultsThis,
      paste0(names(bsim_bb$simulated_cifti[seq(3)]), "_bsim", bb, ".png")
    ),
    nrow=1, title="Comp",
    fname = file.path(dir_resultsThis, paste0("comp_", bb, ".png"))
  )

}
