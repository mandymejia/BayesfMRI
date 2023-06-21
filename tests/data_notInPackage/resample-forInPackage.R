# Build --> Install and Restart

# [Edit these]
# ## path to your PARDISO license
# my_pardiso <- "../INLA/pardiso.lic"
## path to your Connectome Workbench
my_wb <- "~/Desktop/workbench"
## path to test data
dir_data <- "data"
dir_data2 <- "data_notInPackage"
## path to results from tests
dir_results <- file.path(dir_data, "results")

library(testthat)
# library(INLA)
# # [Edit this] path to your PARDISO license
# if (interactive()) {
#   inla.setOption(pardiso.license = my_pardiso)
#   inla.pardiso.check()
# }
library(brainSim)
library(ciftiTools)
if (interactive()) { ciftiTools.setOption('wb_path', my_wb) }
library(BayesfMRI)

xii <- read_cifti(file.path("tests", dir_data2,
  "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-01_bold.dtseries.nii"
), resamp_res=3000, brainstructures="left")
saveRDS(xii, file.path("tests", dir_data, "motor1-3k-LH.rds"))
xii <- read_cifti(file.path("tests", dir_data2,
  "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-02_bold.dtseries.nii"
), resamp_res=3000, brainstructures="left")
saveRDS(xii, file.path("tests", dir_data, "motor2-3k-LH.rds"))
