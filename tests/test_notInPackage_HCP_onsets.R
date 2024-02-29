# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
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
  events_1w = "151526/win_event.txt",
  events_1l = "151526/loss_event.txt",
  events_1n = "151526/neut_event.txt",
  events_2w = "151526_retest/win_event.txt",
  events_2l = "151526_retest/loss_event.txt",
  events_2n = "151526_retest/neut_event.txt"
)
fnames <- lapply(fnames, function(x){file.path(dir_data, x)})

# Make `events`.
events <- lapply(fnames[grepl("events", names(fnames))], read.table, header=FALSE)
events <- lapply(events, function(x){
  x <- x[,seq(2)]
  colnames(x) <- c("onset", "duration")
  x
})
names(events) <- rep(c("win", "loss", "neut"), 2)
# Make version of `events` with no stimuli for some tasks, for testing.
eventsB <- events
eventsB[c(2, 6)] <- NA

nTime <- 253

# Test `make_design` -----------------------------------------------------------

# From `onsets`.
des_o1 <- make_design(
  onsets = events[seq(3)],
  TR = 0.72, nTime=nTime
)
des_o2 <- make_design(
  onsets = list(events[seq(3)], events[seq(4,6)]),
  TR = 0.72, nTime=nTime
)
testthat::expect_equal(des_o1$design[[1]], des_o2$design[[1]])
testthat::expect_equal(des_o1$HRF_info[[1]], des_o2$HRF_info[[1]])
testthat::expect_equal(des_o1$dims[c(2,3,4),], des_o2$dims[c(2,3,4),])
testthat::expect_equal(des_o1$valid_cols[1,], des_o2$valid_cols[1,])
des_o3 <- testthat::expect_warning(make_design(
  onsets = list(eventsB[seq(3)], eventsB[seq(4,6)]),
  TR = 0.72, nTime=c(nTime, 280),
))
testthat::expect_equal(names(des_o2), names(des_o3))
testthat::expect_equal(lapply(des_o2, names), lapply(des_o3, names))
testthat::expect_equal(lapply(des_o2, dim), lapply(des_o3, dim))

# From `design`.
des_d1 <- make_design(des_o1$design)
testthat::expect_equal(des_d1, make_design(des_o1$design[[1]]))
testthat::expect_equal(des_d1$design, des_o1$design)
testthat::expect_equal(des_d1$dims, des_o1$dims)
testthat::expect_equal(des_d1$valid_cols, des_o1$valid_cols)

# From `design_compare`.
des_c1 <- make_design(design_compare=array(rep(c(des_o1$design[[1]], 5)), dim=c(dim(des_o1$design[[1]]), 5)))

# From `design_per_location`.
# [TO DO]

# Summary functions.
des_o3
des_d1
des_c1
