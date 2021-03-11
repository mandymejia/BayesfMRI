check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("Miscellaneous functions are working", {
  check_wb()

  # Change this to the HCP main folder
  dir_data_test <- '/N/dcwan/projects/hcp'
  dir_archive_retest <- file.path(dir_data_test, "retest")
  # Change this to a folder to write the retest data in
  dir_data_retest <- '/N/u/ddpham/Carbonate/Desktop/Data/HCP_Retest'
  subject <- "151526"
  gambling_dir <- file.path(subject, "MNINonLinear/Results/tfMRI_GAMBLING_LR")
  data_zip <- file.path(dir_archive_retest, paste0(subject, "_3T_tfMRI_GAMBLING_preproc.zip"))
  #
  write_dir <- ".."

  fnames <- list(
    cifti = "tfMRI_GAMBLING_LR_Atlas.dtseries.nii",
    e_loss = "EVs/loss_event.txt",
    e_win = "EVs/win_event.txt",
    e_neut = "EVs/neut_event.txt",
    rps = "Movement_Regressors.txt"
  )
  fnames <- lapply(fnames, function(x){file.path(gambling_dir, x)})
  # Load retest data
  for (fname in fnames) {
    if (!file.exists(file.path(dir_data_retest, fname))) {
      cmd <- paste("unzip", data_zip, fname, "-d", dir_data_retest)
      system(cmd)
      stopifnot(file.exists(file.path(dir_data_retest, fname)))
    }
  }
  fnames <- lapply(
    fnames,
    function(x){
      list(test=file.path(dir_data_test, x), retest=file.path(dir_data_retest, x))
    }
  )

  # cii <- lapply(fnames$cifti, function(x){read_xift(x)})
  events <- lapply(fnames[c("e_loss", "e_win", "e_neut")], function(x){lapply(x, function(y){read.table(y, header=FALSE)})})
  RPs <- lapply(fnames$rps, function(x){as.matrix(read.table(x, header=FALSE))})

  surfL_fname <- ciftiTools::demo_files()$surf["left"]
  surfR_fname <- ciftiTools::demo_files()$surf["right"]

  events <- list(
    test = lapply(events, function(x){x$test[,seq(2)]}),
    retest = lapply(events, function(x){x$retest[,seq(2)]})
  )
  events_all <- lapply(events, function(x){do.call(rbind, x)})
  events_all <- lapply(events_all, function(x){x[order(x[,1]),]})

  TR <- .72

  #Test cases:

  #Single session
  # - No prewhitening
  #   - single task -- DONE
  #   - three tasks -- DONE
  # - With prewhitening
  #   - single task - DONE
  #   - three tasks

  #Two sessions, no averaging
  # - No prewhitening
  #   - single task
  #   - three tasks
  # - With prewhitening
  #   - single task
  #   - three tasks

  #Two sessions, with averaging
  # - No prewhitening
  #   - single task
  #   - three tasks
  # - With prewhitening
  #   - single task
  #   - three tasks

  # Default parameters
  system.time(bglm <- BayesGLM_cifti(
    cifti_fname = fnames$cifti$test,
    surfL_fname = surfL_fname,
    surfR_fname = surfR_fname,
    brainstructures = c("left", "right"),
    onsets = events$test,
    TR = TR,
    nuisance = RPs$test,
    resamp_res = 5000,
    verbose = TRUE,
    return_INLA_result = TRUE,
    avg_sessions = FALSE
  ))
  saveRDS(bglm, file=file.path(write_dir, 'BayesGLM_test1.RDS'))

  # Different parameters
  system.time(bglm <- BayesGLM_cifti(
    cifti_fname = fnames$cifti$test,
    surfR_fname = surfR_fname,
    brainstructures = "right",
    onsets = events$test,
    TR = TR,
    nuisance = RPs$test,
    GLM_method = "classical",
    ar_order = 5,
    ar_smooth = 4,
    resamp_res = 6000,
    num.threads=2,
    verbose = TRUE,
    return_INLA_result = TRUE,
    avg_sessions = TRUE
  ))
  saveRDS(bglm, file=file.path(write_dir, 'BayesGLM_test2.RDS'))

  # ...
  system.time(bglm <- BayesGLM_cifti(
    cifti_fname = as.character(fnames$cifti),
    surfL_fname = surfL_fname,
    brainstructures = "left",
    onsets = events, # events_all: Error in onsets[[k]][, 1] : incorrect number of dimensions
    TR = TR,
    nuisance = RPs,
    resamp_res = 5000,
    verbose = TRUE,
    return_INLA_result = TRUE,
    avg_sessions = FALSE
  ))
  saveRDS(bglm, file=file.path(write_dir, 'BayesGLM_test3.RDS'))

  # Combining to one task would be faster

  #thresold = 0 means looking for activations > 0
  #threshold represents percent signal change. Could set threshold=0.5 or 1 as a biologically meaningful level.
  # system.time(act3_Bayes <- id_activations_cifti(
  #   results_pw_three, threshold=0, method='Bayesian', alpha=0.05)
  # )
  # system.time(act3_classical1 <- id_activations_cifti(
  #   results_pw_three, threshold=0, method='classical', correction='FWER', alpha=0.05)
  # )
  # system.time(act3_classical2 <- id_activations_cifti(
  #   results_pw_three, threshold=0, method='classical', correction='FDR', alpha=0.05)
  # )
  #save(act3_Bayes, act3_classical1, act3_classical2, file=file.path(main_dir, 'act_bayesianModel3.RData'))
#
#   #visualize Bayesian estimates
#   view_xifti_surface(results_pw_three$betas_Bayesian$single_session,
#                     surfL = surfL_fname, zlim=c(-1.2, 1.2), title='Bayesian', idx=1)
#
#   #visualize Bayesian estimates and save to PNG
#   view_xifti_surface(results_pw_three$betas_Bayesian$single_session,
#                     surfL = surfL_fname, zlim=c(-1.2, 1.2), title='Bayesian', idx=1,
#                     fname = 'Bayesian_est')
#
#   #visualize classical estimates
#   view_xifti_surface(results_pw_three$betas_classical$single_session, surfL = surfL_fname, zlim=c(-1.2, 1.2), title='classical', idx=1)
#
#
#   view_xifti_surface(act_classical$activations_xifti, surfL = surfL_fname, title='classical (FWER)', idx=1)
#   view_xifti_surface(act_classical1$activations_xifti, surfL = surfL_fname, title='classical (FDR)', idx=1)
#   view_xifti_surface(act_Bayes$activations_xifti, surfL = surfL_fname, title='Bayesian', idx=1)
#
#   #threshold estimates and visualize
#   estBayes_thr <- (results_pw_one$betas_Bayesian$single_session)*(act_Bayes$activations_xifti)
#   estBayes_thr <- transform_xifti(estBayes_thr, function(x) {x[x==0] <- NA; return(x)})
#   view_xifti_surface(estBayes_thr, surfL = surfL_fname, title='Bayesian estimates thresholded', idx=1, zlim=c(0,1.2))
#
#
#   #### RETEST
#
#   system.time(results_pw_one_retest <- BayesGLM_cifti(
#     cifti_fname = cifti_fname2,
#     surfL_fname = surfL_fname,
#     surfR_fname = surfR_fname,
#     brainstructures = 'left',
#     #brainstructures = c("left", "right"),
#     design = NULL,
#     onsets = events_all2,
#     #onsets = events_separate2,  #uncomment to analyze win/loss/neutral separately
#     TR = TR,
#     nuisance = nuisance2,
#     nuisance_include = c("drift", "dHRF"),
#     scale_BOLD = TRUE,
#     scale_design = TRUE,
#     GLM_method = "both",
#     ar_order = 6,
#     ar_smooth = 5,
#     resamp_res = 5000,
#     num.threads = 4,
#     verbose = TRUE,
#     outfile = NULL,
#     return_INLA_result = TRUE,
#     avg_sessions = FALSE,
#     trim_INLA = TRUE
#   ))

  #save(results_pw_one_retest, file=file.path(main_dir, 'bayesianModel_retest.RData'))
  #load(file=file.path(main_dir, 'bayesianModel_retest.RData'))

})
