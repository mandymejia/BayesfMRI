check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("Miscellaneous functions are working", {
  check_wb()

  #main_dir <- "/Users/analysis/Desktop/ASINag/BayesfMRILab"
  main_dir <- '~/Google Drive/My Drive/R25TrainingCourse/Year2021/NeuroStatsBootCampLabs/W1D04-Stats-2/Lab-03/'

  Subject <- "151526"
  data_dir <- file.path(main_dir,Subject)
  data_dir_retest <- file.path(main_dir,paste0(Subject,"_retest"))

  #surfL_fname <-file.path(data_dir,"151526.L.midthickness.32k_fs_LR.surf.gii")
  #surfR_fname <-file.path(data_dir,"151526.R.midthickness.32k_fs_LR.surf.gii")
  surfL_fname <- file.path(main_dir,'../Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii')
  surfR_fname <- file.path(main_dir,'../Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii')

  cifti_fname <-file.path(data_dir,"tfMRI_GAMBLING_LR_Atlas.dtseries.nii")
  cifti_fname2 <-file.path(data_dir_retest,"tfMRI_GAMBLING_LR_Atlas.dtseries.nii")

  nuisance <- as.matrix(read.table(file.path(data_dir,"Movement_Regressors.txt"),header = FALSE))
  nuisance2 <- as.matrix(read.table(file.path(data_dir_retest,"Movement_Regressors.txt"),header = FALSE))

  events_loss <- read.table(file.path(data_dir,"loss_event.txt"),header = FALSE)
  events_win <- read.table(file.path(data_dir,"win_event.txt"),header = FALSE)
  events_neut <- read.table(file.path(data_dir,"neut_event.txt"),header = FALSE)
  events_loss2 <- read.table(file.path(data_dir_retest,"loss_event.txt"),header = FALSE)
  events_win2 <- read.table(file.path(data_dir_retest,"win_event.txt"),header = FALSE)
  events_neut2 <- read.table(file.path(data_dir_retest,"neut_event.txt"),header = FALSE)

  #for modeling all events as a single task
  events_all <- rbind(events_loss,events_win,events_neut)[,1:2]
  theOrder <- order(events_all[,1])
  events_all <- list(task=events_all[theOrder,])
  events_all2 <- rbind(events_loss,events_win,events_neut)[,1:2]
  theOrder2 <- order(events_all2[,1])
  events_all2 <- list(task=events_all2[theOrder2,])

  #for modeling loss, win and neutral events separately
  events_separate <- list(loss=events_loss[,1:2], win=events_win[,1:2], neut=events_neut[,1:2])
  events_separate2 <- list(loss=events_loss2[,1:2], win=events_win2[,1:2], neut=events_neut2[,1:2])

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

  #for class
  system.time(results_pw_three <- BayesGLM_cifti(
    cifti_fname = cifti_fname,
    surfL_fname = surfL_fname,
    surfR_fname = surfR_fname,
    #brainstructures = 'left',
    brainstructures = c("left", "right"),
    design = NULL,
    #onsets = events_all,
    onsets = events_separate,  #uncomment to analyze win/loss/neutral separately
    TR = TR,
    nuisance = nuisance,
    nuisance_include = c("drift", "dHRF"),
    scale_BOLD = TRUE,
    scale_design = TRUE,
    GLM_method = "both",
    #ar_order = 6,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 5000,
    num.threads = 4,
    verbose = TRUE,
    outfile = NULL,
    return_INLA_result = TRUE,
    avg_sessions = FALSE,
    trim_INLA = TRUE
  ))
  save(results_pw_three, file=file.path(main_dir, 'bayesianModel3.RData'))

  #thresold = 0 means looking for activations > 0
  #threshold represents percent signal change. Could set threshold=0.5 or 1 as a biologically meaningful level.
  system.time(act3_Bayes <- id_activations_cifti(results_pw_three, threshold=0, method='Bayesian', alpha=0.05))
  system.time(act3_classical1 <- id_activations_cifti(results_pw_three, threshold=0, method='classical', correction='FWER', alpha=0.05))
  system.time(act3_classical2 <- id_activations_cifti(results_pw_three, threshold=0, method='classical', correction='FDR', alpha=0.05))
  save(act3_Bayes, act3_classical1, act3_classical2, file=file.path(main_dir, 'act_bayesianModel3.RData'))

  system.time(results_pw_three_retest <- BayesGLM_cifti(
    cifti_fname = cifti_fname2,
    surfL_fname = surfL_fname,
    surfR_fname = surfR_fname,
    #brainstructures = 'left',
    brainstructures = c("left", "right"),
    design = NULL,
    #onsets = events_all,
    onsets = events_separate2,  #uncomment to analyze win/loss/neutral separately
    TR = TR,
    nuisance = nuisance2,
    nuisance_include = c("drift", "dHRF"),
    scale_BOLD = TRUE,
    scale_design = TRUE,
    GLM_method = "both",
    #ar_order = 6,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 5000,
    num.threads = 4,
    verbose = FALSE,
    outfile = NULL,
    return_INLA_result = TRUE,
    avg_sessions = FALSE,
    trim_INLA = TRUE
  ))
  save(results_pw_three_retest, file=file.path(main_dir, 'bayesianModel3_retest.RData'))

  system.time(act3_Bayes_retest <- id_activations_cifti(results_pw_three_retest, threshold=0, method='Bayesian', alpha=0.05))
  system.time(act3_classical1_retest <- id_activations_cifti(results_pw_three_retest, threshold=0, method='classical', correction='FWER', alpha=0.05))
  system.time(act3_classical2_retest <- id_activations_cifti(results_pw_three_retest, threshold=0, method='classical', correction='FDR', alpha=0.05))
  save(act3_Bayes_retest, act3_classical1_retest, act3_classical2_retest, file=file.path(main_dir, 'act_bayesianModel3_retest.RData'))


  #visualize Bayesian estimates
  view_xifti_surface(results_pw_three$betas_Bayesian$single_session,
                    surfL = surfL_fname, zlim=c(-1.2, 1.2), title='Bayesian', idx=1)

  #visualize Bayesian estimates and save to PNG
  view_xifti_surface(results_pw_three$betas_Bayesian$single_session,
                    surfL = surfL_fname, zlim=c(-1.2, 1.2), title='Bayesian', idx=1,
                    fname = 'Bayesian_est')

  #visualize classical estimates
  view_xifti_surface(results_pw_three$betas_classical$single_session, surfL = surfL_fname, zlim=c(-1.2, 1.2), title='classical', idx=1)


  view_xifti_surface(act_classical$activations_xifti, surfL = surfL_fname, title='classical (FWER)', idx=1)
  view_xifti_surface(act_classical1$activations_xifti, surfL = surfL_fname, title='classical (FDR)', idx=1)
  view_xifti_surface(act_Bayes$activations_xifti, surfL = surfL_fname, title='Bayesian', idx=1)

  #threshold estimates and visualize
  estBayes_thr <- (results_pw_one$betas_Bayesian$single_session)*(act_Bayes$activations_xifti)
  estBayes_thr <- transform_xifti(estBayes_thr, function(x) {x[x==0] <- NA; return(x)})
  view_xifti_surface(estBayes_thr, surfL = surfL_fname, title='Bayesian estimates thresholded', idx=1, zlim=c(0,1.2))


  #### RETEST

  system.time(results_pw_one_retest <- BayesGLM_cifti(
    cifti_fname = cifti_fname2,
    surfL_fname = surfL_fname,
    surfR_fname = surfR_fname,
    brainstructures = 'left',
    #brainstructures = c("left", "right"),
    design = NULL,
    onsets = events_all2,
    #onsets = events_separate2,  #uncomment to analyze win/loss/neutral separately
    TR = TR,
    nuisance = nuisance2,
    nuisance_include = c("drift", "dHRF"),
    scale_BOLD = TRUE,
    scale_design = TRUE,
    GLM_method = "both",
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 5000,
    num.threads = 4,
    verbose = TRUE,
    outfile = NULL,
    return_INLA_result = TRUE,
    avg_sessions = FALSE,
    trim_INLA = TRUE
  ))

  save(results_pw_one_retest, file=file.path(main_dir, 'bayesianModel_retest.RData'))
  #load(file=file.path(main_dir, 'bayesianModel_retest.RData'))

})