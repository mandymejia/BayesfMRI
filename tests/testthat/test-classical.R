check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("Classical modeling working", {
  check_wb()

  # Folder for the data to use for testing BayesfMRI
  if (!exists("dir_data")) { dir_data <- "data" }
  if (!exists("dir_results")) { dir_results <- file.path(dir_data, "results") }

  # Get file names.
  fnames <- list(
    #cifti = "derivatives surface_pipeline sub-MSC01 processed_task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_bold_32k_fsLR.dtseries.nii",
    #tmask = "derivatives surface_pipeline sub-MSC01 processed_task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_bold_32k_fsLR_tmask.txt"
    cifti_1 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-01_bold.dtseries.nii",
    events_1 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-01_events.tsv",
    cifti_2 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-02_bold.dtseries.nii",
    events_2 = "derivatives surface_pipeline sub-MSC01 task_timecourses ses-func01 sub-MSC01_ses-func01_task-motor_run-02_events.tsv"
  )
  fnames <- lapply(fnames, function(x){file.path(dir_data, x)})

  events <- lapply(fnames[c("events_1", "events_2")], read.table, header=TRUE)
  task_names <- c("RHand", "LHand", "RFoot", "LFoot", "Tongue")
  events <- lapply(events, function(x){
    setNames(lapply(task_names, function(y){
      cbind(
        onset = x$onset[x$trial_type==y],
        duration = x$duration[x$trial_type==y]
      )
    }), task_names)
  })

  # data.frame where each row describes a combination of arguments to test.
  sess_avg <- c("1_TRUE", "2_TRUE", "2_FALSE")
  prewhiten <- c(TRUE, FALSE)
  params <- expand.grid(sess_avg=sess_avg, prewhiten=prewhiten, stringsAsFactors=FALSE)
  params$sess <- as.numeric(gsub("_.*", "", sess_avg))
  params$avg <- as.logical(gsub(".*_", "", sess_avg))
  params$bs <- "left"
  params$smooth <- 5
  # Add non-combinatorial test: both hemispheres; no smoothing
  params <- rbind(
    params,
    data.frame(
      sess_avg="1_TRUE",
      prewhiten=TRUE,
      sess=1,
      avg=TRUE,
      bs=c("both", "right"),
      smooth=c(5, 0)
    )
  )
  # [TO DO] Add non-combinatorial test: no nuisance regressors

  # Test each combination.
  for (ii in seq(nrow(params))) {
    # Print a description of the combination to test.
    cat(
      params$sess[ii],
      ifelse(
        params$sess[ii]>1,
        ifelse(params$avg[ii], "sessions (w/ averaging),", "sessions (not avg),"),
        "session,"
      ), "and",
      ifelse(params$prewhiten[ii], "with prewhitening", "without prewhitening"), "\n\n"
    )

    # Do BayesGLM_cifti.
    # --- left hemisphere only
    # --- AR_smooth is 5
    # --- resolution is 2000
    exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
      cifti_fname = c(fnames$cifti_1, fnames$cifti_2)[seq(params$sess[ii])],
      surfL_fname=ciftiTools.files()$surf["left"],
      surfR_fname=ciftiTools.files()$surf["right"],
      brainstructures = params$bs[ii],
      onsets = switch(params$sess[ii], events[[1]], events[seq(2)]),
      TR = 2.2,
      Bayes = FALSE,
      ar_order = ifelse(params$prewhiten[ii], 6, 0),
      ar_smooth = params$smooth[ii],
      resamp_res = 1000,
      verbose = FALSE,
      return_INLA_result = TRUE,
      outfile = file.path(dir_results, "bfmri_out"),
      avg_sessions = params$avg[ii]
    ))

    # Identify the activations.
    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='classical', alpha=0.05)
    # plot(act_ii$activations_xifti)
    if (ii == 1) {
      # Test the other arguments too.
      act_ii <- id_activations_cifti(
        bfmri_ii, threshold=.1, method='classical', alpha=0.1, correction='FWER'
        #, excur_method='QC'
      )
    }

    # Save/plot?

    print(exec_time)
  }
})
