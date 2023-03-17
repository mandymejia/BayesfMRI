check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("The classical GLM is working", {

  check_wb()

  data_res <- 3000

  # Get file names.
  fnames <- list(
    cifti_1 = "motor1-3k-LH.rds",
    events_1 = "motor1-events.tsv",
    cifti_2 = "motor2-3k-LH.rds",
    events_2 = "motor2-events.tsv"
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
  # Add non-combinatorial test: no smoothing
  params <- rbind(
    params,
    data.frame(
      sess_avg="1_TRUE",
      prewhiten=TRUE,
      sess=1,
      avg=TRUE,
      bs="left",
      smooth=0
    )
  )

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

    # Do BayesGLM_cifti
    xii_ii <- lapply(c(fnames$cifti_1, fnames$cifti_2)[seq(params$sess[ii])], readRDS)
    exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
      cifti_fname = xii_ii,
      surfL_fname=ciftiTools.files()$surf["left"],
      surfR_fname=ciftiTools.files()$surf["right"],
      resamp_res=data_res,
      brainstructures = params$bs[ii],
      onsets = switch(params$sess[ii], events[[1]], events[seq(2)]),
      TR = 2.2,
      Bayes = FALSE,
      ar_order = ifelse(params$prewhiten[ii], 6, 0),
      ar_smooth = params$smooth[ii],
      verbose = FALSE,
      return_INLA_result = TRUE,
      outfile = NULL,
      avg_sessions = params$avg[ii]
    ))
    print(exec_time)

    expect_no_error(inherits(bfmri_ii, "BayesGLM_cifti"))

    # Identify the activations.
    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='classical', alpha=0.05)
  }
})
