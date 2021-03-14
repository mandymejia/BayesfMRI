check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("Miscellaneous functions are working", {
  check_wb()

  # Folder for the data to use for testing BayesfMRI
  if (!exists("dir_data")) { dir_data <- "../BayesfMRI_testData" }

  # Get file names.
  fnames <- list(
    cifti = "tfMRI_GAMBLING_LR_Atlas.dtseries.nii",
    e_loss = "EVs/loss_event.txt",
    e_win = "EVs/win_event.txt",
    e_neut = "EVs/neut_event.txt",
    rps = "Movement_Regressors.txt"
  )
  fnames <- lapply(fnames,
    function(x){
      list(
        test=file.path(dir_data, "test", x), 
        retest=file.path(dir_data, "retest", x)
      )
    }
  )

  # Read files.
  events <- lapply(fnames[c("e_loss", "e_win", "e_neut")], 
    function(x){
      lapply(x, function(y){ read.table(y, header=FALSE) })
    }
  )
  RPs <- lapply(fnames$rps, 
    function(x){
      as.matrix(read.table(x, header=FALSE))
    }
  )
  surfL_fname <- ciftiTools::demo_files()$surf["left"]
  surfR_fname <- ciftiTools::demo_files()$surf["right"]
  events <- list(
    test = lapply(events, function(x){x$test[,seq(2)]}),
    retest = lapply(events, function(x){x$retest[,seq(2)]})
  )
  # Merge the 3 gambling events.
  events_all <- lapply(events, function(x){do.call(rbind, x)})
  events_all <- lapply(events_all, function(x){list(gamble=x[order(x[,1]),])})

  TR <- .72

  # data.frame where each row describes a combination of arguments to test.
  sess_avg <- c("1_TRUE", "2_TRUE", "2_FALSE")
  ntask <- c(1,3)
  prewhiten <- c(TRUE, FALSE)
  params <- expand.grid(sess_avg=sess_avg, ntask=ntask, prewhiten=prewhiten, stringsAsFactors=FALSE)
  params$sess <- as.numeric(gsub("_.*", "", sess_avg))
  params$avg <- as.logical(gsub(".*_", "", sess_avg))
  params$bs <- "left"
  params$smooth <- 5
  # Add non-combinatorial test: both hemispheres
  params <- rbind(params, c("1_TRUE", 1, TRUE, 1, TRUE, "both", 5))
  # Add non-combinatorial test: no AR smoothing
  params <- rbind(params, c("1_TRUE", 1, TRUE, 1, TRUE, "both", 0))
  # Add non-combinatorial test: no nuisance regressors
  params <- rbind(params, c("1_TRUE", 1, TRUE, 1, TRUE, "both", 5))

  # Test each combination.
  for (ii in seq(nrow(params))) {
    # Print a description of the combination to test.
    cat(
      params$sess[ii],
      ifelse(
        params$sess[ii]>1,
        ifelse(params$avg[ii], "sessions (w/ averaging),", "sessions (not avg),"),
        "session,"
      ),
      ifelse(params$ntask[ii]==1, "1 task,", "3 tasks,"), "and",
      ifelse(params$prewhiten[ii], "with prewhitening", "without prewhitening"), "\n\n"
    )

    # Format the arguments
    if (params$ntask[ii] == 3) {
      onsets_ii <- events
    } else {
      onsets_ii <- events_all
    }
    cii_ii <- fnames$cifti
    nuisance_ii <- RPs
    if (params$sess[ii] == 1) {
      onsets_ii <- onsets_ii$test
      cii_ii <- cii_ii$test
      nuisance_ii <- nuisance_ii$test
    } else {
      cii_ii <- as.character(cii_ii)
    }
    if (params$bs[ii] == "left") {
      bs_ii <- "left"
    } else {
      bs_ii <- c("left", "right")
    }

    # Do BayesGLM_cifti.
    # --- left hemisphere only
    # --- AR_smooth is 5
    # --- resolution is 2000
    exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
      cifti_fname = cii_ii,
      surfL_fname = surfL_fname,
      surfR_fname = surfR_fname,
      brainstructures = bs_ii,
      onsets = onsets_ii,
      TR = TR,
      nuisance = nuisance_ii,
      GLM_method = "classical",
      ar_order = ifelse(params$prewhiten[ii], 6, 0),
      ar_smooth = params$smooth[ii],
      resamp_res = 2000,
      verbose = FALSE,
      return_INLA_result = TRUE,
      outfile = file.path(dir_write, "bfmri_out"),
      avg_sessions = params$avg[ii]
    ))

    # Identify the activations.
    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='Bayesian', alpha=0.05)
    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='classical', correction='FWER', alpha=0.05)

    # Save/plot?

    print(exec_time)
  }
})
