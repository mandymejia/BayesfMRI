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
  # Change this to a folder to write the BayesGLM_cifti outfile
  dir_write <- "../DeleteLater"

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
  events_all <- lapply(events_all, function(x){list(gamble=x[order(x[,1]),])})

  TR <- .72

  sess_avg <- c("1_TRUE", "2_TRUE", "2_FALSE")
  ntask <- c(1,3)
  prewhiten <- c(TRUE, FALSE)
  params <- expand.grid(sess_avg=sess_avg, ntask=ntask, prewhiten=prewhiten, stringsAsFactors=FALSE)
  params$sess <- as.numeric(gsub("_.*", "", sess_avg))
  params$avg <- as.logical(gsub(".*_", "", sess_avg))

  for (ii in seq(nrow(params))) {
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

    exec_time <- system.time(bfmri_ii <- BayesGLM_cifti(
      cifti_fname = cii_ii,
      surfL_fname = surfL_fname,
      surfR_fname = surfR_fname,
      brainstructures = "left",
      onsets = onsets_ii,
      TR = TR,
      nuisance = nuisance_ii,
      ar_order = ifelse(params$prewhiten[ii], 6, 0),
      ar_smooth = 5,
      resamp_res = 2000,
      verbose = FALSE,
      return_INLA_result = TRUE,
      outfile = file.path(dir_write, "bfmri_out"),
      avg_sessions = params$avg[ii]
    ))

    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='Bayesian', alpha=0.05)
    act_ii <- id_activations_cifti(bfmri_ii, threshold=0, method='classical', correction='FWER', alpha=0.05)

    print(exec_time)
  }
})
