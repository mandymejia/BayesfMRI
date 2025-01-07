#' Make design matrix
#'
#' Make the design matrix for the GLM, from the task information.

#' @param EVs The explanatory variables i.e. the task stimulus information,
#'  from which a design matrix will be constructed. This is a list where each
#'  entry represents a task as a matrix of onsets (first column) and durations
#'  (second column) for each stimuli (each row) of the task, in seconds. List
#'  names should be the task names. \code{nTime} and \code{TR} are required.
#'
#'  An example of a properly-formatted \code{EVs} is:
#'  \code{on_s1 <- list(taskA=cbind(on=c(1,9,17), dr=rep(1,3)),
#'  taskB=cbind(on=c(3,27), dr=rep(5,2)))}.
#'  In this example, there are two tasks: the first has three 1s-long stimuli,
#'  while the second has two 5s-long stimuli.
#'  with \code{on_s2} formatted similarly to \code{on_s1}.
#' @param nTime the number of timepoints (volumes) in the task fMRI data.
#' @param TR the temporal resolution of the data, in seconds.
#' @param dHRF Controls the extent of HRF derivatives modeling.
#'
#'  Set to \code{0} to only model the main HRF regressor (default), and not include its
#'  derivatives; set to \code{1} to model the temporal derivative too;
#'  or, set to \code{2} to model both the temporal and dispersion derivatives.
#'  If \code{dHRF==0}, there is one design column (field) per task. If
#'  \code{dHRF==1}, there are two fields per task. And if \code{dHRF==2}, there
#'  are three fields per task.
#'
#'  If there are several tasks and \code{dHRF>0}, the total number of design
#'  matrix columns may exceed five, which may require large computation times
#'  with INLA. The analysis can be adjusted by modeling the derivatives as
#'  nuisance signals rather than as fields. To do so, move the corresponding
#'  columns from the design matrix to the \code{nuisance} argument for
#'  \code{BayesGLM}.
#' @param upsample Upsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param onset,offset Add task regressors indicating the onset and/or offset of
#'  each event block? Provide the names of the tasks as a character vector. All
#'  onsets (or offsets) across the specified tasks will be represented by one
#'  additional column in the design matrix. The task names must match the names
#'  of \code{EVs}. Can also be \code{"all"} to use all tasks.
#'
#'  Onsets/offset modeling is only compatible with a block design experiment.
#'  An error will be raised if the events in \code{EVs} do not have duration
#'  greater than one second.
#' @param scale_design Scale the columns of the design matrix? Default:
#'  \code{TRUE}.
#' @param onsets_sep,offsets_sep Model the onsets (\code{onsets_sep}) or offsets
#'  (\code{offsets_sep}) separately for each task? Default: \code{FALSE}, to
#'  model all onsets together, or all offsets together, as a single field in the
#'  design.
#' @param verbose Print diagnostic messages? Default: \code{TRUE}.
#' @param ... Additional arguments to \code{\link{HRF_calc}}.
#'
#' @importFrom car vif
#' @importFrom fMRItools is_1
#' @importFrom stats convolve cor lm
#'
#' @return A \code{"BfMRI_design"} object: a list with elements
#'  \describe{
#'    \item{design}{The volumes by fields design matrix. Column names are field names.}
#'    \item{field_names}{The name of each task from the provided onsets.}
#'    \item{dHRF}{The input \code{dHRF} parameter.}
#'    \item{HRF_info}{Additional HRF modeling results.}
#' }
#' @export
make_design <- function(
  EVs, nTime, TR, dHRF=0, upsample=100,
  onset=NULL, offset=NULL,
  scale_design=TRUE,
  onsets_sep=FALSE, offsets_sep=FALSE,
  verbose=TRUE, ...
){

  ortho_block <- FALSE

  # Check inputs. --------------------------------------------------------------
  if(!is.list(EVs)) {
    if(is.data.frame(EVs) | is.matrix(EVs)) EVs <- list(task = EVs)
  }
  stopifnot(is.list(EVs))

  EVs <- lapply(EVs, format_EV)
  EVs_NA <- vapply(EVs, function(q){identical(q, NA)}, FALSE)
  if (any(EVs_NA)) {
    if (all(EVs_NA)) { stop("Not a single event has an onset.") }
    warning("Removing ", sum(EVs_NA), " events with missing data (no onsets).")
    EVs <- EVs[!EVs_NA]
  }
  stopifnot(!any(vapply(EVs, isFALSE, FALSE)))
  stopifnot(fMRItools::is_1(nTime, "numeric"))
  stopifnot(all(nTime > 0))
  stopifnot(all(nTime == round(nTime)))
  stopifnot(fMRItools::is_1(TR, "numeric"))
  stopifnot(TR > 0)
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))
  upsample <- round(TR*upsample)/TR # TR*upsample must be an integer
  stopifnot(fMRItools::is_1(ortho_block, "logical"))
  stopifnot(is_1(onsets_sep, "logical"))
  stopifnot(is_1(offsets_sep, "logical"))

  # In the future, this might be an argument.
  FIR_nSec <- 0
  do_FIR <- FIR_nSec > 0
  if (do_FIR) {
    FIR_nSec <- round(FIR_nSec/TR)*TR # FIR_nSec/TR must be an integer to match indices of upsampled stimulus, must round if not
    FIR_nTimeR <- FIR_nSec/TR # (sec)/(sec/vol) = vol  #number of TRs to model with FIR
  }

  # Add `onset` and `offset` to `EVs`. --------------------------------------
  task_names <- names(EVs) # will be updated
  nJ0 <- length(task_names) # not including `onset` or `offset`

  if (!is.null(onset) || !is.null(offset)) {
    min_dur <- min(do.call(c, lapply(EVs, '[[', "duration")))
    if (min_dur == 0) { warning(
      "`onset` and `offset` are only compatible with block-design experiments, ",
      "but at least one event has duration length zero. Proceeding anyway."
    ) }
  }

  # [TO DO] check design is actually block design

  if (!is.null(onset)) {
    stopifnot(is.character(onset))
    if (length(onset)==1 && onset=="all") { onset <- task_names }
    stopifnot(all(onset %in% task_names))
    onset <- EVs[onset]
    # Same onset, but set new duration at zero
    onset <- lapply(onset, function(q){q$duration <- 0; q})
    if (!onsets_sep) {
      nJ_on <- 1
      onset <- do.call(rbind, onset)
      onset$onset <- sort(onset$onset)
      EVs <- c(EVs, list(onset=onset))
    } else {
      nJ_on <- length(onset)
      names(onset) <- paste0("onset_", names(onset))
      EVs <- c(EVs, onset)
    }
  } else {
    nJ_on <- 0
  }

  if (!is.null(offset)) {
    stopifnot(is.character(offset))
    if (length(offset)==1 && offset=="all") { offset <- task_names }
    stopifnot(all(offset %in% task_names))
    offset <- EVs[offset]
    # New onset is onset+duration; new duration is zero
    offset <- lapply(offset, function(q){q$onset <- q$onset + q$duration; q$duration <- 0; q})
    if (!offsets_sep) {
      nJ_off <- 1
      offset <- do.call(rbind, offset)
      offset$onset <- sort(offset$onset)
      EVs <- c(EVs, list(offset=offset))
    } else {
      nJ_off <- length(offset)
      names(offset) <- paste0("offset_", names(offset))
      EVs <- c(EVs, offset)
    }
  } else {
    nJ_off <- 0
  }

  J_on_idx <- if (is.null(onset)) { NULL } else { seq(nJ0+1, nJ0+nJ_on) }
  J_off_idx <- if (is.null(offset)) { NULL } else { seq(nJ0+nJ_on+1, nJ0+nJ_on+nJ_off) }
  all_onset_onset <- do.call(c, lapply(EVs[J_on_idx], '[[', "onset")) # NULL if none
  all_offset_onset <- do.call(c, lapply(EVs[J_off_idx], '[[', "onset")) # NULL if none

  if ((!is.null(onset)) && (!is.null(offset))) {
    if (any(all_onset_onset %in% all_offset_onset)) {
      warning("At least one `onset` coincides with an `offset`.")
    }
  }

  # Define dimension variables and names again (after `onset` & `offset`).
  task_names <- names(EVs)
  nJ <- length(task_names)

  # Compute HRFs. --------------------------------------------------------------
  # Prepare to upsample the stimulus function.
  nSec <- nTime*TR; # Total time of experiment in seconds
  inds <- seq(TR*upsample, nTime*TR*upsample, by = TR*upsample) # Extract EVs by TR after convolution
  #note: since nTime is an integer, the inds will always be integers, as long as TR*upsample is also an integer
  if(max(abs(inds - round(inds))) > 1e-6) stop('Contact developer: detected non-integer indices for downsampling')
  inds <- round(inds) #to fix tiny deviations from zero

  # For FIR.
  if (do_FIR) {
    inds_FIR <- seq(1, by = TR*upsample, length.out = FIR_nTimeR) #first index for each FIR basis function
    if(max(abs(inds_FIR - round(inds_FIR))) > 1e-6) {
      stop('Contact developer: detected non-integer indices for FIR. ',
        'Set FIR_nSec = 0 to skip FIR basis set calculation.')
    }
  }

  # Canonical HRF to use in convolution as a function of time (in seconds).
  u <- seq(0, 30, by=1/upsample) #go out 30 sec
  HRF <- list(c=HRF_calc(t = u, deriv=0, ...))
  if (dHRF > 0) { HRF$d = HRF_calc(t = u, deriv=1, ...)}
  if (dHRF > 1) { HRF$dd = HRF_calc(t = u, deriv=2, ...)}

  ### Iterate over tasks. ------------------------------------------------------
  nK_block <- (dHRF+1)*nJ0
  nK <- nK_block + nJ_on + nJ_off

  if ((!is.null(onset)) || (!is.null(offset))) {
    if (dHRF>0) { cat("Not making fields for HRF derivatives of onsets and/or",
      "offsets. If these are desired, provide with `EVs`.\n")
    }
  }

  stimulus  <- vector("list", nJ)
  design <- matrix(NA, nrow=nTime, ncol=nK)
  field_names <- vector("character", nK)

  for (jj in seq(nJ)) {
    is_onset_or_offset <- jj > nJ0

    # Get design matrix indices and field names for this task.
    if (is_onset_or_offset) {
      field_idx <- jj + dHRF*nJ0
      field_names[field_idx] <- task_names[jj]
    } else {
      field_idx <- jj + seq(0, dHRF)*nJ0
      field_names[field_idx] <- paste0(
        task_names[jj],
        c("", "_dHRF", "_ddHRF")[seq(dHRF+1)]
      )
    }

    # [NOTE] Instead of inputting a zero column, we now delete the column.
    # # Skip if no stimulus for this task.
    # if (length(EVs[[jj]])==1 && is.na(EVs[[jj]])) { #NA is placeholder for missing tasks
    #   warning(
    #     "For task '", task_names[jj],
    #     "': inputting constant-zero column in design matrix for lack of stimulus. ",
    #     "Proceed with caution."
    #   )
    #   next()
    # }

    # Define `ots` (onset times) and `dur` (durations) for this task.
    ots_jj <- EVs[[jj]][,1]
    dur_jj <- EVs[[jj]][,2]

    ### Round `ots_jj` to TR for event-related designs.
    #ots_jj <- ifelse(dur_jj==0, round(ots_jj/TR)*TR, ots_jj) # [NOTE] Damon commented out (see Slack messages)

    ##### Get `stimulus`. -----------------------------------------------------
    stim_jj <- rep(0, nSec*upsample)

    if (scale_design) {
      # Get the stimulus for just one event, to use for scaling.
      stim_jj_one <- stim_jj
      one_idx <- which.min(dur_jj)
      # [NOTE] If the min-len event occurs close to the end of the session, the peak will be missed.
      if (length(unique(dur_jj)) > 1) {
        warning("Task '", task_names[jj], "' has events with differing ",
          "durations. Using the minimum duration, ", dur_jj[one_idx], ", to ",
          "determine the scaling factor for scaling the design matrix columns."
        )
      }
      start_ii <- round(ots_jj[one_idx]*upsample)
      if (start_ii==0) { start_ii <- 1 }
      end_ii <- round(ots_jj[one_idx]*upsample + dur_jj[one_idx]*upsample)
      stim_jj_one[start_ii:end_ii] <- 1
    }

    for (ii in seq(length(ots_jj))) {
      start_ii <- round(ots_jj[ii]*upsample)
      if (start_ii==0) { start_ii <- 1 }
      end_ii <- round(ots_jj[ii]*upsample + dur_jj[ii]*upsample)
      stim_jj[start_ii:end_ii] <- 1
    }
    stimulus[[jj]] <- c(stim_jj,0)[inds]

    ##### Get `design` by convolving stimulus and HRFs. ------------------------
    # Convolve.
    HRF_jj <- if (is_onset_or_offset) { HRF[1] } else { HRF }

    HRF_conv <- lapply(HRF_jj, function(q){
      convolve(stim_jj, rev(q), type="open")
    })

    # Normalize each HRF by dividing by the one-event HRF max.
    # Note that this occurs prior to downsampling.
    if (scale_design) {
      HRF_one_conv <- lapply(HRF_jj, function(q){
        convolve(stim_jj_one, rev(q), type="open")
      })

      for (hh in seq(length(HRF_conv))) {
        HRF_conv[[hh]] <- HRF_conv[[hh]] / max(HRF_one_conv[[hh]])
      }
    }

    # Downsample and save to `design`.
    design[,field_idx] <- do.call(cbind, HRF_conv)[inds,]

    ##### Get `FIR`. -----------------------------------------------------------
    if (FIR_nSec > 0) {
      for (tt in seq(FIR_nTimeR)) {
        dummy_tt <- rep(0, TR*upsample*FIR_nTimeR)
        start_tt <- inds_FIR[tt]
        end_tt <- start_tt + TR*upsample - 1 #each FIR basis spans a single TR
        dummy_tt[start_tt:end_tt] <- 1
        FIR_tj <- convolve(stim_jj, rev(dummy_tt), type='open')
        FIR_tj <- round(FIR_tj, 5)
        # out$FIR[,jj,tt] <- FIR_tj[inds]
      } #end loop over FIR bases
    } #end FIR basis set construction
  } #end loop over tasks

  # Ortogonalize block fields wrt onset/offset ---------------------------------
  if (ortho_block && nK_block < nK) {
    design[,seq(nK_block)] <- nuisance_regression(
      design[,seq(nK_block)], cbind(1, design[,seq(nK_block+1, nK)])
    )
  }

  # Diagnostics ----------------------------------------------------------------

  if (nK == 1) {
    des_cor_max <- des_vif_max <- tdiff_min <- NA
  } else {
    ### Correlation ------------------------------------------------------------
    des_cor <- checkX_cor(design)
    des_cor_max <- max(abs(des_cor), na.rm=TRUE)
    if (verbose) { cat("Maximum corr.: ", round(des_cor_max, 3), "\n") }
    if (des_cor_max > .9) {
      warning("Maximum corr. between design matrix columns is high (",
        round(des_cor_max, 3), "). The design may ",
        "be too collinear, causing issues for model estimation. ",
        "Please visually inspect the design matrix with `plot` and modify ",
        "as necessary. ")
    }

    ### VIF --------------------------------------------------------------------
    des_vif <-  checkX_VIF(design)
    if (inherits(des_vif, "try-error")) {
      des_vif_max <- des_vif
    } else {
      des_vif_max <- max(des_vif, na.rm=TRUE)
      if (verbose) { cat("Maximum VIF:   ", round(des_vif_max, 3), "\n") }
      if (des_vif_max > 10) {
        warning("Maximum VIF is high (", round(des_vif_max, 3), "). The ",
          "design may be too collinear, causing issues for model estimation. ",
          "Please visually inspect the design matrix with `plot` and modify ",
          "as necessary. ")
      }
    }

    ### Min time btwn onset/offset ---------------------------------------------
    if (!is.null(onset) || !is.null(offset)) {
      tdiff_min <- min(diff(sort(c(all_onset_onset, all_offset_onset))))
      if (verbose) { cat("Min. time btwn onset/offset (s): ", round(tdiff_min, 3), "\n") }
      if (tdiff_min < 1) {
        warning("Min. time btwn onset/offset is small (", round(tdiff_min, 3), " s). ")
      }
    } else {
      tdiff_min <- NA
    }
    cat("\n")
  }

  # Format results and return. -------------------------------------------------

  # `FIR` stuff.
  # Binarize FIR matrix.
  # if (do_FIR) out$FIR <- 1*(out$FIR > 0)
  # FIR <- NULL # lapply(x, '[[', "FIR") #T x K x nFIR -- just return this for now, don't use in modeling
  # design_FIR <- if (is.null(FIR)) { NULL } else {
  #   lapply(FIR, function(q){
  #     # [NOTE] Untested
  #     q <- as.matrix(q, ncol=nJ)
  #     colnames(q) <- c(outer(paste0(task_names, "_FIR"), seq(nJ), paste0))
  #     q
  #   })
  # }

  colnames(design) <- field_names
  stimulus <- do.call(cbind, stimulus)

  out <- list(
    design=design,
    field_names=field_names,
    dHRF=dHRF,
    HRF_info=list(
      u=u,
      HRF=HRF,
      stimulus=stimulus
    ),
    diagnostics = list(
      cor = des_cor_max,
      vif = des_vif_max,
      tdiff = tdiff_min
    )
  )
  class(out) <- "BfMRI_design"
  out
}

#' Is this a valid entry in `EVs`?
#'
#' Is this valid data for a single task's EVs? Expects a data.frame or
#'  numeric matrix with two numeric columns, EVs and durations, and at least
#'  one row. Or, the value \code{NA} for an empty task.
#'
#' @param EV The putative EVs matrix or data frame
#' @keywords internal
#' @return Length-one logical vector.
#'
format_EV <- function(EV){
  # First check if EVs is NA, which can be the case in multi-session analysis
  #   where not all fields are present in all sessions.
  if(length(EV)==1 && is.na(EV)) { return(NA) }

  # Otherwise: For each task, expect a two-column numeric matrix; second >= 0.
  #   Except missing tasks, which should just be the value `NA`.
  is_nummat <- is.numeric(EV) && is.matrix(EV)
  is_df <- is.data.frame(EV) && all(vapply(EV, class, "") %in% c("integer", "numeric"))
  if (!(is_nummat || is_df)) { stop("The EVs are not a numeric matrix or data.frame.") }

  if (nrow(EV)<1) { stop("The EVs must have at least one row.") }

  if (ncol(EV) < 2) {
    stop("The EVs should have two columns named `onset` and `duration`. But there's less than two columns.")
  } else if (ncol(EV) > 2) {
    #message("The EVs should have two columns named `onset` and `duration`. Using the first two, and dropping the other columns.")
    EV <- EV[,seq(2)]
  }
  if (is.null(colnames(EV)) || !all(colnames(EV) == c("onset", "duration"))) {
    #message("The EVs should have two columns named `onset` and `duration`.")
    colnames(EV) <- c("onset", "duration")
  }

  as.data.frame(EV)
}
