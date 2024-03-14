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
#'  Set to \code{0} to only model the main HRF regressor, and not include its
#'  derivatives; set to \code{1} (default) to model the temporal derivative too;
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
#'  \code{BayesGLM_cifti}.
#' @param upsample Upsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param onset,offset Add task regressors indicating the onset and/or offset of
#'  each event? Provide the names of the tasks as a character vector. All
#'  onsets (or offsets) across the specified tasks will be represented by one
#'  additional column in the design matrix. The task names must match the names
#'  of \code{EVs}. Can also be \code{"all"} to use all tasks.
#' @param scale_design Scale the columns of the design matrix? Default: \code{TRUE}.
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
  EVs, nTime, TR, dHRF=c(1, 0, 2), upsample=100,
  onset=NULL, offset=NULL,
  scale_design=TRUE, ...
){

  # Check inputs. --------------------------------------------------------------
  stopifnot(is.list(EVs))
  EVs <- lapply(EVs, format_EV)
  stopifnot(!any(lapply(EVs, isFALSE)))
  stopifnot(fMRItools::is_1(nTime, "numeric"))
  stopifnot(all(nTime > 0))
  stopifnot(all(nTime == round(nTime)))
  stopifnot(fMRItools::is_1(TR, "numeric"))
  stopifnot(TR > 0)
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))
  upsample <- round(TR*upsample)/TR # TR*upsample must be an integer

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

  if (!is.null(onset)) {
    stopifnot(is.character(onset))
    if (length(onset)==1 && onset=="all") { onset <- task_names }
    stopifnot(all(onset %in% task_names))
    onset <- EVs[onset]
    onset <- lapply(onset, function(q){q$duration <- 0; q})
    onset <- do.call(rbind, onset)
    onset$onset <- sort(onset$onset)
    EVs <- c(EVs, list(onset=onset))
  }

  if (!is.null(offset)) {
    stopifnot(is.character(offset))
    if (length(offset)==1 && offset=="all") { offset <- task_names }
    stopifnot(all(offset %in% task_names))
    offset <- EVs[offset]
    offset <- lapply(offset, function(q){q$onset <- q$onset + q$duration; q$duration <- 0; q})
    offset <- do.call(rbind, offset)
    offset$onset <- sort(offset$onset)
    EVs <- c(EVs, list(offset=offset))
  }

  if ((!is.null(onset)) && (!is.null(offset))) {
    if (any(onset$onset %in% offset$onset)) {
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
  nK <- (dHRF+1)*nJ0 + (!is.null(onset)) + (!is.null(offset))

  if ((!is.null(onset)) || (!is.null(offset))) {
    if (dHRF>0) { cat("Not making fields for HRF derivatives of `onset` or ",
      "`offset`. If these are desired, provide with `EVs`.\n")
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

    # Skip if no stimulus for this task.
    if (length(EVs[[jj]])==1 && is.na(EVs[[jj]])) { #NA is placeholder for missing tasks
      warning(
        "For task '", task_names[jj],
        "': inputting constant-zero column in design matrix for lack of stimulus. ",
        "Proceed with caution."
      )
      next()
    }

    # Define `ots` (onset times) and `dur` (durations) for this task.
    ots_sj <- EVs[[jj]][,1]
    dur_sj <- EVs[[jj]][,2]

    ### Round `ots_sj` to TR for event-related designs.
    ots_sj <- ifelse(dur_sj==0, round(ots_sj/TR)*TR, ots_sj)

    ##### Get `stimulus`. -----------------------------------------------------
    stim_sj <- rep(0, nSec*upsample)
    for (ii in seq(length(ots_sj))) {
      start_ii <- round(ots_sj[ii]*upsample)
      end_ii <- round(ots_sj[ii]*upsample + dur_sj[ii]*upsample)
      stim_sj[start_ii:end_ii] <- 1
    }
    stimulus[[jj]] <- c(stim_sj,0)[inds]

    ##### Get `design` by convolving stimulus and HRFs. ------------------------
    # Convolve.
    HRF_jj <- if (is_onset_or_offset) { HRF[1] } else { HRF }
    HRF_conv <- lapply(HRF_jj, function(q){
      convolve(stim_sj, rev(q), type="open")
    })

    # Normalize each by dividing by its maximum, so the peak = 1.
    # Note that this occurs prior to downsampling.
    HRF_conv <- lapply(HRF_conv, function(q){ q / max(q) })

    # Downsample and save to `design`.
    design[,field_idx] <- do.call(cbind, HRF_conv)[inds,]

    ##### Get `FIR`. -----------------------------------------------------------
    if (FIR_nSec > 0) {
      for (tt in seq(FIR_nTimeR)) {
        dummy_tt <- rep(0, TR*upsample*FIR_nTimeR)
        start_tt <- inds_FIR[tt]
        end_tt <- start_tt + TR*upsample - 1 #each FIR basis spans a single TR
        dummy_tt[start_tt:end_tt] <- 1
        FIR_tj <- convolve(stim_sj, rev(dummy_tt), type='open')
        FIR_tj <- round(FIR_tj, 5)
        # out$FIR[,jj,tt] <- FIR_tj[inds]
      } #end loop over FIR bases
    } #end FIR basis set construction
  } #end loop over tasks

  # Diagnostics ----------------------------------------------------------------

  # Correlation
  des_cor <- cor(design)
  des_cor_max <- max(des_cor[upper.tri(des_cor)])
  if (des_cor_max < .9) {
    cat("Maximum corr.: ", round(des_cor_max, 3), "\n")
  } else {
    warning("Maximum corr. between design matrix columns is high (", 
      round(des_cor_max, 3), "). The design may ",
      "be too collinear, causing issues for model estimation. Consider ",
      "modeling some fields as nuisance instead of task, if possible.")
  }

  # VIF
  f_rhs <- paste(field_names, collapse = ' + ')
  des2 <- design
  colnames(des2) <- field_names
  des2 <- as.data.frame(des2)
  des2$y <- rnorm(nTime) #add fake y variable, has no influence
  f <- as.formula(paste0('y ~ ',f_rhs))
  des_vif <-  car::vif(lm(f, data = des2))
  des_vif_max <- max(des_vif)
  if (des_vif_max < 5) {
    cat("Maximum VIF:   ", round(des_vif_max, 3), "\n")
  } else {
    warning("Maximum VIF is high (", round(des_vif_max, 3), "). The design may ",
      "be too collinear, causing issues for model estimation. Consider ",
      "modeling some fields as nuisance instead of task, if possible.")
  }

  # onset-offset minimum time between
  if (!is.null(onset) && !is.null(offset)) {
    tdiffs <- outer(onset$onset, offset$onset, '-')
    tdiff_min <- min(abs(tdiffs))
    if (tdiff_min > 1) {
      cat("Minimum onset-offset time diff: ", round(tdiff_min, 3), "\n")
    } else {
      warning("Minimum onset-offset time diff is small (", round(tdiff_min, 3), "). ")
    }
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

  design <- if(scale_design) {
    scale_design_mat(design)
  } else {
    scale(design, scale=FALSE)
  }

  out <- list(
    design=design,
    field_names=field_names,
    dHRF=dHRF,
    HRF_info=list(
      u=u,
      HRF=HRF,
      stimulus=stimulus
    ),
    per_location_modeling=FALSE
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
  if(length(EV)==1 && is.na(EV)) { return(TRUE) }

  # Otherwise: For each task, expect a two-column numeric matrix; second >= 0.
  #   Except missing tasks, which should just be the value `NA`.
  is_nummat <- is.numeric(EV) && is.matrix(EV)
  is_df <- is.data.frame(EV) && all(vapply(EV, class, "") %in% c("integer", "numeric"))
  if (!(is_nummat || is_df)) { warning("The EVs are not a numeric matrix or data.frame."); return(FALSE) }

  if (nrow(EV)<1) { warning("The EVs must have at least one row."); return(FALSE) }

  if (ncol(EV) < 2) { 
    warning("The EVs should have two columns named `onset` and `duration`.")
    return(FALSE)
  } else if (ncol(EV) > 2) {
    message("The EVs should have two columns named `onset` and `duration`. Using the first two, and dropping the other columns.")
    EV <- EV[,seq(2)]
  }
  if (is.null(colnames(EV)) || !all(colnames(EV) == c("onset", "duration"))) {
    message("The EVs should have two columns named `onset` and `duration`.")
    colnames(EV) <- c("onset", "duration")
  }

  EV
}
