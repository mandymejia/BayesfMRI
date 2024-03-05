#' Make design matrix
#'
#' Make the design matrix for the GLM, from the task information.

#' @param onsets Task stimulus information, from which a design matrix will be 
#'  constructed. This is a list where each entry represents a task as a matrix
#'  of onsets (first column) and durations (second column) for each stimuli 
#'  (each row) of the task, in seconds. List names should be the task names. 
#'  \code{nTime} and \code{TR} are required with \code{onsets}.
#'
#'  An example of a properly-formatted \code{onsets} is:
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
#' @param downsample Downsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param ... Additional arguments to \code{\link{HRF_calc}}.
#'
#' @importFrom fMRItools is_1
#' @importFrom stats convolve
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
  onsets, nTime, TR, dHRF=c(1, 0, 2), downsample=100, ...
){

  # Check inputs. --------------------------------------------------------------
  stopifnot(is.list(onsets))
  stopifnot(all(vapply(onsets, is_onsets, FALSE)))
  stopifnot(fMRItools::is_1(nTime, "numeric"))
  stopifnot(all(nTime > 0))
  stopifnot(all(nTime == round(nTime)))
  stopifnot(fMRItools::is_1(TR, "numeric"))
  stopifnot(TR > 0)
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))
  downsample <- round(TR*downsample)/TR # TR*downsample must be an integer

  # In the future, this might be an argument.
  FIR_nSec <- 0
  do_FIR <- FIR_nSec > 0
  if (do_FIR) {
    FIR_nSec <- round(FIR_nSec/TR)*TR # FIR_nSec/TR must be an integer to match indices of downsampled stimulus, must round if not
    FIR_nTimeR <- FIR_nSec/TR # (sec)/(sec/vol) = vol  #number of TRs to model with FIR
  }

  # Define dimension variables and names.
  task_names <- names(onsets)
  nJ <- length(task_names)

  # Compute HRFs. --------------------------------------------------------------

  # Prepare to downsample the stimulus function.
  nSec <- nTime*TR; # Total time of experiment in seconds
  inds <- seq(TR*downsample, nTime*TR*downsample, by = TR*downsample) # Extract EVs by TR after convolution
  #note: since nTime is an integer, the inds will always be integers, as long as TR*downsample is also an integer
  if(max(abs(inds - round(inds))) > 1e-6) stop('Contact developer: detected non-integer indices for downsampling')
  inds <- round(inds) #to fix tiny deviations from zero

  # For FIR.
  if (do_FIR) {
    inds_FIR <- seq(1, by = TR*downsample, length.out = FIR_nTimeR) #first index for each FIR basis function
    if(max(abs(inds_FIR - round(inds_FIR))) > 1e-6) stop('Contact developer: detected non-integer indices for FIR. Set FIR_nSec = 0 to skip FIR basis set calculation.')
  }

  # Canonical HRF to use in convolution as a function of time (in seconds).
  u <- seq(0, 30, by=1/downsample) #go out 30 sec
  HRF <- list(c=HRF_calc(t = u, deriv=0, ...)[-1])
  if (dHRF > 0) { HRF$d = HRF_calc(t = u, deriv=1, ...)[-1] }
  if (dHRF > 1) { HRF$dd = HRF_calc(t = u, deriv=2, ...)[-1] }

  ### Iterate over tasks. ------------------------------------------------------
  stimulus  <- vector("list", nJ)
  design <- matrix(NA, nrow=nTime, ncol=(dHRF+1)*nJ)
  field_names <- vector("character", (dHRF+1)*nJ)
  for (jj in seq(nJ)) {
    field_idx <- jj + seq(0, dHRF)*nJ
    # Field names for this task.
    field_names[field_idx] <- paste0(
      task_names[jj],
      c("", "_dHRF", "_ddHRF")[seq(dHRF+1)]
    )

    # Skip if no stimulus for this task.
    if (length(onsets[[jj]])==1 && is.na(onsets[[jj]])) { #NA is placeholder for missing tasks
      warning(
        "For task '", task_names[jj],
        "': inputting constant-zero column in design matrix for lack of stimulus. ",
        "Proceed with caution."
      )
      next()
    }

    # Define `ots` (onset times) and `dur` (durations) for this task.
    ots_sj <- onsets[[jj]][,1]
    dur_sj <- onsets[[jj]][,2]

    ### Round `ots_sj` to TR for event-related designs.
    ots_sj <- ifelse(dur_sj==0, round(ots_sj/TR)*TR, ots_sj)

    ##### Get `stimulus`. -----------------------------------------------------
    stim_sj <- rep(0, nSec*downsample)
    for (ii in seq(length(ots_sj))) {
      start_ii <- round(ots_sj[ii]*downsample)
      end_ii <- round(ots_sj[ii]*downsample + dur_sj[ii]*downsample)
      stim_sj[start_ii:end_ii] <- 1
    }
    stimulus[[jj]] <- c(stim_sj,0)[inds]

    ##### Get `design` by convolving stimulus and HRFs. ------------------------
    design[,field_idx] <- do.call(cbind,
      lapply(HRF, function(q){ convolve(stim_sj, rev(q), type="open") })
    )[inds,]

    ##### Get `FIR`. -----------------------------------------------------------
    if (FIR_nSec > 0) {
      for (tt in seq(FIR_nTimeR)) {
        dummy_tt <- rep(0, TR*downsample*FIR_nTimeR)
        start_tt <- inds_FIR[tt]
        end_tt <- start_tt + TR*downsample - 1 #each FIR basis spans a single TR
        dummy_tt[start_tt:end_tt] <- 1
        FIR_tj <- convolve(stim_sj, rev(dummy_tt), type='open')
        FIR_tj <- round(FIR_tj, 5)
        # out$FIR[,jj,tt] <- FIR_tj[inds]
      } #end loop over FIR bases
    } #end FIR basis set construction
  } #end loop over tasks

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
    per_location_modeling=FALSE
  )
  class(out) <- "BfMRI_design"
  out
}

#' Is this a valid entry in `onsets`?
#'
#' Is this valid data for a single task's onsets? Expects a data.frame or
#'  numeric matrix with two numeric columns, onsets and durations, and at least
#'  one row. Or, the value \code{NA} for an empty task.
#'
#' @param x The putative onsets matrix or data frame
#' @keywords internal
#' @return Length-one logical vector.
#'
is_onsets <- function(x){
  # First check if onsets is NA, which can be the case in multi-session analysis
  #   where not all fields are present in all sessions.
  if(length(x)==1 && is.na(x)) { return(TRUE) }

  # Otherwise: For each task, expect a two-column numeric matrix; second >= 0.
  #   Except missing tasks, which should just be the value `NA`.
  is_nummat <- is.numeric(x) & is.matrix(x)
  is_df <- is.data.frame(x) && all(vapply(x, class, "") == "numeric")
  if (!(is_nummat || is_df)) { warning("The onsets are not a numeric matrix or data.frame."); return(FALSE) }

  if (nrow(x)<1) { warning("The onsets must have at least one row."); return(FALSE) }

  if (ncol(x) != 2) { warning("The onsets should have two columns named `onset` and `duration`."); return(FALSE) }
  if (!all(colnames(x) == c("onset", "duration"))) {
    warning("The onsets should have two columns named `onset` and `duration`."); return(FALSE)
  }

  TRUE
}
