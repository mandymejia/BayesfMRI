#' Make HRFs
#' 
#' Make HRFs from onsets. This is a helper function to \code{\link{make_design}}.
#'
#' @param onsets,nTime,TR,dHRF See \code{\link{make_design}}.
# @param FIR_nSec The number of seconds to cover in the FIR basis set. Default: \code{0} (do not return FIR bases)
#' @param upsample Upsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param a1 (Optional) delay of response. Default: \code{6}
#' @param b1 Response dispersion. Default: \code{1}
#' @param a2 Delay of undershoot. Default: \code{16/6 * a1 * sqrt(b1) == 16}
#' @param b2 Dispersion of undershoot. Default: \code{b1 == 1}
#' @param c Scale of undershoot. Default: \code{1/6}
#' @return A session-length list. Each entry is a sublist with components:
#'  \describe{
#'    \item{HRF}{The upsampled HRFs.}
#'    \item{stimulus}{The task stimulus.}
#'    \item{design}{The corresponding design matrix.}
#' }
#' @importFrom stats convolve
#'
#' @examples
#' # onsets: 1 session, 1 task, 3 stimuli, each 4 sec. long
#' onsets <- list(session1=list(taskA=cbind(c(2,17,23),4))) 
#' TR <- .72 # .72 seconds per volume, or (1/.72) Hz
#' nTime <- 300 # session has 300 volumes. Session duration is 300*.72 sec.
#' make_HRFs(onsets, nTime, TR)
#'
#' @export
make_HRFs <- function(
  onsets, nTime, TR, dHRF=c(1, 0, 2),
  upsample=100,
  a1=6,
  b1=1,
  a2=16/6 * a1 * sqrt(b1),
  b2=b1,
  c=1/6){

  # Check input arguments. -----------------------------------------------------
  # `onsets`.
  ### For each task, expect a two-column numeric matrix; second column >= 0.
  ###   Except missing tasks, which should just be the value `NA`.
  stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(r){
    if (length(r)==1 && is.na(r)) {TRUE} else { all(apply(r, 2, is.numeric)) }
  }, FALSE))}, FALSE)))
  stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(r){
    if (length(r)==1 && is.na(r)) {TRUE} else { is.matrix.or.df(r) }
  }, FALSE))}, FALSE)))
  stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(r){
    if (length(r)==1 && is.na(r)) {TRUE} else { ncol(r)==2 }
  }, FALSE))}, FALSE)))
  stopifnot(all(vapply(onsets, function(q){all(vapply(q, function(r){
    if (length(r)==1 && is.na(r)) {TRUE} else { all(r[,2] >= 0) }
  }, FALSE))}, FALSE)))
  ### Expect every session to have the same tasks.
  stopifnot(length(unique(lapply(onsets, length)))==1)
  stopifnot(length(unique(lapply(onsets, names)))==1)
  nS <- length(onsets)

  # `nTime`.
  if (is.null(nTime)) { stop("`nTime` is required if `onsets` is provided.") }
  stopifnot(is.numeric(nTime))
  stopifnot(all(nTime > 0))
  stopifnot(all(nTime == round(nTime)))
  if (length(nTime)==1) {
    nTime <- rep(nTime, nS)
  } else {
    stopifnot(length(nTime)==nS)
  }

  # `TR`.
  if (is.null(TR)) { stop("`TR` is required if `onsets` is provided.") }
  stopifnot(is.numeric(TR))
  stopifnot(TR > 0)

  # `dHRF`.
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("1", "0", "2")))

  # Others.
  stopifnot(is_1(upsample, "numeric"))
  upsample <- round(TR*upsample)/TR # TR*upsample must be an integer
  stopifnot(is_1(a1, "numeric"))
  stopifnot(is_1(b1, "numeric"))
  stopifnot(is_1(a2, "numeric"))
  stopifnot(is_1(b2, "numeric"))
  stopifnot(is_1(c, "numeric"))

  # In the future, this might be an argument.
  FIR_nSec <- 0
  do_FIR <- FIR_nSec > 0
  if (do_FIR) {
    FIR_nSec <- round(FIR_nSec/TR)*TR # FIR_nSec/TR must be an integer to match indices of upsampled stimulus, must round if not
    FIR_nTimeR <- FIR_nSec/TR # (sec)/(sec/vol) = vol  #number of TRs to model with FIR
  }

  # Define dimension variables and names.
  session_names <- names(onsets)
  nS <- length(onsets)
  task_names <- names(onsets[[1]])
  nJ <- length(task_names)
  HRF_names <- c("HRF", "dHRF", "ddHRF")[seq(dHRF+1)]
  field_names <- c(outer(
    task_names,
    c("", "_dHRF", "_ddHRF")[seq(dHRF+1)],
    paste0
  ))

  # Initialize return variables.
  out <- setNames(vector("list", nS), session_names)
  for (ss in seq(nS)) {
    out[[ss]] <- list(
      HRF = NULL,
      stimulus = array(0, dim=c(nTime[ss], nJ),
        dimnames=list(vol=NULL, task=task_names)
      ),
      design = array(0, dim=c(nTime[ss], nJ*(dHRF+1)),
        dimnames=list(vol=NULL, field=field_names)
      )#,
      # FIR = if (do_FIR) {
      #   array(0, dim=c(nTime, nJ, FIR_nTimeR),
      #     dimnames = list(vol=NULL, task=task_names, FIR=seq(FIR_nTimeR))
      #   )
      # } else {
      #   NULL
      # }
    )
  }

  # Iterate over sessions. -----------------------------------------------------
  for (ss in seq(nS)) {
    ### Compute `HRF`. ---------------------------------------------------------
    # For multi-session, avoid re-computing if `nT` is unchanged.
    if (ss==1 || (nT_ss != nTime[ss])) {
      nT_ss <- nTime[ss]

      # Prepare to upsample the stimulus function.
      nSec <- nT_ss*TR; # Total time of experiment in seconds
      inds <- seq(TR*upsample, nT_ss*TR*upsample, by = TR*upsample) # Extract EVs by TR after convolution
      #note: since nT_ss is an integer, the inds will always be integers, as long as TR*upsample is also an integer
      if(max(abs(inds - round(inds))) > 1e-6) stop('Contact developer: detected non-integer indices for downsampling')
      inds <- round(inds) #to fix tiny deviations from zero

      # For FIR.
      if (do_FIR) {
        inds_FIR <- seq(1, by = TR*upsample, length.out = FIR_nTimeR) #first index for each FIR basis function
        if(max(abs(inds_FIR - round(inds_FIR))) > 1e-6) stop('Contact developer: detected non-integer indices for FIR. Set FIR_nSec = 0 to skip FIR basis set calculation.')
      }

      # Canonical HRF to use in convolution as a function of time (in seconds).
      u <- seq(0, 30, by=1/upsample) #go out 30 sec
      out[[ss]]$HRF <- list(
        u = u,
        HRF = HRF(t = u, deriv=0, a1, b1, a2, b2, c)
      )
      if (dHRF > 0) {
        out[[ss]]$HRF <- c(out[[ss]]$HRF, list(
          HRF_d = HRF(t = u, deriv=1, a1, b1, a2, b2, c)
        ))
      }
      if (dHRF > 1) {
        out[[ss]]$HRF <- c(out[[ss]]$HRF, list(
          HRF_dd = HRF(t = u, deriv=2, a1, b1, a2, b2, c)
        ))
      }
    } else {
      out[[ss]]$HRF <- out[[ss-1]]$HRF
    }

    ### Iterate over tasks. ----------------------------------------------------
    for (jj in seq(nJ)) {

      # Skip if no stimulus for this session & task.
      if (all(is.na(onsets[[ss]][[jj]]))) { #NA is placeholder for missing tasks
        warning(
          "For task '", task_names[jj], "' in session '", session_names[ss],
          "': inputting constant-zero column in design matrix for lack of stimulus. ",
          "Proceed with caution."
        )
        next()
      }

      # Define `ots` (onset times) and `dur` (durations) for this sess & task.
      ots_sj <- onsets[[ss]][[jj]][,1]
      dur_sj <- onsets[[ss]][[jj]][,2]

      ### Round `ots_sj` to TR for event-related designs.
      ots_sj <- ifelse(dur_sj==0,
        round(ots_sj/TR)*TR,
        ots_sj
      )

      ##### Get `stimulus`. ---------------------------------------------------
      stim_sj <- rep(0, nSec*upsample)
      for (ii in seq(length(ots_sj))) {
        start_ii <- round(ots_sj[ii]*upsample)
        end_ii <- round(ots_sj[ii]*upsample + dur_sj[ii]*upsample)
        stim_sj[start_ii:end_ii] <- 1
      }
      out[[ss]]$stimulus[,jj] <- c(stim_sj,0)[inds]

      ##### Get `design` by convolving stimulus and HRFs. ----------------------
      #1 = HRF only, 2 = HRF + dHRF, 3 = HRF + dHRF + ddHRF
      for (dd in seq(0, dHRF)) {
        HRF_dd <- out[[ss]]$HRF[[c("HRF", "HRF_d", "HRF_dd")[dd+1]]]
        out[[ss]]$design[,jj+dd*nJ] <- convolve(
          stim_sj, rev(HRF_dd), type='open'
        )[inds]
      }

      ##### Get `FIR`. ---------------------------------------------------------
      if (FIR_nSec > 0) {
        for (tt in seq(FIR_nTimeR)) {
          dummy_tt <- rep(0, TR*upsample*FIR_nTimeR)
          start_tt <- inds_FIR[tt]
          end_tt <- start_tt + TR*upsample - 1 #each FIR basis spans a single TR
          dummy_tt[start_tt:end_tt] <- 1
          FIR_tj <- convolve(stim_sj, rev(dummy_tt), type='open')
          FIR_tj <- round(FIR_tj, 5)
          # out[[ss]]$FIR[,jj,tt] <- FIR_tj[inds]
        } #end loop over FIR bases
      } #end FIR basis set construction
    } #end loop over tasks
  } # end loop over sessions

  # Format results and return. -------------------------------------------------

  # `FIR` stuff.
  # Binarize FIR matrix.
  # if (do_FIR) out[[ss]]$FIR <- 1*(out[[ss]]$FIR > 0)
  # FIR <- NULL # lapply(x, '[[', "FIR") #T x K x nFIR -- just return this for now, don't use in modeling
  # design_FIR <- if (is.null(FIR)) { NULL } else {
  #   lapply(FIR, function(q){
  #     # [NOTE] Untested
  #     q <- as.matrix(q, ncol=nJ)
  #     colnames(q) <- c(outer(paste0(task_names, "_FIR"), seq(nJ), paste0))
  #     q
  #   })
  # }

  out
}

#' Canonical HRF and Derivatives
#'
#' Calculate the HRF from a time vector and parameters, or its derivative with
#' respect to delay or dispersion.
#'
#' @param t time vector (in units of seconds)
# @param TR temporal resolution of the data, in seconds
#' @param deriv \code{0} (default) for the HRF, \code{1} for the delay derivative
#'  of the HRF, or \code{2} for the dispersion derivative of the HRF.
# @param dt If TRUE, take the first and second derivatives of the HRF with
# respect to t, rather than a1 and b1. Default: \code{FALSE}
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{1}
#' @param a2 delay of undershoot. Default: \code{16/6 * a1 * sqrt(b1) = 16}
#' @param b2 dispersion of undershoot. Default: \code{b1 = 1}
#' @param c scale of undershoot. Default: \code{1/6}
#' @param o onset of response. Default: \code{0}
#'
#' @return HRF vector (or dHRF, or d2HRF) corresponding to time vector t
#'
#' @importFrom fMRItools is_1
#' @export
#'
HRF <- function(t, deriv=0, a1=6, b1=1, a2=NULL, b2=NULL, c=1/6, o=0){

  # Arg checks
  stopifnot(is.numeric(t)); stopifnot(min(t) >= 0)
  deriv <- as.numeric(match.arg(as.character(deriv), c("0", "1", "2")))
  stopifnot(is_1(a1, "numeric")); stopifnot(a1 > 0)
  stopifnot(is_1(b1, "numeric")); stopifnot(b1 > 0)
  stopifnot(is_1(c, "numeric")); stopifnot(c >= 0)

  #set a2 and b2 based on a1 and b1
  if(is.null(a2)) a2 <- (16/6)*a1*sqrt(b1)
  if(is.null(b2)) b2 <- b1

  #check that the shape parameters are not less than 1
  if(a1/b1 < 1) stop('The shape parameter of the first Gamma in the HRF (a1/b1) is less than 1, which is invalid. Adjust parameters.')
  if(a2/b2 < 1) stop('The shape parameter of the second Gamma in the HRF (a2/b2) is less than 1, which is invalid. Adjust parameters.')

  if(deriv==0){
    h <- HRF_main(t, a1, b1, a2, b2, c)
  }

  #Temporal derivative
  if(deriv==1){
    delta <- 0.5
    fplus <- HRF_main(t+delta, a1, b1, a2, b2, c)
    fminus <- HRF_main(t-delta, a1, b1, a2, b2, c)
    h <- (fplus - fminus) / (2*delta)
    # if(deriv==2)  h <- (fplus + fminus - 2*HRF_main(t, TR, a1, b1, a2, b1, c)) / (delta^2) #second temporal derivative
  }

  #Dispersion derivative
  if(deriv==2){
    delta <- 0.01*b1 #delta = 1% each direction
    fplus <- HRF_main(t, a1, b1+delta, a2, b2, c)
    fminus <- HRF_main(t, a1, b1-delta, a2, b2, c)
    h <- (fplus - fminus)/(2*delta)
  }

  return(h)
}

#' Canonical (double-gamma) HRF
#'
#' Calculate the HRF from a time vector and parameters. Optionally compute the
#'  first or second derivative of the HRF instead. Form of HRF is similar to SPM
#'  but here the response and undershoot are scaled so the difference of the HRFs
#'  peaks at 1 and -c
#'
#'
#' @param t time vector (in seconds). Must be equally spaced.
# @param TR temporal resolution of the data, in seconds
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{1}
#' @param a2 delay of undershoot. Default: \code{16/6*a1 = 16}
#' @param b2 dispersion of undershoot. Default: \code{b1 = 1}
#' @param c scale of undershoot. Default: \code{1/6}
#' @param o onset of response (in seconds). Default: \code{0}
# @param mt microtime resolution (number of bins to divide each TR into). Default: \code{16}
#'
#' @return HRF vector corresponding to time vector t
#'
#' @importFrom fMRItools is_1
#' @importFrom stats dgamma
#' @export
#'
HRF_main <- function(t, a1=6, b1=1, a2=NULL, b2=NULL, c=1/6, o=0){ #}, mt=16){

  #[TO DO] Check that t is equally spaced
  #[TO DO] If onsets are provided, check that it is a factor of the spacing of t, i.e. o/(t[2]-t[1]) is an integer

  TR <- 1 # this function was originally designed to return the HRF as a function of the TR/volume.
  #this is a workaround to return as a function of time

  # Arg checks
  stopifnot(is.numeric(t))
  stopifnot(is.numeric(TR))
  stopifnot(is_1(a1, "numeric")); stopifnot(a1 > 0)
  stopifnot(is_1(b1, "numeric")); stopifnot(b1 > 0)
  stopifnot(is_1(c, "numeric")); stopifnot(c >= 0)

  #set a2 and b2 based on a1 and b1
  if(is.null(a2)) a2 <- (16/6)*a1*sqrt(b1)
  if(is.null(b2)) b2 <- b1

  stopifnot(is_1(a2, "numeric")); stopifnot(a2 > 0)
  stopifnot(is_1(b2, "numeric")); stopifnot(b2 > 0)

  #check that the shape parameters are not less than 1
  if(a1/b1 < 1) stop('The shape parameter of the first Gamma in the HRF (a1/b1) is less than 1, which is invalid. Adjust parameters.')
  if(a2/b2 < 1) stop('The shape parameter of the second Gamma in the HRF (a2/b2) is less than 1, which is invalid. Adjust parameters.')

  #gamma functions
  shape1 <- a1/b1; rate1 <- TR/b1
  shape2 <- a2/b2; rate2 <- TR/b2
  rm(a1, a2, b1, b2)
  gamma1 <- dgamma(t, shape=shape1, rate = rate1)  #response
  gamma2 <- dgamma(t, shape=shape2, rate = rate2)  #undershoot
  gamma1_max <- dgamma((shape1 - 1)/rate1, shape=shape1, rate = rate1)
  gamma2_max <- dgamma((shape2 - 1)/rate2, shape=shape2, rate = rate2)

  #identify modes of gamma1 - c*gamma2 (may not be equal to mode of each individual gamma)
  times <- seq(0, 100/TR, 0.001) #go out 100 sec to ensure capturing second mode
  diff <- dgamma(times, shape=shape1, rate = rate1)/gamma1_max -
    c*dgamma(times, shape=shape2, rate = rate2)/gamma2_max
  mode1 <- times[which.max(diff)] #identify peak of response
  mode2 <- times[which.min(diff)] #identify peak of undershoot (only used when c>0)

  #solve for scaling factors of response (c1) and undershoot (c2)
  if(c > 0){
    c1_num <- c*dgamma(mode1, shape=shape2, rate = rate2) + dgamma(mode2, shape=shape2, rate = rate2)
    c1_den <- dgamma(mode1, shape=shape1, rate = rate1)*dgamma(mode2, shape=shape2, rate = rate2) -
      dgamma(mode2, shape=shape1, rate = rate1)*dgamma(mode1, shape=shape2, rate = rate2)
    c1 <- c1_num/c1_den
    c2 <- (c1*dgamma(mode1, shape=shape1, rate=rate1) - 1) / dgamma(mode1, shape=shape2, rate=rate2)
  } else {
    c1 <- 1/dgamma(mode1, shape=shape1, rate = rate1)
  }

  #hrf
  gamma1 <- c1*gamma1 #scale to have height = 1
  if(c > 0) gamma2 <- c2*gamma2 else gamma2 <- 0 #scale to have height = c
  h <- (gamma1 - gamma2)

  #onset > 0?
  if(o > 0){
    len_h <- length(h)
    shift <- round(1 + o/(t[2]-t[1])) #onset, in terms of the spacing of t
    h_new <- rep(0, len_h) #start the HRF with zeros, prior to the onset
    len_new <- length(shift:len_h)
    h_new[shift:len_h] <- h[1:len_new]
    h <- h_new
  }

  return(h)

  # #dhrf w.r.t. a1 (delay) -- this assumes that a2=16/6*a1 and b2=b1
  # if(deriv==1){
  #   part1 <- c1 * (t*b1)^(a1) * log(a1) * exp(-b1*t)/t
  #   part2 <- c2 * (t*b1)^(a1*16/6) * log(a1*16/6) * (16/6) * exp(-b1*t)/t
  #   return(part1 - part2)
  # }
}

# #dgamma without the gamma function scaling factor
# dgamma0 <- function(x, shape, rate){
#   if(shape <= 0) stop('for dgamma, shape must be positive')
#   d <- x^(shape-1) * rate^shape * exp(-rate*x)
#   d[x<0] <- 0
#   d
# }

# --------------------------------------------------------------------------

#' Canonical (double-gamma) HRF (old one from SPM96, Glover)
#'
#' Calculate the HRF from a time vector and parameters. Optionally compute the
#'  first or second derivative of the HRF instead.
#'
#' @param t time vector
#' @param deriv \code{0} (default) for the HRF, \code{1} for the first derivative
#'  of the HRF, or \code{2} for the second derivative of the HRF.
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{0.9}
#' @param a2 delay of undershoot. Default: \code{12}
#' @param b2 dispersion of undershoot. Default: \code{0.9}
#' @param c scale of undershoot. Default: \code{0.35}
#'
#' @return HRF vector (or dHRF, or d2HRF) corresponding to time
#'
#' @examples
#' upsample <- 100
#' HRF96(seq(0, 30, by=1/upsample))
#'
#' @importFrom fMRItools is_1
#' @export
#'
HRF96 <- function(t, deriv=0, a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {

  # Arg checks
  stopifnot(is.numeric(t))
  deriv <- as.numeric(match.arg(as.character(deriv), c("0", "1", "2")))
  stopifnot(is_1(a1, "numeric"))
  stopifnot(is_1(b1, "numeric"))
  stopifnot(is_1(a2, "numeric"))
  stopifnot(is_1(b2, "numeric"))
  stopifnot(is_1(c, "numeric"))

  # HRF
  if (deriv==0) {
    out <- ((t/(a1*b1))^a1) * exp(-(t-a1*b1)/b1) - c * ((t/(a2*b2))^a2) * exp(-(t - a2*b2)/b2)
    #out2 <- (t^(a1-1)*b1^a1*exp(-b1*t))/gamma(a1) - c * (t^(a2-1) * b2^a2 * exp(-b2*t))/gamma(a2)

    # dHRF
  } else if (deriv==1) {
    C1 <- (1/(a1*b1))^a1
    C2 <- c*(1/(a2*b2))^a2
    A1 <- a1*t^(a1 - 1)*exp(-(t - a1*b1)/b1)
    A2 <- a2*t^(a2 - 1)*exp(-(t - a2*b2)/b2)
    B1 <- t^a1 / b1 * exp(-(t - a1*b1)/b1)
    B2 <- t^a2 / b2 * exp(-(t - a2*b2)/b2)
    out <- C1*(A1 - B1) - C2 * (A2 - B2)

    # ddHRF
  } else if (deriv==2) {
    C1 <- (1/(a1*b1))^a1
    C2 <- c*(1/(a2*b2))^a2
    dA1 <- a1 * ((a1 - 1) - t/b1) * t^(a1-2) * exp(-(t - a1*b1)/b1)
    dB1 <- (1/b1) * (a1 - (t / b1)) * t^(a1 - 1) * exp(-(t - a1*b1)/b1)
    dA2 <- a2 * ((a2 - 1) - t/b2) * t^(a2-2) * exp(-(t - a2*b2)/b2)
    dB2 <- (1/b2) * (a2 - (t / b2)) * t^(a2 - 1) * exp(-(t - a2*b2)/b2)
    out <- C1 * (dA1 - dB1) - C2 * (dA2 - dB2)
  }

  out
}

#' Central derivative
#'
#' Take the central derivative of numeric vectors by averaging the forward and
#'  backward differences.
#' @param x A numeric matrix, or a vector which will be converted to a
#'  single-column matrix.
#' @return A matrix or vector the same dimensions as \code{x}, with the
#'  derivative taken for each column of \code{x}. The first and last rows may
#'  need to be deleted, depending on the application.
#' @export
#'
#' @examples
#' x <- cderiv(seq(5))
#' stopifnot(all(x == c(.5, 1, 1, 1, .5)))
#'
cderiv <- function(x){
  x <- as.matrix(x)
  dx <- diff(x)
  (rbind(0, dx) + rbind(dx, 0)) / 2
}
