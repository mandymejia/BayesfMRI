#' Make HRFs
#'
#' Based on onsets and durations, compute canonical HRF, its derivatives, and FIR basis set for each task
#'
#' @param onsets \eqn{L}-length list in which the name of each element is the
#'  name of the corresponding task, and the value of each element is a matrix of
#'  onsets (first column) and durations (second column) for each stimuli (each
#'  row) of the corresponding task.
#' @param TR Temporal resolution of the data, in seconds.
#' @param nT The number of volumes in the fMRI data.
#' @param FIR_nsec The number of seconds to cover in the FIR basis set. Default: \code{0} (do not return FIR bases)
#' @param upsample Upsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param dHRF Number of HRF derivatives to include in the model. If \code{dHRF=1}
#' include the temporal derivative of the HRF, if \code{dHRF=2} also include the
#' dispersion derivative. Default: \code{2}.
#' @param a1 (Optional) delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{1}
#' @param a2 delay of undershoot. Default: \code{16/6 * a1 * sqrt(b1) = 16}
#' @param b2 dispersion of undershoot. Default: \code{b1 = 1}
#' @param c scale of undershoot. Default: \code{1/6}
#'
#' @return List with the HRF and derivatives, the stimulus stick function, and the FIR basis set.
#'
#' @importFrom stats convolve
#'
#' @examples
#' onsets <- list(taskA=cbind(c(2,17,23),4)) # one task, 3 four sec-long stimuli
#' TR <- .72 # .72 seconds per volume, or (1/.72) Hz
#' nT <- 300 # session is 300 volumes long (300*.72 seconds long)
#' make_HRFs(onsets, TR, nT)
#'
#' @export
make_HRFs <- function(
  onsets,
  TR,
  nT,
  FIR_nsec=0,
  upsample=100,
  dHRF=2,
  a1=6,
  b1=1,
  a2=16/6 * a1 * sqrt(b1),
  b2=1,
  c=1/6){

  nK <- length(onsets) #number of tasks
  task_names <- if (is.null(names(onsets))) {
    task_names <- paste0('task', 1:nK)
  } else {
    task_names <- names(onsets)
  }

  #prepare to upsample the stimulus function
  upsample <- round(TR*upsample)/TR # TR*upsample must be an integer
  nsec <- nT*TR; # Total time of experiment in seconds
  stimulus <- rep(0, nsec*upsample) # For stick function to be used in convolution
  inds <- seq(TR*upsample, nT*TR*upsample, by = TR*upsample) # Extract EVs by TR after convolution
  #note: since nT is an integer, the inds will always be integers, as long as TR*upsample is also an integer

  if(max(abs(inds - round(inds))) > 1e-6) stop('Contact developer: detected non-integer indices for downsampling')
  inds <- round(inds)

  #for FIR
  if(FIR_nsec > 0){
    FIR_nsec <- round(FIR_nsec/TR)*TR # FIR_nsec/TR must be an integer to match indices of upsampled stimulus, must round if not
    FIR_nTR <- FIR_nsec/TR # (sec)/(sec/vol) = vol  #number of TRs to model with FIR
    inds_FIR <- seq(1, by = TR*upsample, length.out = FIR_nTR) #first index for each FIR basis function
    if(max(abs(inds_FIR - round(inds_FIR))) > 1e-6) stop('Contact developer: detected non-integer indices for FIR. Set FIR_nsec = 0 to skip FIR basis set calculation.')
  }

  #initialize arrays
  theStimulus <- array(NA, dim=c(nT, nK),
                       dimnames = list(volume=1:nT, task=task_names))
  theHRFs <- array(NA, dim=c(nT, nK, 3),
                   dimnames = list(volume=1:nT, task=task_names, HRF=c("HRF", "dHRF", "ddHRF")))
  if(FIR_nsec > 0) theFIR <- array(NA, dim=c(nT, nK, FIR_nTR),
                  dimnames = list(volume=1:nT, task=task_names, FIR=1:FIR_nTR))
  if(!(FIR_nsec > 0)) theFIR <- NULL

  # canonical HRF to be used in convolution as a function of time (in seconds)
  u <- seq(0, 30, by=1/upsample) #go out 30 sec
  HRFs <- list(HRF = HRF(t = u, deriv=0, a1, b1, a2, b2, c),
               dHRF = HRF(t = u, deriv=1, a1, b1, a2, b2, c),
               ddHRF = HRF(t = u, deriv=2, a1, b1, a2, b2, c))

  for (kk in seq(nK)) {

    onsets_k <- onsets[[kk]][,1] #onset times in scans
    durations_k <- onsets[[kk]][,2] #durations in scans

    # Round onsets to TR for event-related designs
    round2TR <- function(x, TR){ round(x/TR)*TR }
    onsets_k <- ifelse(durations_k == 0, round2TR(onsets_k, TR), onsets_k)

    # Define stimulus function
    stimulus_k <- stimulus
    for (ii in seq(length(onsets_k))) {
      start_ii <- round(onsets_k[ii]*upsample)
      end_ii <- round(onsets_k[ii]*upsample + durations_k[ii]*upsample)
      stimulus_k[start_ii:end_ii] <- 1
    }
    theStimulus[,kk] <- c(stimulus_k,0)[inds]

    ### Canonical HRF and Derivatives

    #1 = HRF only, 2 = HRF + dHRF, 3 = HRF + dHRF + ddHRF
    for (dd in 1:(dHRF+1)) {
      HRF_k <- convolve(stimulus_k, rev(HRFs[[dd]]), type='open')
      theHRFs[,kk,dd] <- HRF_k[inds]
    }

    ### FIR Basis Set

    if(FIR_nsec > 0){
      for(tt in 1:FIR_nTR){
        dummy_tt <- rep(0, TR*upsample*FIR_nTR)
        start_tt <- inds_FIR[tt]
        end_tt <- start_tt + TR*upsample - 1 #each FIR basis spans a single TR
        dummy_tt[start_tt:end_tt] <- 1
        FIR_tk <- convolve(stimulus_k, rev(dummy_tt), type='open')
        FIR_tk <- round(FIR_tk, 5)
        theFIR[,kk,tt] <- FIR_tk[inds]
        } #end loop over FIR bases
      } #end FIR basis set construction

    } #end loop over tasks

  ### Format results

  #starting point for design matrix (may add derivatives later)
  design <- matrix(theHRFs[,,1], ncol=nK) #reform into a matrix in case nK=1
  colnames(design) <- task_names

  #binarize FIR matrix
  if(!is.null(theFIR)) theFIR <- 1*(theFIR > 0)

  #add time to HRFs
  HRFs_df <- data.frame(time = u, HRF = HRFs[[1]], dHRF = HRFs[[2]], ddHRF = HRFs[[3]])

  list(
    stimulus=theStimulus,
    HRF_convolved=theHRFs,
    HRF=HRFs_df, #HRFs and derivatives to use in convolution
    design=design, #starting point for design matrix
    FIR=theFIR
  )
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
#'
#' @return HRF vector (or dHRF, or d2HRF) corresponding to time vector t
#'
#' @importFrom fMRItools is_1
#' @export
#'
HRF <- function(t, deriv=0, a1=6, b1=1, a2=NULL, b2=NULL, c=1/6){

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

  #Temporal derivative [TO DO] emulate SPM to make this truly an onset derivative
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
#' @param t time vector (in seconds)
# @param TR temporal resolution of the data, in seconds
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{1}
#' @param a2 delay of undershoot. Default: \code{16/6*a1 = 16}
#' @param b2 dispersion of undershoot. Default: \code{b1 = 1}
#' @param c scale of undershoot. Default: \code{1/6}
# @param mt microtime resolution (number of bins to divide each TR into). Default: \code{16}
#'
#' @return HRF vector corresponding to time vector t
#'
#' @importFrom fMRItools is_1
#' @importFrom stats dgamma
#' @export
#'
HRF_main <- function(t, a1=6, b1=1, a2=NULL, b2=NULL, c=1/6){ #}, mt=16){

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
