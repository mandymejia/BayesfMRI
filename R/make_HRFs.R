#' Make HRFs
#'
#' Create HRF design matrix columns from onsets and durations
#'
#' @param onsets \eqn{K}-length list in which the name of each element is
#'  the name of the corresponding task, and the value of each element is a
#'  matrix of onsets (first column) and durations (second column) for each
#'  stimuli (each row) of the corresponding task, in SECONDS. 
#' @param TR Temporal resolution of fMRI data, in SECONDS.
#' @param duration Length of fMRI timeseries, in SCANS.
#' @param downsample Downsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param deriv \code{0} (default), \code{1}, or \code{2}, to use the HRF
#'   function, the first derivative of the HRF, or the second derivative of the
#'   HRF, respectively.
#'
#' @return Design matrix containing one HRF column for each task
#'
#' @importFrom stats convolve
#' 
#' @examples 
#' onsets <- list(taskA=cbind(c(2,17,23),4)) # one task, 3 four sec-long stimuli
#' TR <- .72 # .72 seconds per volume, or (1/.72) Hz
#' duration <- 300 # session is 300 volumes long (300*.72 seconds long)
#' make_HRFs(onsets, TR, duration)
#'
#' @export
make_HRFs <- function(onsets, TR, duration, downsample=100, deriv = 0){
  if(!deriv %in% c(0,1,2)) stop('Argument "deriv" must take the value 0, 1, or 2.')

  K <- length(onsets) #number of tasks
  if(is.null(names(onsets))) task_names <- paste0('task', 1:K) else task_names <- names(onsets)

  nsec <- duration*TR; # Total time of experiment in seconds
  stimulus <- rep(0, nsec*downsample) # For stick function to be used in convolution
  inds <- seq(TR*downsample, nsec*downsample, by = TR*downsample) # Extract EVs in a function of TR

  design <- vector("list", length=deriv+1)

  for (dd in seq(deriv+1)) {
    dname_dd <- switch(dd, "HRF", "dHRF", "d2HRF")
    design[[dd]] <- matrix(NA, nrow=duration, ncol=K)
    colnames(design[[dd]]) <- paste0(task_names, "_", dname_dd)

    HRF_fn <- switch(dd, HRF, dHRF, d2HRF)
    HRF_dd <- HRF_fn(seq(0, 30, by=1/downsample))[-1] #canonical HRF to be used in convolution

    for(k in 1:K){
      onsets_k <- onsets[[k]][,1] #onset times in scans
      durations_k <- onsets[[k]][,2] #durations in scans

      # Define stimulus function
      stimulus_k <- stimulus
      for(ii in 1:length(onsets_k)){
        start_ii <- round(onsets_k[ii]*downsample)
        end_ii <- round(onsets_k[ii]*downsample + durations_k[ii]*downsample)
        stimulus_k[start_ii:end_ii] <- 1
      }

      # Convolve boxcar with canonical HRF & add to design[[ii]] matrix
      HRF_k <- convolve(stimulus_k, rev(HRF_dd), type='open')
      design[[dd]][,k] <- HRF_k[inds]
    }
  }
  design <- do.call(cbind, design)

  return(design)
}

#' Calculate the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{0.9}
#' @param a2 delay of undershoot. Default: \code{12}
#' @param b2 dispersion of undershoot. Default: \code{0.9}
#' @param c scale of undershoot. Default: \code{0.35}
#'
#' @return HRF vector corresponding to time
#' 
#' @examples
#' downsample <- 100
#' HRF(seq(0, 30, by=1/downsample))
#' 
#' @export
HRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  ((t/(a1*b1))^a1) * exp(-(t-a1*b1)/b1) - c * ((t/(a2*b2))^a2) * exp(-(t - a2*b2)/b2)
}

#' Calculate first derivative of the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{0.9}
#' @param a2 delay of undershoot. Default: \code{12}
#' @param b2 dispersion of undershoot. Default: \code{0.9}
#' @param c scale of undershoot. Default: \code{0.35}
#'
#' @return dHRF vector corresponding to time
#' 
#' @examples
#' downsample <- 100
#' dHRF(seq(0, 30, by=1/downsample))
#' 
#' @export
dHRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  C1 <- (1/(a1*b1))^a1
  C2 <- c*(1/(a2*b2))^a2
  A1 <- a1*t^(a1 - 1)*exp(-(t - a1*b1)/b1)
  A2 <- a2*t^(a2 - 1)*exp(-(t - a2*b2)/b2)
  B1 <- t^a1 / b1 * exp(-(t - a1*b1)/b1)
  B2 <- t^a2 / b2 * exp(-(t - a2*b2)/b2)
  C1*(A1 - B1) - C2 * (A2 - B2)
}

#' Calculate second derivative of the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response. Default: \code{6}
#' @param b1 response dispersion. Default: \code{0.9}
#' @param a2 delay of undershoot. Default: \code{12}
#' @param b2 dispersion of undershoot. Default: \code{0.9}
#' @param c scale of undershoot. Default: \code{0.35}
#'
#' @return d2HRF vector corresponding to time
#' 
#' @examples
#' downsample <- 100
#' d2HRF(seq(0, 30, by=1/downsample))
#' 
#' @export
d2HRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  C1 <- (1/(a1*b1))^a1
  C2 <- c*(1/(a2*b2))^a2
  dA1 <- a1 * ((a1 - 1) - t/b1) * t^(a1-2) * exp(-(t - a1*b1)/b1)
  dB1 <- (1/b1) * (a1 - (t / b1)) * t^(a1 - 1) * exp(-(t - a1*b1)/b1)
  dA2 <- a2 * ((a2 - 1) - t/b2) * t^(a2-2) * exp(-(t - a2*b2)/b2)
  dB2 <- (1/b2) * (a2 - (t / b2)) * t^(a2 - 1) * exp(-(t - a2*b2)/b2)
  C1 * (dA1 - dB1) - C2 * (dA2 - dB2)
}
