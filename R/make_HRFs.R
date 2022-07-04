#' Make HRFs
#'
#' Create HRF design matrix columns from onsets and durations
#'
#' @param onsets A matrix of onsets (first column) and durations (second column)
#'  for each task in SECONDS (set duration to zero for event related design),
#'  organized as a list where each element of the list corresponds to one task.
#'  Names of list indicate task names.
#' @param TR Temporal resolution of fMRI data, in SECONDS.
#' @param duration Length of fMRI timeseries, in SCANS.
#' @param downsample Downsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF
#' @param deriv This can take the value of 0, 1, or 2, and will use the HRF
#'   function, the first derivative of the HRF, or the second derivative of the
#'   HRF, respectively.
#'
#' @return Design matrix containing one HRF column for each task
#'
#' @importFrom stats convolve
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

    HRF_fn <- switch((dd+1), HRF, dHRF, d2HRF)
    HRF <- HRF_fn(seq(0, 30, by=1/downsample))[-1] #canonical HRF to be used in convolution

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
      HRF_k <- convolve(stimulus_k, rev(HRF), type='open')
      design[[dd]][,k] <- HRF_k[inds]
    }
  }
  design <- do.call(cbind, design)

  return(design)
}

#' Calculate the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response, default is 6
#' @param b1 response dispersion, default is 0.9
#' @param a2 delay of undershoot, default is 12
#' @param b2 dispersion of undershoot, default is 0.9
#' @param c scale of undershoot, default is 0.35
#'
#' @return HRF vector corresponding to time
#' @export
HRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  return(((t/(a1*b1))^a1) * exp(-(t-a1*b1)/b1) - c * ((t/(a2*b2))^a2) * exp(-(t - a2*b2)/b2))
}

#' Calculate first derivative of the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response, default is 6
#' @param b1 response dispersion, default is 0.9
#' @param a2 delay of undershoot, default is 12
#' @param b2 dispersion of undershoot, default is 0.9
#' @param c scale of undershoot, default is 0.35
#'
#' @return dHRF vector corresponding to time
#' @export
dHRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  C1 = (1/(a1*b1))^a1; C2 = c*(1/(a2*b2))^a2
  A1 = a1*t^(a1 - 1)*exp(-(t - a1*b1)/b1)
  A2 = a2*t^(a2 - 1)*exp(-(t - a2*b2)/b2)
  B1 = t^a1 / b1 * exp(-(t - a1*b1)/b1)
  B2 = t^a2 / b2 * exp(-(t - a2*b2)/b2)
  out <- C1*(A1 - B1) - C2 * (A2 - B2)
  return(out)
}

#' Calculate second derivative of the canonical (double-gamma) HRF
#'
#' @param t time vector
#' @param a1 delay of response, default is 6
#' @param b1 response dispersion, default is 0.9
#' @param a2 delay of undershoot, default is 12
#' @param b2 dispersion of undershoot, default is 0.9
#' @param c scale of undershoot, default is 0.35
#'
#' @return d2HRF vector corresponding to time
#' @export
d2HRF <- function(t,a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {
  C1 = (1/(a1*b1))^a1; C2 = c*(1/(a2*b2))^a2
  dA1 = a1 * ((a1 - 1) - t/b1) * t^(a1-2) * exp(-(t - a1*b1)/b1)
  dB1 = (1/b1) * (a1 - (t / b1)) * t^(a1 - 1) * exp(-(t - a1*b1)/b1)
  dA2 = a2 * ((a2 - 1) - t/b2) * t^(a2-2) * exp(-(t - a2*b2)/b2)
  dB2 = (1/b2) * (a2 - (t / b2)) * t^(a2 - 1) * exp(-(t - a2*b2)/b2)
  out <- C1 * (dA1 - dB1) - C2 * (dA2 - dB2)
  return(out)
}
