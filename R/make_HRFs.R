#' Canonical HRF
#' 
#' The canonical hemodynamic response function for a given time vector. The code
#'  is adapted from \code{neuRosim::canonicalHRF}.
#' 
#' @param x The time vector (seconds)
#' 
#' @return The canonical HRF
#' 
#' @keywords internal
canonical_HRF <- function(x) {
  a1 <- 6; a2 <- 12
  b1 <- b2 <- 0.9
  c <- 0.35
  q1 <- a1 * b1
  q2 <- a2 * b2
  (x/q1)^a1 * exp(-(x - q1)/b1) - c*(x/q2)^a2 * exp(-(x - q2)/b2)
}

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
#'
#' @return Design matrix containing one HRF column for each task
#' 
#' @importFrom stats convolve
#' 
#' @export
make_HRFs <- function(onsets, TR, duration, downsample=100){

  K <- length(onsets) #number of tasks
  if(is.null(names(onsets))) task_names <- paste0('task', 1:K) else task_names <- names(onsets)

  nsec <- duration*TR; # Total time of experiment in seconds
  stimulus <- rep(0, nsec*downsample) # For stick function to be used in convolution
  HRF <- canonical_HRF(seq(0, 30, by=1/downsample))[-1] #canonical HRF to be used in convolution
  inds <- seq(TR*downsample, nsec*downsample, by = TR*downsample) # Extract EVs in a function of TR

  design <- matrix(NA, nrow=duration, ncol=K)
  colnames(design) <- task_names
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

    # Convolve boxcar with canonical HRF & add to design matrix
    HRF_k <- convolve(stimulus_k, rev(HRF), type='open')
    design[,k] <- HRF_k[inds]
  }

  return(design)
}
