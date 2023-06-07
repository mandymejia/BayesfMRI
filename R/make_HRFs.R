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

#' Make HRFs
#'
#' Create HRF design matrix columns from onsets and durations
#'
#' @param onsets \eqn{L}-length list in which the name of each element is the 
#'  name of the corresponding task, and the value of each element is a matrix of
#'  onsets (first column) and durations (second column) for each stimuli (each 
#'  row) of the corresponding task.
#'
#' @param TR Temporal resolution of the data, in seconds.
#' @param duration The number of volumes in the fMRI data. 
#' @param dHRF Set to \code{1} to add the temporal derivative of each column
#'  in the design matrix, \code{2} to add the second derivatives too, or
#'  \code{0} to not add any columns. Default: \code{1}.
#' @param dHRF_as Only applies if \code{dHRF > 0}. Model the temporal 
#'  derivatives as \code{"nuisance"} signals to regress out, \code{"tasks"}, or 
#'  \code{"auto"} to treat them as tasks unless the total number of columns in
#'  the design matrix (i.e. the total number of tasks, times `dHRF+1`), would be 
#'  \code{>=10}, the limit for INLA.
#' @param downsample Downsample factor for convolving stimulus boxcar or stick
#'  function with canonical HRF. Default: \code{100}.
#' @param verbose If applicable, print a message saying how the HRF derivatives
#'  will be modeled? Default: \code{FALSE}.
#'
#' @return List with the design matrix and/or the nuisance matrix containing the
#'  HRF-convolved stimuli as columns, depending on \code{dHRF_as}.
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
make_HRFs <- function(
  onsets, TR, duration, 
  dHRF=c(0, 1, 2),
  dHRF_as=c("auto", "nuisance", "task"),
  downsample=100, 
  verbose=FALSE){
  
  dHRF <- as.numeric(match.arg(as.character(dHRF), c("0", "1", "2")))
  if (dHRF == 0) {
    if (identical(dHRF_as, "nuisance") || identical(dHRF_as, "task")) {
      warning("`dHRF_as` is only applicable if `dHRF > 0`. If `dHRF == 0`, there's no need to specify `dHRF_as`.")
    }
  }
  dHRF_as <- match.arg(dHRF_as, c("auto", "nuisance", "task"))

  nK <- length(onsets) #number of tasks

  if (dHRF > 0 && dHRF_as=="auto") {
    nJ <- (dHRF+1) * nK # number of design matrix columns
    if (nJ > 5) {
      if (verbose) { 
        message("Modeling the HRF derivatives as nuisance signals.")
      }
      dHRF_as <- "nuisance"
    } else {
      if (verbose) {
        message("Modeling the HRF derivatives as tasks signals.")
      }
      dHRF_as <- "task"
    }
  }

  task_names <- if (is.null(names(onsets))) {
    task_names <- paste0('task', 1:nK)
  } else {
    names(onsets)
  }

  nsec <- duration*TR; # Total time of experiment in seconds
  stimulus <- rep(0, nsec*downsample) # For stick function to be used in convolution
  inds <- seq(TR*downsample, nsec*downsample, by = TR*downsample) # Extract EVs in a function of TR

  design <- nuisance <- vector("list", length=dHRF+1)

  for (dd in seq(0, dHRF)) {
    dname_dd <- switch(dd+1, "HRF", "dHRF", "d2HRF")
    theHRF_dd <- matrix(NA, nrow=duration, ncol=nK)
    colnames(theHRF_dd) <- paste0(task_names, "_", dname_dd)

    # canonical HRF to be used in convolution
    HRF_dd <- HRF(seq(0, 30, by=1/downsample), deriv=dd)[-1]

    for (kk in seq(nK)) {
      onsets_k <- onsets[[kk]][,1] #onset times in scans
      durations_k <- onsets[[kk]][,2] #durations in scans

      # Define stimulus function
      stimulus_k <- stimulus
      for (ii in seq(length(onsets_k))) {
        start_ii <- round(onsets_k[ii]*downsample)
        end_ii <- round(onsets_k[ii]*downsample + durations_k[ii]*downsample)
        stimulus_k[start_ii:end_ii] <- 1
      }

      # Convolve boxcar with canonical HRF & add to design2[[ii]] matrix
      HRF_k <- convolve(stimulus_k, rev(HRF_dd), type='open')
      theHRF_dd[,kk] <- HRF_k[inds]
    }

    if (dd > 0 && dHRF_as == "nuisance") {
      nuisance[[dd+1]] <- theHRF_dd
    } else {
      design[[dd+1]] <- theHRF_dd
    }
  }

  list(
    design=do.call(cbind, design), 
    nuisance=do.call(cbind, nuisance)
  )
}

#' Canonical (double-gamma) HRF
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
#' downsample <- 100
#' HRF(seq(0, 30, by=1/downsample))
#' 
#' @importFrom fMRItools is_1
#' @export
HRF <- function(t, deriv=0, a1 = 6,b1 = 0.9,a2 = 12,b2 = 0.9,c = 0.35) {

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