#' Estimate Template for Template or Diagnostic ICA
#' 
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries 
#'  (*.dtseries.nii) for template estimation.
#' @param cifti_fnames2 (Optional) Vector of file paths of "retest" 
#'  CIFTI-format fMRI timeseries (*.dtseries.nii) for template estimation.
#'  Must be from the same subjects and in the same order as \code{cifti_fnames}.
#'  If none specified, will create pseudo test-retest data from single session.
#' @param GICA File path or a "cifti" object for the CIFTI-format group ICA 
#'  maps (ending in .d*.nii).
#' @param inds Indicators of which group ICs to include in template.
#'  \code{NULL} (default) use all group ICs.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface), and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Default: \code{c("left","right")} (brain surface only).
#' @param verbose If \code{TRUE}, display progress updates.
#'
#' @return A list 
#' @export
#' @importFrom ciftiTools read_cifti_flat match_input is.fname
#'
template_estimate_cifti <- function(
  cifti_fnames, cifti_fnames2=NULL, GICA, 
  inds=NULL, brainstructures=c('left','right'), verbose=TRUE){

  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )

  # ----------------------------------------------------------------------------
  # Read and flatten GICA file. ------------------------------------------------
  # ----------------------------------------------------------------------------

  # [TO DO]: Remove this line
  GICA <- '~/Dropbox/RESEARCH/TemplateICA/templates/HCP_surface_ICA25/HCP25.dtseries.nii'

  if (is.fname(GICA)) {
    GICA <- read_cifti(GICA, brainstructures=brainstructures)
  }

  # May need to save some info from GICA CIFTI to map templates back to 
  #   surfaces for visualization.
  GICA_cifti <- GICA

  # Flatten GICA (any component may be NULL).
  GICA <- flatten_cifti(GICA, brainstructures=brainstructures)
  L0 <- ncol(GICA); V <- nrow(GICA)
  # [TO DO]: medial wall indicator function in ciftiTools
  #   rowSums < EPS or == 0 exactly?
  #   handling of single- or few-column CIFTIs that may have zero values.
  wall_mask <- (is.nan(GICA[,1]) | rowSums(GICA)==0) #medial wall indicator

  if (is.null(inds)) {
    inds <- 1:L0
  } else {
    if(any(!(inds %in% 1:L0))) stop('Invalid entries in inds argument.')
  }
  L <- length(inds)

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\nNumber of...'))
    cat(paste0('\n\tdata locations: ', V))
    cat(paste0('\n\tgroup ICs: ', L0))
    cat(paste0('\n\ttemplate ICs: ', L))
    cat(paste0('\n\ttraining subjects: ', N))
  }

  retest <- !is.null(cifti_fnames2)
  if(retest){
    if(length(cifti_fnames) != length(cifti_fnames2)) {
      stop(paste0(
        'If provided, cifti_fnames2 must have same length as',
        'cifti_fnames and be in the same subject order.'
      ))
    }
  }

  # ----------------------------------------------------------------------------
  # Perform dual regression on (pseudo) test-retest data for each subject. -----
  # ----------------------------------------------------------------------------

  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  for(ii in 1:N){

    # Read in flattened BOLD.
    BOLD1_ii <- flatten_cifti(
      read_cifti(cifti_fnames[ii], brainstructures=brainstructures), brainstructures=brainstructures
    )

    ## Damon, can you see if this actually removes BOLD_cifti_ii from memory?
    ## Damon: as far as I know, rm() and gc() are the best you can do to try to 
    ##  free up memory, without running code in a separate R script. It's not
    ##  guaranteed that will work, but supposedly, if R needs the memory,
    ##  it will delete it.
    gc()

    if (nrow(BOLD1_ii) != V) {
      stop(paste0(
        'The number of data locations in GICA (', nrow(BOLD1_ii), ') ', 
        'and timeseries data from subject ', ii, ' (', V, ') ', 
        'do not match.'
      ))
    }
    # Temporary comment: replaced ntime with T_, similar to clever.
    #   Don't use T because that can mean TRUE.
    T_ <- ncol(BOLD1_ii) 

    # Read in BOLD retest data OR create pseudo test-retest data.
    if(!retest){
      part1 <- 1:round(T_/2)
      part2 <- setdiff(1:T_, part1)
      BOLD2_ii <- BOLD1_ii[,part2]
      BOLD1_ii <- BOLD1_ii[,part1]
    } else {
      # Read in flattened BOLD from retest.
      BOLD2_ii <- flatten_cifti(
        read_cifti(cifti_fnames2[ii], brainstructures=brainstructures), brainstructures=brainstructures
      )
    }

    gc()

    # Perform dual regression on test and retest data.
    DR1_ii <- dual_reg(t(BOLD1_ii[!wall_mask,]), t(GICA[!wall_mask,]))$S
    DR2_ii <- dual_reg(t(BOLD2_ii[!wall_mask,]), t(GICA[!wall_mask,]))$S
    DR1[ii,,!wall_mask] <- DR1_ii[inds,]
    DR2[ii,,!wall_mask] <- DR2_ii[inds,]
  }

  # ----------------------------------------------------------------------------
  # Estimate mean and variance. ------------------------------------------------
  # ----------------------------------------------------------------------------

  template_estimate(DR1, DR2)
}

#' Estimate template mean and variance from test-retest data
#' 
#' Estimates the template mean and variance based on repeated 
#'  estimates of subject-level ICs. These can be noisy as long
#'  as enough subjects are included, since this function estimates 
#'  and removes the noise variance.  
#' 
#' The estimates x1 and x2 should be independent, i.e. from repeated 
#'  fMRI runs or from the same fMRI run, split down the middle.
#' 
#' The L ICs of x1 and x2 must be matched. The easiest way to 
#'  do this is to use a common set of group ICs for both sets of visits or 
#'  sessions.
#' 
#' @param x1 Test data/first measurement: 
#'  an N (subjects) x L (ICs) x V (data locations/brainordinates) array.
#' @param x2 Retest data/second measurement: 
#'  an N (subjects) x L (ICs) x V (data locations/brainordinates) array.
#' 
#' @return A list of two components: "mean" and "variance", each storing the
#'  corresponding L x V template.
#' 
#' @export
template_estimate <- function(x1, x2){

  template <- list(mean=NULL, variance=NULL)

  # Mean.
  mean1 <- apply(x1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(x2, c(2,3), mean, na.rm=TRUE)
  template$mean <- (mean1 + mean2) / 2

  # Total variance.
  var_tot1 <- apply(x1, c(2,3), var, na.rm=TRUE)
  var_tot2 <- apply(x2, c(2,3), var, na.rm=TRUE)
  var_tot <- (var_tot1 + var_tot2)/2

  # Noise (within-subject) variance.
  mean_diff <- x1 - x2
  var_noise <- apply(mean_diff, c(2,3), var, na.rm=TRUE) / 2

  # Signal (between-subject) variance.
  template$variance <- var_tot - var_noise
  template$variance[template$variance < 0] <- 0

  template
}