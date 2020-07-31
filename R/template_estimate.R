#' Estimate Template for Template or Diagnostic ICA
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries (*.dtseries.nii) for template estimation
#' @param cifti_fnames2 Vector of file paths of "retest" CIFTI-format fMRI timeseries (*.dtseries.nii) for template estimation.  Must be from the same subjects and in the same order as cifti_fnames.  If none specified, will create pseudo test-retest data from single session.
#' @param GICA_fname File path of CIFTI-format group ICA maps (ending in .d*.nii)
#' @param inds Indicators of which group ICs to include in template. If NULL, use all group ICs.
#' @param brainstructures Vector of brainstructures to include ('left','right','surface').  Default is c('left','right').
#' @param verbose If TRUE, display progress updates
#'
#' @return
#' @export
#' @importFrom ciftiTools read_cifti flatten_cifti
#'
estimate_template.cifti <- function(cifti_fnames, cifti_fnames2=NULL, GICA_fname, inds=NULL, brainstructures=c('left','right'), verbose=TRUE){

  if(verbose) cat('\n Reading in GICA result')
  GICA_cifti <- read_cifti(GICA_fname, brainstructures=brainstructures)
  #may need to save some info from GICA_cifti to map templates back to surfaces for visualization

  if(verbose) cat('\n Flattening GICA')
  GICA <- flatten_cifti(GICA_cifti, brainstructures='everything') #Damon I know this might need to be changed if the brainstructures argument is not everything
  L0 <- ncol(GICA)
  V <- nrow(GICA)
  wall_mask <- rep(FALSE, V)
  wall_mask[is.nan(GICA[,1])] <- TRUE
  wall_mask[is.na(GICA[,1])] <- TRUE
  wall_mask[rowSums(GICA^2)==0] <- TRUE
  GICA2 <- GICA[!wall_mask,]

  if(verbose){
    cat(paste0('\n Number of data locations: ',V))
    cat(paste0('\n Number of group ICs: ',L0))
  }

  L <- L0
  if(!is.null(inds)){
    if(any(!(inds %in% 1:L0))) stop('Invalid entries in inds argument.')
    L <- length(inds)
  } else {
    inds <- 1:L0
  }

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ',L))
    cat(paste0('\n Number of training subjects: ',N))
  }

  retest <- !is.null(cifti_fnames2)
  if(retest){
    if(length(cifti_fnames) != length(cifti_fnames2)) stop('If provided, cifti_fnames2 must have same length as cifti_fnames and be in the same subject order.')
  }

  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  missing_data <- 0
  for(ii in 1:N){

    if(verbose) cat(paste0('\n Reading in data for subject ',ii,' of ',N))

    #read in BOLD
    fname_ii <- cifti_fnames[ii]
    if(!file.exists(fname_ii)) {
      missing_data <- missing_data + 1
      if(verbose) cat(paste0('\n Data not available'))
      next
    }
    BOLD1_ii <- read_cifti(fname_ii, brainstructures=brainstructures, flat=TRUE)

    if(nrow(BOLD1_ii) != nrow(GICA2)) stop(paste0('The number of data locations in GICA and timeseries data from subject ',ii,' do not match.'))
    ntime <- ncol(BOLD1_ii)

    #read in BOLD retest data OR create pseudo test-retest data
    if(!retest){
      part1 <- 1:round(ntime/2)
      part2 <- setdiff(1:ntime, part1)
      BOLD2_ii <- BOLD1_ii[,part2]
      BOLD1_ii <- BOLD1_ii[,part1]
    } else {
      #read in BOLD from retest
      fname_ii <- cifti_fnames2[ii]
      if(!file.exists(fname_ii)) {
        missing_data <- missing_data + 1
        if(verbose) cat(paste0('\n Data not available'))
        next
      }
      BOLD2_ii <- read_cifti(fname_ii, brainstructures=brainstructures, flat=TRUE)
    }

    #perform dual regression on test and retest data
    DR1_ii <- dual_reg(t(BOLD1_ii), t(GICA2))$S
    DR2_ii <- dual_reg(t(BOLD2_ii), t(GICA2))$S
    DR1[ii,,!wall_mask] <- DR1_ii[inds,]
    DR2[ii,,!wall_mask] <- DR2_ii[inds,]
  }

  # ESTIMATE MEAN

  if(verbose) cat('\n Estimating Template Mean')
  mean1 <- apply(DR1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(DR2, c(2,3), mean, na.rm=TRUE)
  template_mean <- (mean1 + mean2)/2;

  # ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE

  # total variance
  if(verbose) cat('\n Estimating Total Variance')
  var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  var_tot <- (var_tot1 + var_tot2)/2

  # noise (within-subject) variance
  if(verbose) cat('\n Estimating Within-Subject Variance')
  DR_diff = DR1 - DR2;
  var_noise <- (1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE)

  # signal (between-subject) variance
  if(verbose) cat('\n Estimating Template (Between-Subject) Variance')
  template_var <- var_tot - var_noise
  template_var[template_var < 0] <- 0

  template <- list(mean = template_mean, var = template_var)

  return(template)


}
