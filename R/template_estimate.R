#' Estimate Template for Template or Diagnostic ICA
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries (*.dtseries.nii) for template estimation.
#' @param cifti_fnames2 (Optional) Vector of file paths of "retest" CIFTI-format fMRI timeseries (*.dtseries.nii) for template estimation.
#'  Must be from the same subjects and in the same order as \code{cifti_fnames}.
#'  If none specified, will create pseudo test-retest data from single session.
#' @param GICA_fname File path of CIFTI-format group ICA maps (ending in .d*.nii)
#' @param inds Indicators of which group ICs to include in template.
#'  If NULL, use all group ICs.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface), and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Default: \code{c("left","right")} (brain surface only).
#' @param verbose If \code{TRUE}, display progress updates.
#'
#' @return
#' @export
#' @importFrom ciftiTools read_cifti_flat match_input
#'
template_estimate_cifti <- function(cifti_fnames, cifti_fnames2=NULL, GICA_fname, inds=NULL, brainstructures=c('left','right'), verbose=TRUE){

  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical"),
    user_value_label="brainstructures"
  )

  GICA_fname <- '~/Dropbox/RESEARCH/TemplateICA/templates/HCP_surface_ICA25/HCP25.dtseries.nii'
  GICA_cifti <- read_cifti(GICA_fname, brainstructures=brainstructures)
  #may need to save some info from GICA_cifti to map templates back to surfaces for visualization

  # Flatten GICA (any component may be NULL).
  GICA <- rbind(GICA_cifti$CORTEX_LEFT, GICA_cifti$CORTEX_RIGHT)#, GICA_cifti$SUBCORT$VOL)
  L0 <- ncol(GICA)
  V <- nrow(GICA)
  wall_mask <- (is.nan(GICA[,1]) | rowSums(GICA)==0) #medial wall indicator

  if(verbose){
    cat(paste0('\n Number of data locations: ',V))
    cat(paste0('\n Number of group ICs: ',L0))
  }

  if (is.null(inds)) {
    inds <- 1:L0
  } else {
    if(any(!(inds %in% 1:L0))) stop('Invalid entries in inds argument.')
  }
  L <- length(inds)

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ', L))
    cat(paste0('\n Number of training subjects: ', N))
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

  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  for(ii in 1:N){

    #read in BOLD
    fname_ii <- cifti_fnames[ii]
    BOLD1_ii <- read_cifti(fname_ii, brainstructures=brainstructures)

    #flatten (any component may be NULL
    BOLD1_ii <- rbind(BOLD1_ii$CORTEX_LEFT, BOLD1_ii$CORTEX_RIGHT)
    #BOLD1_ii <- rbind(BOLD1_ii, BOLD_cifti_ii$VOL) #pending version of ciftiTools with flat volumetric data

    ##Damon, can you see if this actually removes BOLD_cifti_ii from memory?
    ##Damon: as far as I know, rm() and gc() are the best you can do to try to free up memory.
    ##But
    #rm(BOLD_cifti_ii)
    #gc()

    if(nrow(BOLD1_ii) != V) stop(paste0('The number of data locations in GICA and timeseries data from subject ',ii,' do not match.'))
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
      BOLD_cifti_ii <- read_cifti(fname_ii, brainstructures=brainstructures)

      #flatten
      BOLD2_ii <- BOLD_cifti_ii$CORTEX_LEFT
      BOLD2_ii <- rbind(BOLD2_ii, BOLD_cifti_ii$CORTEX_RIGHT)
      #BOLD2_ii <- rbind(BOLD2_ii, BOLD_cifti_ii$VOL) #pending version of ciftiTools with flat volumetric data
      rm(BOLD_cifti_ii) #Damon, can you see if this actually removes BOLD_cifti_ii from memory?
      gc()
    }

    #perform dual regression on test and retest data
    DR1_ii <- dual_reg(t(BOLD1_ii[!wall_mask,]), t(GICA[!wall_mask,]))$S
    DR2_ii <- dual_reg(t(BOLD2_ii[!wall_mask,]), t(GICA[!wall_mask,]))$S
    DR1[ii,,!wall_mask] <- DR1_ii[inds,]
    DR2[ii,,!wall_mask] <- DR2_ii[inds,]
  }

  # ESTIMATE MEAN

  mean1 <- apply(DR1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(DR2, c(2,3), mean, na.rm=TRUE)

  mean_meas1 = squeeze(nanmean(meas1,1));
  mean_meas2 = squeeze(nanmean(meas2,1));
  template_mean = (mean_meas1 + mean_meas2)/2;

  # ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE

  # total variance
  var_tot1 = squeeze(nanvar(meas1,0,1));
  var_tot2 = squeeze(nanvar(meas2,0,1));
  var_tot = (var_tot1 + var_tot2)/2;

  # noise (within-subject) variance
  meas_diff = meas1 - meas2;
  var_noise = (1/2)*squeeze(nanvar(meas_diff,0,1));

  # signal (between-subject) variance
  template_var = var_tot - var_noise;
  template_var(template_var < 0) = 0;




}
