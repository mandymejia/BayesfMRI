#' Perform diagnostic independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (A list of G matrices, each VxL) template mean estimates for each group 1 to G
#' @param template_var (A list of G matrices, each VxL) template variance estimates for each group 1 to G
#' @param BOLD (VxT matrix) BOLD fMRI data matrix, where T is the number of volumes (time points) and V is the number of brain locations
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If TRUE, display progress of algorithm
#'
#' @return A list containing the posterior probabilities of group membership, the estimated independent components S (a VxL matrix), their mixing matrix A (a TxL matrix), the number of nuisance ICs estimated (Q_nuis)
#' @export
#' @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @import pesel
#' @importFrom stats optim
#' @importFrom abind abind
#' @importFrom matrixStats rowVars
#'
diagnosticICA <- function(template_mean,
                        template_var,
                        BOLD,
                        scale=TRUE,
                        maxQ=NULL,
                        maxiter=100,
                        epsilon=0.001,
                        verbose=TRUE){

  if(!is.list(template_mean)) stop('template_mean must be a list')
  if(!is.list(template_var)) stop('template_var must be a list')
  if(length(template_mean) != length(template_var)) stop('template_mean and template_var must have the same length')
  G <- length(template_mean)

  if(!is.matrix(BOLD)) stop('BOLD must be a matrix')
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of data locations
  L <- ncol(template_mean[[1]]) #number of ICs

  if(verbose) cat(paste0('Length of timeseries: T = ',ntime,'\n'))
  if(verbose) cat(paste0('Number of voxels/vertices: V = ',nvox,'\n'))
  if(verbose) cat(paste0('Number of ICs: L = ',L,'\n'))
  if(verbose) cat(paste0('Number of groups: G = ',G,'\n'))

  #check that each element of template_mean and template_var is a matrix
  #check that the dimensions of template_mean and template_var elements are ok
  if(any(sapply(template_mean, is.matrix)==FALSE)) stop('Each element of template_mean must be an VxL matrix')
  if(any(sapply(template_var, is.matrix)==FALSE)) stop('Each element of template_var must be an VxL matrix')
  if(any(sapply(template_mean, dim)[2,] != L) | any(sapply(template_mean, dim)[1,] != nvox)) stop('Each element of template_mean must be an VxL matrix')
  if(any(sapply(template_var, dim)[2,] != L) | any(sapply(template_var, dim)[1,] != nvox)) stop('Each element of template_mean must be an VxL matrix')

  #check that the number of data locations (nvox), time points (ntime) and ICs (L) makes sense
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(L > nvox) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if(L > ntime) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  if(round(maxiter) != maxiter | maxiter <= 0) stop('maxiter must be a positive integer')

  #check that maxQ makes sense
  if(!is.null(maxQ)){ if(round(maxQ) != maxQ | maxQ <= 0) stop('maxQ must be NULL or a round positive number') }
  if(is.null(maxQ)) maxQ <- ntime
  if(maxQ < L){
    warning('maxQ must be at least L.  Setting maxQ=L.')
    maxQ <- L
  }
  if(maxQ > ntime){
    warning('maxQ must be no more than T.  Setting maxQ = T.')
    maxQ <- ntime
  }

  if(class(scale) != 'logical' | length(scale) != 1) stop('scale must be a logical value')


  ### IDENTIFY AND REMOVE ANY BAD VOXELS/VERTICES

  keep <- rep(TRUE, nvox)
  for(g in 1:G){
    keep[rowSums(is.nan(template_mean[[g]])) > 0] <- FALSE
    keep[rowSums(is.na(template_mean[[g]])) > 0] <- FALSE
    keep[rowSums(is.nan(template_var[[g]])) > 0] <- FALSE
    keep[rowSums(is.na(template_var[[g]])) > 0] <- FALSE
  }
  keep[rowSums(is.nan(BOLD)) > 0] <- FALSE
  keep[rowSums(is.na(BOLD)) > 0] <- FALSE
  keep[rowVars(BOLD) == 0] <- FALSE

  if(sum(!keep) > 0){
    template_mean_orig <- template_mean
    template_var_orig <- template_var
    nvox <- sum(keep)
    if(verbose) cat(paste0('Excluding ',sum(!keep),' bad (NA, NaN or flat) voxels/vertices from analysis.\n'))
    for(g in 1:G){
      template_mean[[g]] <- template_mean[[g]][keep,]
      template_var[[g]] <- template_var[[g]][keep,]
    }
    BOLD <- BOLD[keep,]
  }


  ### 1. ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)

  template_mean_avg <- apply(abind(template_mean, along=3), c(1,2), mean)
  # template_var1_avg <- apply(abind(template_var, along=3), c(1,2), mean) #avg within-group var
  # template_var2_avg <- apply(abind(template_mean, along=3), c(1,2), var) #between-group var (the unbiased N-1 version), equals 1/2*(y1-y2)^2 for N=2
  # template_var_avg <- template_var1_avg + template_var2_avg

  if(maxQ > L){

    #i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
    BOLD1 <- scale_BOLD(BOLD, scale=scale)
    DR1 <- dual_reg(BOLD1, template_mean_avg)

    #ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
    fit <- DR1$A %*% DR1$S
    BOLD2 <- BOLD1 - t(fit) #data without template ICs

    #iii. ESTIMATE THE NUMBER OF REMAINING ICS
    #pesel function expects nxp data and will determine asymptotic framework
    #here, we consider n=T (volumes) and p=V (vertices), and will use p-asymptotic framework
    if(verbose) cat(paste0('DETERMINING NUMBER OF NUISANCE COMPONENTS\n'))
    pesel_BOLD2 <- pesel(BOLD2, npc.max=maxQ-L, method='homo')
    Q2_hat <- pesel_BOLD2$nPCs #estimated number of nuisance ICs

    #iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX
    #if(verbose) cat(paste0('ESTIMATING AND REMOVING ',Q2_hat,' NUISANCE COMPONENTS\n'))
    #ICA_BOLD2 <- icaimax(BOLD2, nc=Q2_hat, center=TRUE)

    #iv. INSTEAD OF ESTIMATING ICS, JUST ESTIMATE PCS!
    #THE RESIDUAL (BOLD3) IS THE EXACT SAME BECAUSE THE ICS ARE JUST A ROTATION OF THE PCS
    #IF THE NUISANCE ICS ARE NOT OF INTEREST, CAN TAKE THIS APPROACH
    svd_BOLD2 <- svd(t(BOLD2) %*% BOLD2, nu=Q2_hat, nv=0)
    vmat <- diag(1/svd_BOLD2$d[1:Q2_hat]) %*% t(svd_BOLD2$u) %*% t(BOLD2)
    fit <- svd_BOLD2$u %*% diag(svd_BOLD2$d[1:Q2_hat]) %*% vmat

    #v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD3
    #fit <- ICA_BOLD2$M %*% t(ICA_BOLD2$S)
    BOLD3 <- BOLD1 - t(fit) #original data without nuisance ICs

  } else {

    # USE ORIGINAL DATA, SCALED, SINCE WE ARE ASSUMING NO NUISANCE COMPONENTS
    BOLD3 <- scale_BOLD(BOLD, scale=scale) #center, and if scale=TRUE, scale

  }

  ### 2. PERFORM DIMENSION REDUCTION --> BOLD4

  dat_list <- dim_reduce(BOLD3, L)

  BOLD4 <- dat_list$data_reduced
  H <- dat_list$H
  Hinv <- dat_list$H_inv
  C_diag <- dat_list$C_diag
  C_diag <- C_diag * (dat_list$sigma_sq) #(nu^2)HH' in paper

  ### 3. SET INITIAL VALUES FOR PARAMETERS

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD3, template_mean_avg)
  HA <- H %*% dat_DR$A #apply dimension reduction
  theta0 <- list(A = HA)
  theta00 <- theta0
  theta00$nu0_sq = dat_list$sigma_sq

  ### 4. RUN EM ALGORITHM!


  # if(verbose) cat('RUNNING STANDARD TEMPLATE ICA TO INITIALIZE\n')
  # resultEM0 <- EM_templateICA.independent(t(template_mean_avg), t(template_var_avg), t(BOLD4), theta00, C_diag, maxiter=maxiter, epsilon=epsilon, verbose=verbose)
  # theta0$A <- resultEM0$theta_MLE$A

  if(verbose) cat('BEGINNING EM ALGORITHM\n')
  resultEM <- EM_diagnosticICA(template_mean, template_var, BOLD4, theta0, C_diag, maxiter=maxiter, epsilon=epsilon, verbose=verbose)
  resultEM$A <- Hinv %*% resultEM$theta_MLE$A
  class(resultEM) <- 'dICA'

  # #regression-based estimate of A
  # tmp <- dual_reg(BOLD3, resultEM$subjICmean)
  # resultEM$A_reg <- tmp$A

  resultEM$keep <- keep

  #map estimates & templates back to original locations
  if(sum(!keep)>0){
    #estimates
    subjICmean <- subjICvar <- matrix(nrow=length(keep), ncol=L)
    subjICmean[keep,] <- resultEM$subjICmean
    subjICvar[keep,] <- resultEM$subjICvar
    resultEM$subjICmean <- subjICmean
    resultEM$subjICvar <- subjICvar
    #templates
    resultEM$template_mean <- template_mean_orig
    resultEM$template_var <- template_var_orig
  }

  return(resultEM)

}




