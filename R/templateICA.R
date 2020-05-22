#' Perform template independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (LxV matrix) template mean estimates, i.e. mean of empirical population prior for each of L independent components
#' @param template_var (LxV matrix) template variance estimates, i.e. between-subject variance of empirical population prior for each of L ICs
#' @param BOLD (TxV matrix) BOLD fMRI data matrix, where T is the number of volumes (time points) and V is the number of brain locations
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param mesh Either NULL (assume spatial independence) or an object of type \code{templateICA_mesh} created by \code{make_templateICA_mesh()} (spatial priors are assumed on each independent component)
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If TRUE, display progress of algorithm
#' @param kappa_init Starting value for kappa.  If NULL, starting value will be determined automatically.
#' @param dim_reduce_flag If FALSE, do not perform dimension reduction in spatial template ICA. Not applicable for standard template ICA.
#'
#' @return A list containing the estimated independent components S (a LxV matrix), their mixing matrix A (a TxL matrix), and the number of nuisance ICs estimated (Q_nuis)
#' @export
#' @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @import pesel
#' @importFrom ica icaimax
#' @importFrom stats optim
#'
templateICA <- function(template_mean, template_var, BOLD, scale=TRUE, mesh=NULL, maxQ=NULL, common_smoothness=TRUE, maxiter=100, epsilon=0.001, verbose=TRUE, kappa_init=NULL, dim_reduce_flag=TRUE){

  flag <- inla.pardiso.check()
  if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO is required for computational efficiency. See inla.pardiso().')
  inla.setOption(smtp='pardiso')

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of data locations
  L <- nrow(template_mean) #number of ICs

  #check that the number of data locations (nvox), time points (ntime) and ICs (L) makes sense
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(L > nvox) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if(L > ntime) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  #check that all arguments have consistent number of data locations (nvox) and ICs (L)
  if(ncol(template_mean) != nvox | ncol(template_var) != nvox) stop('template_mean, template_var and BOLD must have same number of data locations (columns), but they do not.')
  if(nrow(template_var) != L) stop('template_mean and template_var must have the same number of ICs (rows), but they do not.')

  #check that the supplied mesh object is of type templateICA_mesh
  do_spatial <- !is.null(mesh)
  if(!do_spatial){
    message('No mesh supplied: Using standard template ICA model, which assumes spatial independence. If this is not what you want, stop and supply a valid mesh. See help(make_templateICA_mesh).')
  } else if(class(mesh) != 'templateICA_mesh'){
    stop('mesh argument should be of class templateICA_mesh. See help(make_templateICA_mesh).')
  }

  if(!is.null(maxQ)){
    if(round(maxQ) != maxQ | maxQ <= 0) stop('maxQ must be NULL or a round positive number')
  }

  if(round(maxiter) != maxiter | maxiter <= 0) stop('maxiter must be a positive integer')

  if(!is.null(kappa_init)){
    if(length(kappa_init) != 1 | kappa_init <= 0) stop('kappa_init must be a positive scalar or NULL')
  }

  #check that maxQ makes sense
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
  if(class(common_smoothness) != 'logical' | length(common_smoothness) != 1) stop('common_smoothness must be a logical value')

  ### 1. ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)

  if(maxQ > L){

    #i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
    BOLD1 <- scale_BOLD(BOLD, scale=scale)
    DR1 <- dual_reg(BOLD1, template_mean)

    #ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
    fit <- DR1$A %*% DR1$S
    BOLD2 <- BOLD1 - fit #data without template ICs

    #iii. ESTIMATE THE NUMBER OF REMAINING ICS
    #pesel function expects nxp data and will determine asymptotic framework
    #here, we consider n=T (volumes) and p=V (vertices), and will use p-asymptotic framework
    pesel_BOLD2 <- pesel(BOLD2, npc.max=100, method='homogenous')
    Q2_hat <- pesel_BOLD2$nPCs #estimated number of nuisance ICs
    if(verbose) cat(paste0('ESTIMATING AND REMOVING ',Q2_hat,' NUISANCE COMPONENTS'))

    #iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX
    ICA_BOLD2 <- icaimax(t(BOLD2), nc=Q2_hat, center=TRUE)

    #v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD3
    fit <- ICA_BOLD2$M %*% t(ICA_BOLD2$S)
    BOLD3 <- BOLD1 - fit #original data without nuisance ICs

  } else {

  # USE ORIGINAL DATA, SCALED, SINCE WE ARE ASSUMING NO NUISANCE COMPONENTS
   BOLD3 <- scale_BOLD(BOLD, scale=scale) #center, and if scale=TRUE, scale

  }

  ### 2. PERFORM DIMENSION REDUCTION --> BOLD4

  if(dim_reduce_flag) if(verbose) cat('PERFORMING DIMENSION REDUCTION \n')

  dat_list <- dim_reduce(BOLD3, L)

  BOLD4 <- dat_list$data_reduced
  H <- dat_list$H
  Hinv <- dat_list$H_inv
  C_diag <- dat_list$C_diag

  ### 3. SET INITIAL VALUES FOR PARAMETERS

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD3, template_mean)

  # Keep?
  HA <- H %*% dat_DR$A #apply dimension reduction
  #HA <- orthonorm(HA)  #orthogonalize
  # sd_A <- sqrt(colVars(Hinv %*% HA)) #get scale of A (after reverse-prewhitening)
  # HA <- HA %*% diag(1/sd_A) #standardize scale of A
  theta0 <- list(A = HA)


  #initialize residual variance
  theta0$nu0_sq = dat_list$sigma_sq
  if(verbose) print(paste0('nu0_sq = ',round(theta0$nu0_sq,1)))


  ### 4. RUN EM ALGORITHM!

  #TEMPLATE ICA
  if(do_spatial) print('INITIATING WITH STANDARD TEMPLATE ICA')
  resultEM <- EM_templateICA.independent(template_mean, template_var, BOLD4, theta0, C_diag, maxiter=maxiter, epsilon=epsilon)
  resultEM$A <- Hinv %*% resultEM$theta_MLE$A
  class(resultEM) <- 'tICA'

  #SPATIAL TEMPLATE ICA
  if(do_spatial){

    #include regression-based estimate of A for tICA & save results
    tmp <- dual_reg(BOLD3, t(resultEM$subjICmean))
    resultEM$A_reg <- tmp$A
    resultEM_tICA <- resultEM

    if(dim_reduce_flag == FALSE){
      BOLD4 <- BOLD3
      C_diag <- rep(1, ntime)
      theta0$A <- dat_DR$A
    }

    #starting value for kappas
    if(is.null(kappa_init)){
      if(verbose) print('Using ML on tICA estimates to determine starting value for kappa')
      locs <- mesh$mesh$idx$loc[!is.na(mesh$mesh$idx$loc)]

      #organize data and replicates
      for(q in 1:L){
        #print(paste0('IC ',q,' of ',L))
        d_q <- resultEM$subjICmean[q,] - template_mean[q,]
        #d_q <- tmp[,q] - template_mean[q,]
        rep_q <- rep(q, length(d_q))
        D_diag_q <- sqrt(template_var[q,])
        if(q==1) {
          dev <- d_q
          rep <- rep_q
          D_diag <- D_diag_q
        } else {
          dev <- c(dev, d_q)
          rep <- c(rep, rep_q)
          D_diag <- c(D_diag, D_diag_q)
        }
      }

      #determine MLE of kappa
      print(system.time(opt <- optim(par=c(0,-20),
                     fn=loglik_kappa_est,
                     method='L-BFGS-B',
                     lower=c(-5,-20),
                     upper=c(1,Inf), #kappa usually less than 1, log(1)=0
                     delta=dev,
                     D_diag=D_diag,
                     mesh=mesh,
                     Q=L)))
      kappa_init <- exp(opt$par[1])

      # data_inla <- list(y = dev, x = rep(locs, L), repl=rep)
      # formula <- y ~ -1 + f(x, model = mesh$spde, replicate = repl)
      # result <- inla(formula, data = data_inla, verbose = FALSE)
      # result_spde <- inla.spde.result(result, name='x', spde=mesh$spde)
      # kappa_init <- exp(result_spde$summary.log.kappa$mean)
      if(verbose) print(paste0('Starting value for kappa = ',paste(round(kappa_init,3), collapse=', ')))
    }

    theta0$kappa <- rep(kappa_init, L)

    #project BOLD and templates to mesh locations
    Amat <- mesh$A # n_orig x n_mesh matrix
    nmesh <- ncol(Amat)
    if(nrow(Amat) != nvox) stop('Mesh projection matrix (mesh$A) must have nvox rows (nvox is the number of data locations, the columns of BOLD, template_mean and template_var)')

    print('RUNNING SPATIAL TEMPLATE ICA')
    t000 <- Sys.time()
    resultEM <- EM_templateICA.spatial(template_mean, template_var, mesh, BOLD=BOLD4, theta0, C_diag, common_smoothness=common_smoothness, maxiter=maxiter, verbose=verbose, dim_reduce_flag=dim_reduce_flag)
    print(Sys.time() - t000)

    #project estimates back to data locations
    resultEM$subjICmean_mat <- matrix(resultEM$subjICmean, ncol=L)
  }


  if(dim_reduce_flag) {
    A <- Hinv %*% resultEM$theta_MLE$A
    resultEM$A <- A
  }

  #regression-based estimate of A
  tmp <- dual_reg(BOLD3, t(resultEM$subjICmean_mat))
  resultEM$A_reg <- tmp$A

  #for stICA, return tICA estimates also
  if(do_spatial){
    resultEM$result_tICA <- resultEM_tICA
    class(resultEM) <- 'stICA'
  }

  return(resultEM)

}




#' Identify areas of activation in each independent component map
#'
#' @param result Fitted stICA or tICA model object (of class stICA or tICA)
#' @param u Activation threshold, default = 0
#' @param alpha Significance level for joint PPM, default = 0.1
#' @param type Type of region.  Default is '>' (positive excursion region).
#' @param method_p If result is type tICA, the type of multiple comparisons correction to use for p-values, or NULL for no correction.  See \code{help(p.adjust)}.
#' @param verbose If TRUE, display progress of algorithm
#' @param which.ICs Indices of ICs for which to identify activations.  If NULL, use all ICs.
#' @param deviation If TRUE, identify significant deviations from the template mean
#'
#' @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#' @export
#' @import excursions
#' @importFrom stats pnorm p.adjust
#'
activations <- function(result, u=0, alpha=0.01, type=">", method_p='BH', verbose=FALSE, which.ICs=NULL, deviation=FALSE){

  if(!(type %in% c('>','<','!='))) stop("type must be one of: '>', '<', '!='")
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')
  if(!(class(result) %in% c('stICA','tICA'))) stop("result must be of class stICA or tICA")

  L <- ncol(result$A)
  if(is.null(which.ICs)) which.ICs <- 1:L
  if(min((which.ICs) %in% (1:L))==0) stop('Invalid entries in which.ICs')

  template_var <- result$template$var
  template_mean <- result$template$mean

  if(class(result) == 'stICA'){

    if(verbose) cat('Determining areas of activations based on joint posterior distribution of latent fields\n')

    nvox <- nrow(result$subjICmean_mat)
    L <- ncol(result$subjICmean_mat)

    #identify areas of activation in each IC
    active <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nvox, ncol=L)

     for(q in which.ICs){
      if(verbose) cat(paste0('.. IC ',q,' (',which(which.ICs==q),' of ',length(which.ICs),') \n'))
      inds_q <- (1:nvox) + (q-1)*nvox
      if(deviation){
        Dinv_mu_s <-  (result$subjICmean - as.vector(template_mean) - u)/as.vector(sqrt(template_var))
      } else {
        Dinv_mu_s <- (result$subjICmean - u)/as.vector(sqrt(template_var))
      }

      if(q==which.ICs[1]) {
        #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
        #we subtract u first since rescaling by D^(-1) would affect u too
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = result$Omega, type = type, u = 0, ind = inds_q))
        if(verbose) print(tmp)
        rho <- res_q$rho
      } else {
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = result$Omega, type = type, u = 0, ind = inds_q, rho=rho))
        if(verbose) print(tmp)
      }
      active[,q] <- res_q$E[inds_q]
      jointPPM[,q] <- res_q$F[inds_q]
      marginalPPM[,q] <- res_q$rho[inds_q]
      vars[,q] <- res_q$vars[inds_q]
    }

    result <- list(active=active, jointPPM=jointPPM, marginalPPM=marginalPPM, vars=vars, u = u, alpha = alpha, type = type, deviation=deviation)
  }

  if(class(result) == 'tICA'){

    if(verbose) cat('Determining areas of activations based on hypothesis testing at each location\n')

    nvox <- nrow(result$subjICmean)
    L <- ncol(result$subjICmean)

    if(deviation){
      t_stat <- as.matrix((result$subjICmean - template_mean - u)/sqrt(result$subjICvar))
    } else {
      t_stat <- as.matrix((result$subjICmean - u)/sqrt(result$subjICvar))
    }

    if(type=='>') pvals <- 1-pnorm(t_stat)
    if(type=='<') pvals <- pnorm(t_stat)
    if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

    if(verbose) cat(paste0('Correcting for multiple comparisons with method ', method_p))

    pvals_adj <- active <- matrix(NA, nrow=nvox, ncol=L)
    if(is.null(method_p)) method_p <- 'none'
    for(q in which.ICs){
      pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
      active[,q] <- (pvals_adj[,q] < alpha)
    }

    result <- list(active = active, pvals = pvals, pvals_adj = pvals_adj, tstats = t_stat, vars = result$subjICvar, u = u, alpha = alpha, method_p = method_p, deviation=deviation)
  }

  return(result)
}


