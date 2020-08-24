#' @name EM_templateICA
#' @rdname EM_templateICA
#'
#' @title EM Algorithms for Template ICA Models
#'
#' @param template_mean (QxV matrix) mean maps for each IC in template, where Q is the number of ICs, V is the number of data or mesh locations.
#' @param template_var  (QxV matrix) between-subject variance maps for each IC in template
#' @param mesh NULL for spatial independence model, otherwise an object of class "templateICA_mesh" containing the triangular mesh (see `help(make_mesh)`)
#' @param BOLD  (QxV matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at parameter values: A (QxQ mixing matrix), nu0_sq (residual variance from first level) and (for spatial model only) kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param maxiter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#' @param verbose If TRUE, display progress of algorithm
#'
#' @return  A list: theta (list of final parameter estimates), subICmean (estimates of subject-level ICs), subICvar (variance of subject-level ICs, for non-spatial model) or subjICcov (covariance matrix of subject-level ICs, for spatial model -- note that only diagonal and values for neighbors are computed), and success (flag indicating convergence (\code{TRUE}) or not (\code{FALSE}))
#'
#' @details \code{EM_templateICA.spatial} implements the expectation-maximization (EM) algorithm described in Mejia et al. (2019+) for estimating the subject-level ICs and unknown parameters in the template ICA model with spatial priors on subject effects.
#'
#' In both models, if original fMRI timeseries has covariance \eqn{\sigma^2 I_T}, the prewhitened timeseries achieved by premultiplying by (QxT) matrix \eqn{H} from PCA has diagonal covariance \eqn{\sigma^2HH'}, so C_diag is \eqn{diag(HH')}.
#'
#'
NULL

#' @rdname EM_templateICA
#' @export
#' @importFrom INLA inla.spde2.matern inla.qsolve
#' @importFrom Matrix Diagonal
#' @import SQUAREM
#'
EM_templateICA.spatial = function(template_mean, template_var, mesh, BOLD, theta0, C_diag, common_smoothness=TRUE, maxiter=100, epsilon=0.001, verbose=FALSE){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of data locations
  if(ntime > nvox) warning('More time points than data locations. Are you sure?')
  if(ncol(template_mean) != nvox) stop('Templates and BOLD must have the same number of data locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than data locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  if(class(mesh) != 'templateICA_mesh') stop('mesh argument should be of class templateICA_mesh. See help(make_mesh).')

	iter = 1
	theta = theta0
	success = 1

	template_var[template_var < 1e-6] = 1e-6 #to prevent problems when inverting covariance

	#pre-compute s0, D and D^{-1}*s0
	V <- ncol(template_mean)
	s0_vec = as.vector(t(template_mean))
	D_vec <- as.vector(sqrt(t(template_var))) #template_var is QxV
	D = Diagonal(V*Q, D_vec)
	Dinv_s0 <- inla.qsolve(Q = D, B=matrix(s0_vec, ncol=1), method='solve')

  ### REFINE STARTING VALUE FOR KAPPA

	if(verbose) cat('Refining starting value for kappa \n')

	# Determine direction of change:
	# Positive change --> search for kappa_max, set kappa_min to kappa1.
	# Negative change --> search for kappa_min, set kappa_max to kappa1.
	kappa_min <- kappa_max <- theta0$kappa[1]
	theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta0, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=verbose, update='kappa')
	kappa_diff0 <- theta1$kappa[1] - theta0$kappa[1]
	theta <- theta0

	kappa_diff <- kappa_diff0
	if(kappa_diff0 < 0){

	  if(verbose) cat('...Kappa decreasing, finding lower bound for kappa search \n ')

	  kappa_min <- kappa_min/2
	  while(kappa_diff < 0){
	    if(verbose) cat(paste0('... testing kappa = ',round(kappa_min,3),'\n '))
	    theta$kappa <- rep(kappa_min, Q)
	    theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=verbose, update='kappa')
	    kappa_diff <- theta1$kappa[1] - theta$kappa[1]
	    if(kappa_diff > 0) {
	      #set minimum and stop here
	      kappa_min <- theta1$kappa[1]
	      break
	    } else {
	      #set new value for kappa
	      kappa_min <- kappa_min/2
	    }
	  }
	} else if(kappa_diff0 > 0){

	  if(verbose) cat('...Kappa increasing, finding upper bound for kappa search \n ')

	  kappa_max <- kappa_max*2
	  while(kappa_diff > 0){
	    if(verbose) cat(paste0('... testing kappa = ',round(kappa_max, 3),'\n '))
	    theta$kappa <- rep(kappa_max, Q)
	    theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, update='kappa')
	    kappa_diff <- theta1$kappa[1] - theta$kappa[1]
	    if(kappa_diff < 0) {
	      #set maximum and stop here
	      kappa_max <- theta1$kappa[1]
	      break
	    } else {
	      #set new value for kappa
	      kappa_max <- kappa_max*2
	    }
	  }
	}

	#use binary search until convergence
	if(verbose) cat('...Starting binary search for starting value of kappa \n ')
	kappa_test <- (kappa_min + kappa_max)/2
	kappa_change <- 1
	while(kappa_change > epsilon){
	  if(verbose) cat(paste0('... testing kappa = ',round(kappa_test, 3),'\n '))
	  theta$kappa <- rep(kappa_test, Q)
	  theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, update='kappa')
	  kappa_diff <- theta1$kappa[1] - theta$kappa[1] #which direction is the estimate of kappa moving in?
	  if(kappa_diff > 0) {
	    kappa_min <- theta1$kappa[1]  #reset minimum to current value
	    kappa_test <- (theta1$kappa[1] + kappa_max)/2 #go halfway to max
	  } else {
	    kappa_max <- theta1$kappa[1]  #reset maximum to current value
	    kappa_test <- (theta1$kappa[1] + kappa_min)/2 #go halfway to min
	  }

	  #how much different is the next value of kappa to be tested versus the current one?
	  kappa_change <- abs((kappa_test - theta1$kappa[1])/theta1$kappa[1])
	}


	### RUN SQUAREM ALGORITHM UNTIL CONVERGENCE

	theta0 <- theta1 #last tested value of kappa0
	theta0$LL <- c(0,0)
	theta0_vec <- unlist(theta0[1:3]) #everything but LL
	names(theta0_vec)[1] <- 0 #store LL value in names of theta0_vec (required for squarem)

	t00000 <- Sys.time()
	result_squarem <- squarem(par=theta0_vec, fixptfn = UpdateThetaSQUAREM, objfn=LL_SQUAREM, control=list(trace=verbose, intermed=TRUE, tol=epsilon, maxiter=maxiter), template_mean, template_var, mesh, BOLD, C_diag, s0_vec, D, Dinv_s0, common_smoothness, verbose=FALSE)
	if(verbose) print(Sys.time() - t00000)

	path_A <- result_squarem$p.inter[,1:(Q^2)]
	path_kappa <- result_squarem$p.inter[,(Q^2+1)+(1:Q)]
	path_LL <- result_squarem$p.inter[,ncol(result_squarem$p.inter)]
	theta_path <- list(A=path_A, kappa=path_kappa, LL=path_LL)

	theta_MLE <- theta0
	theta_MLE$A <- matrix(result_squarem$par[1:(Q^2)], Q, Q)
	theta_MLE$kappa <- result_squarem$par[(Q^2+1)+(1:Q)]
	theta_MLE$LL <- as.numeric(names(result_squarem$par)[1])

	#success <- (result_squarem$convergence==0) #0 indicates convergence, 1 indicates failure to converge within maxiter
	numiter <- result_squarem$fpevals #number of parameter update steps (approximately 3x the number of SQUAREM iterations)

	### Compute final posterior mean of subject ICs
	if(verbose) cat('Computing final posterior mean of subject ICs \n')
	mu_Omega_s = UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta_MLE, C_diag, s0_vec, D, Dinv_s0, common_smoothness=common_smoothness, verbose=verbose, return_MAP=TRUE)

	subjICcov=(D %*% inla.qinv(mu_Omega_s$Omega_s) %*% D)

	result <- list(subjICmean=mu_Omega_s$mu_s,
	               subjICcov=subjICcov,
	               Omega = mu_Omega_s$Omega_s,
	               theta_MLE=theta_MLE,
	               theta_path=theta_path,
	               numiter=numiter,
	               squarem = result_squarem,
	               template_mean = template_mean,
	               template_var=template_var)
	return(result)
}

#' @rdname EM_templateICA
#' @export
EM_templateICA.independent = function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001, verbose){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of brain locations
  if(ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if(ncol(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  iter = 1
  theta = theta0
  success = 1
  template_var[template_var < 1e-6] = 1e-6 #to prevent problems when inverting covariance

  err = 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))

    t00 <- Sys.time()
    theta_new = UpdateTheta.independent(template_mean, template_var, BOLD, theta, C_diag, verbose=verbose)
    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    A_old = theta$A
    A_new = theta_new$A
    #2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
    A_change = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")

    nu0_sq_old = theta$nu0_sq
    nu0_sq_new = theta_new$nu0_sq
    nu0_sq_change = abs(nu0_sq_new - nu0_sq_old)/nu0_sq_old

    change = c(A_change, nu0_sq_change)
    err = max(change)
    change = format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for nu0_sq \n'))

    ### Move to next iteration
    theta <- theta_new
    iter = iter + 1
    if(iter > maxiter){
      success = 0;
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }

  ### Compute final posterior mean of subject ICs

  A = theta$A
  At_nu0Cinv = t(theta$A) %*% diag(1/(C_diag*theta$nu0_sq))
  At_nu0Cinv_A = At_nu0Cinv %*% theta$A
  miu_s = matrix(NA, nrow=Q, ncol=nvox)
  var_s = matrix(NA, nrow=Q, ncol=nvox)
  for(v in 1:nvox){
    y_v <- BOLD[,v]
    s0_v <- template_mean[,v]
    E_v_inv <- diag(1/template_var[,v])
    Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
    miu_s[,v] <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    var_s[,v] <- diag(Sigma_s_v)
  }

  result <- list(subjICmean=t(miu_s),
                 subjICvar=t(var_s),
                 theta_MLE=theta,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean = template_mean,
                 template_var = template_var)
  #names(result) <- c('subjICmean', 'subjICvar', 'theta_MLE', 'success_flag')
  return(result)
}


#' @name UpdateTheta
#' @rdname UpdateTheta
#'
#' @title Parameter Estimates in EM Algorithm for Template ICA Model
#'
#' @param template_mean (QxV matrix) mean maps for each IC in template
#' @param template_var (QxV matrix) between-subject variance maps for each IC in template
#' @param BOLD dimension-reduced fMRI data
#' @param mesh NULL for spatial independence model, otherwise an object of class "templateICA_mesh" containing the triangular mesh (see `help(make_templateICA_mesh)`)
#' @param theta A list of current parameter estimates (mixing matrix A, noise variance nu0_sq and (for spatial model) SPDE parameters kappa)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param s0_vec Vectorized template means
#' @param D Sparse diagonal matrix of template standard deviations
#' @param Dinv_s0 The inverse of D times s0_vec
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param verbose If TRUE, print progress updates for slow steps.
#' @param return_MAP If TRUE, returns the posterior mean and precision of the latent fields instead of the parameter estimates
#' @param update Which parameters to update. Either "all", "A" or "kappa".
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields
#'
NULL

#' @rdname UpdateTheta
#' @export
#' @importFrom stats optimize
#' @importFrom INLA inla.qsolve inla.qinv inla.setOption
#' @import Matrix
UpdateTheta.spatial = function(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, return_MAP=FALSE, update=c('all','kappa','A')){

  Q = nrow(template_mean)
  V = ncol(BOLD)
  ntime = nrow(BOLD)
  spde = mesh$spde

  #initialize new parameter values
  theta_new = theta

  #which parameters to update
  update_A <- (update[1] == 'all' | update[1] =='A')
  update_kappa <- (update[1] == 'all' | update[1] =='kappa')
  if(update_A) theta_new$A <- NA
  if(update_kappa) theta_new$kappa <- NA
  if(!update_A & !update_kappa & !return_MAP) stop('Please indicate at least one parameter to update (e.g. update="A") or set return_MAP=TRUE to return MAP estimates of S based on current parameter values.')

  #1. Compute mu_sy by solving for x in the system of equations Omega*x = m
  #   Compute m and Omega (using previous parameter estimates)
  #   Ingredients for Omega and m:
  #     s0_vec = "s_0" (long vector organized by ICs)
  #     y (long vector organized by locations)
  #     D, big diagonal matrix of template standard deviations
  #     R_inv, prior inverse correlation matrix for data locations (big block diagonal matrix, each block is sparse-ish)
  #     P (permutation matrix)

  if(verbose) cat('Computing Posterior Moments of S \n')

  y_vec = as.vector(BOLD)

  if(verbose) cat('...posterior precision \n') # less than 1 sec
  #Compute SPDE matrices (F, G, GFinvG) and Sigma_inv (QVxQV), a sparse block diagonal matrix
  stuff <- compute_R_inv(mesh, kappa=theta$kappa, C1=1/(4*pi))
  R_inv <- stuff$R_inv
  Fmat <- stuff$Fmat
  Gmat <- stuff$Gmat
  GFinvG <- stuff$GFinvG

  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(Q, V)

  #1. Compute mu_s
  if(verbose) cat('...posterior mean \n') #20 seconds with pardiso! (1 hour without)
  stuff <- compute_mu_s(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega
  Omega_inv_m <- stuff$Omega_inv_m

  if(verbose) print(summary(as.vector(mu_s)))

  if(return_MAP){
    return(list(mu_s = mu_s, Omega_s = Omega))
    stop()
  }

  if(update_A){
    if(verbose) cat('Updating A \n')
    t00 <- Sys.time()

    #Compute first matrix appearing in A-hat
    Pmu_vec <- P %*% mu_s
    yPmu <- matrix(0, nrow=ntime, ncol=Q)
    for(v in 1:V){
      inds_y <- (1:ntime) + (v-1)*ntime
      inds_Pmu <- (1:Q) + (v-1)*Q
      y_v <- y_vec[inds_y]
      Pmu_v <- Pmu_vec[inds_Pmu]
      yPmu_v <- outer(y_v, Pmu_v)
      yPmu <- yPmu + yPmu_v #sum up to later divide by V to average
    }

    #Compute second matrix appearing in A-hat
    Omega_PP <- P %*% Omega %*% t(P)
    D_PP <- P %*% D %*% t(P)

    if(verbose) cat('..Calculating only non-sparse covariance terms \n')

    #set up sparse indicator matrix for needed elements of Omega_PP^{-1}
    oneblock <- matrix(1, nrow=Q, ncol=Q)
    ind_mat <- bdiag_m(rep(list(oneblock), V))

    #augment Omega_PP with additional non-sparse locations
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    Omega_PP_aug <- Omega_PP + ind_mat
    Omega_PP_inv <- inla.qinv(Omega_PP_aug)

    T_mat <- matrix(0, Q, Q)
    for(v in 1:V){
      inds <- (1:Q) + (v-1)*Q

      #T1 matrix
      Omega_PP_inv_vv <- Omega_PP_inv[inds,inds]
      T1_vv <- D_PP[inds,inds] %*% Omega_PP_inv_vv %*% D_PP[inds,inds]

      #T2 matrix
      Pmu_v <- Pmu_vec[inds]
      T2_vv <- outer(Pmu_v, Pmu_v)

      #T(v,v)
      T_vv <- T1_vv + T2_vv
      T_mat <- T_mat + T_vv
    }

    A_hat <- t(solve(t(T_mat), t(yPmu))) # = yPmu %*% T_mat_inv
    theta_new$A <- A_hat

    LL1 <- 0
    for(v in 1:V){
      inds_y <- (1:ntime) + (v-1)*ntime
      inds_Pmu <- (1:Q) + (v-1)*Q
      y_v <- y_vec[inds_y]
      y_v_Cinv <- y_v * (1/C_diag)
      Pmu_v <- Pmu_vec[inds_Pmu] #t(v) vector in paper

      LL1_v <- 2 * t(y_v_Cinv) %*% A_hat %*% Pmu_v
      LL1 <- LL1 + as.numeric(LL1_v)
    }
    theta_new$LL[1] <- LL1

  }


  #keep value of nu0sq_hat from PCA-based dimension reduction
  nu0sq_hat <- theta$nu0_sq
  theta_new$nu0_sq <- nu0sq_hat[1]

  ##########################################
  ### E-STEP for kappa
  ##########################################

  if(update_kappa){

    if(verbose) cat('Updating kappa  \n')

    # Likelihood in terms of kappa_q's.
    # LL2_part1 = log(det(R_q_inv))
    # LL2_part2 = Tr(R_q_inv * Omega_inv_q) + Tr(R_q_inv * W_hat_q)
    # LL2_part3 = u_q' * R_q_inv * v_hat_q
    #             u_q = Dinv * s0
    #             v_q = 2 Omega_inv * m - Dinv * s0

    # Tr(R_q_inv * Omega_inv_q) --> Use inla.inv to compute necessary elements (non-zeros in R_q_inv) of Omega_inv_q for q=1,...,Q
    # Tr(R_q_inv * W_hat_q) --> W_hat_q = outer(Omega_inv*m,Omega_inv*m)_qq, where Omega_inv*m is known. Just calculate necessary elements (non-zeros in R_q_inv) of W_hat_q for q=1,...,Q

    if(verbose) cat('..Computing trace terms in log-likelihood of kappa \n')

    ### SET UP FOR PART 2

    #1. Generate Monte Carlo samples for estimation of Omega_hat^{-1}_qq

    # if(verbose) cat(paste0('....drawing ',nsamp,' Monte Carlo samples for covariance estimation \n'))
    #
    # nsamp <- 5000 #number of samples
    # musamp <- inla.qsample(nsamp, Q = Omega, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=8) #sample from N(0, Omega^(-1)) to estimate Omega^(-1) (diagonal blocks only)
    #9 sec for 500 samples with V*Q = 2530*3 = 12,642

    #2. Determine non-zero terms of R_q_inv

    if(verbose) cat('....calculating only non-sparse covariance terms \n')

    #set up estimation of only terms necessary for trace computation
    nonzero_Rinv <- as.matrix(1*(R_inv[1:V,1:V] != 0))
    #diag(nonzero_Rinv) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
    nonzero_cols <- which(R_inv[1:V,1:V] != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations

    #augment Omega with additional non-sparse locations
    ind_mat <- 1*(R_inv[1:V,1:V] != 0)
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    Omega_aug <- Omega + bdiag(rep(list(ind_mat), Q)) #technically not needed since Omega and R_inv have same sparsity structure

    #compute elements of Omega_inv that are non-sparse in Omega
    Omega_inv <- inla.qinv(Omega_aug)

    #if(verbose) cat('....setting up helper objects for trace computation \n') #15 sec (Q=16, V=5500)

    # #3. Set up matrices needed for computation of only necessary elements of Omega_hat^{-1}
    #
    # X <- matrix(1, nrow=nsamp, ncol=V) #placeholder for X in Omega^(-1) = XX'
    # bigX_left <- KhatriRao(diag(1, V), X) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
    # bigX_right <- KhatriRao(nonzero_Rinv, X) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
    # inds_left_X <- which(bigX_left != 0, arr.ind=TRUE)
    # inds_right_X <- which(bigX_right != 0, arr.ind=TRUE)
    #
    #4. Set up vectors needed for computation of only necessary elements of W_hat

    x <- matrix(1, nrow=1, ncol=V) #placeholder for x in M = xx'
    bigx_left <- KhatriRao(diag(1, V), x) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
    bigx_right <- KhatriRao(nonzero_Rinv, x) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
    inds_left_x <- which(bigx_left != 0, arr.ind=TRUE)
    inds_right_x <- which(bigx_right != 0, arr.ind=TRUE)

    if(verbose) cat('....computing necessary elements of RHS matrices in trace terms \n')
    #15 sec for Q=16, V=5500

    if(!common_smoothness) OplusW <- vector('list', length=Q) #Omega_inv_qq + W_hat_qq

    for(q in 1:Q){
      #if(verbose) cat(paste('......block ',q,' of ',Q,' \n'))
      inds_q <- (1:V) + (q-1)*V

      # COMPUTE OMEGA_INV[q,q] (NECESSARY ENTRIES ONLY)

      # Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE)
      #
      # #left-multiplication matrix
      # vals_left_X <- as.vector(Xctr_q)
      # bigX_left_q <- sparseMatrix(i = inds_left_X[,1], j = inds_left_X[,2], x = vals_left_X)
      # #right-multiplication matrix
      # X_repcols <- Xctr_q[,nonzero_cols]
      # vals_right_X <- as.vector(X_repcols)
      # bigX_right_q <- sparseMatrix(i = inds_right_X[,1], j = inds_right_X[,2], x = vals_right_X)
      # #multiply together
      # Omega_inv_qq <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

      Omega_inv_qq <- Omega_inv[inds_q,inds_q]

      # COMPUTE W_hat[q,q] = Omega_inv_m[q] * Omega_inv_m[q' (NECESSARY ENTRIES ONLY)

      # T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
      e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
      e_q[1,q] <- 1
      T_q <- kronecker(e_q, Diagonal(V))
      Omega_inv_m_q <- T_q %*% Omega_inv_m

      #left-multiplication matrix
      vals_left_x <- as.vector(Omega_inv_m_q)
      bigx_left_q <- sparseMatrix(i = inds_left_x[,1], j = inds_left_x[,2], x = vals_left_x)
      #right-multiplication matrix
      x_repcols <- t(Omega_inv_m_q)[,nonzero_cols]
      vals_right_x <- as.vector(x_repcols)
      bigx_right_q <- sparseMatrix(i = inds_right_x[,1], j = inds_right_x[,2], x = vals_right_x)
      #multiply together
      W_hat_qq <- t(bigx_left_q) %*% bigx_right_q

      # COMBINE OMEGA_INV[q,q] AND W_hat[q,q] --> two trace terms involving R_q_inv can be combined

      OplusW_qq <- Omega_inv_qq + W_hat_qq
      if(common_smoothness){
        if(q==1) OplusW <- OplusW_qq else OplusW <- OplusW + OplusW_qq
      } else {
        OplusW[[q]] <- OplusW_qq
      }
    }


    ### SET UP FOR PART 3

    u <- Dinv_s0
    v <- 2*Omega_inv_m - u

    # NUMERICALLY SOLVE FOR MLE OF KAPPA

    if(verbose) cat("..Performing numerical optimization \n")
    #~8 seconds for V=5500, Q=3 (comon smoothness)

    if(common_smoothness){
      kappa_opt <- optimize(LL2_kappa, lower=0, upper=5, maximum=TRUE,
                            mesh=mesh, OplusW=OplusW, u=u, v=v, Q=Q)  #Q=Q to indicate common smoothness model to the LL2_kappa function
      LL2 <- kappa_opt$objective
      kappa_opt <- rep(kappa_opt$maximum, Q)
    } else {
      kappa_opt <- LL2 <- rep(NA, Q)
      for(q in 1:Q){
        if(verbose) cat(paste('Optimization ',q,' of ',Q,' \n'))
        inds_q <- (1:V) + (q-1)*V
        kappa_opt_q <- optimize(LL2_kappa, lower=0, upper=5, maximum=TRUE,
                                spde=spde, OplusW=OplusW[[q]], u=u[inds_q], v=v[inds_q], Q=NULL)
        LL2[q] <- kappa_opt_q$objective
        kappa_opt[q] <- (kappa_opt_q$maximum)
      }
      LL2 <- sum(LL2)
    }

    theta_new$kappa <- kappa_opt
    theta_new$LL[2] <- LL2

  }

  # RETURN NEW PARAMETER ESTIMATES

  return(theta_new)

}
#' @rdname UpdateTheta
#' @export
UpdateTheta.independent = function(template_mean, template_var, BOLD, theta, C_diag, verbose){

  Q = nrow(BOLD)
  V = ncol(BOLD)

  #initialize new objects
  theta_new = list(A = matrix(NA, Q, Q), nu0_sq = NA)
  A_part1 = A_part2 = matrix(0, Q, Q) #two parts of product for A-hat (construct each looping over voxels)

  A = theta$A
  nu0_sq = theta$nu0_sq
  nu0C_inv = diag(1/(C_diag*nu0_sq)) #Sigma0_inv in matlab code
  At_nu0Cinv = t(A) %*% nu0C_inv
  At_nu0Cinv_A = At_nu0Cinv %*% A

  if(verbose) cat('Updating A \n')

  #store posterior moments for M-step of nu0_sq
  miu_s = matrix(NA, nrow=Q, ncol=V)
  miu_ssT = array(NA, dim=c(Q, Q, V))

  for(v in 1:V){

    y_v = BOLD[,v]
    s0_v = template_mean[,v]

    ##########################################
    ### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
    ##########################################

    E_v_inv = diag(1/template_var[,v])
    Sigma_s_v = solve(E_v_inv + At_nu0Cinv_A)
    miu_s_v = Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    miu_ssT_v = (miu_s_v %*% t(miu_s_v)) + Sigma_s_v #QxQ
    miu_s[,v] = miu_s_v #save for M-step of nu0_sq
    miu_ssT[,,v] = miu_ssT_v #save for M-step of nu0_sq

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 = A_part1 + y_v %*% t(miu_s_v) #QxQ
    A_part2 = A_part2 + miu_ssT_v #QxQ

  }

  #A_hat = orthonorm(A_part1 %*% solve(A_part2))
  A_hat = (A_part1 %*% solve(A_part2))

  ##########################################
  ### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
  ##########################################

  # cat('Updating Error Variance nu0_sq \n')
  #
  # #use A-hat or A?
  #
  # Cinv = diag(1/C_diag)
  # Cinv_A = Cinv %*% A_hat
  # At_Cinv_A = t(A_hat) %*% Cinv %*% A_hat
  # nu0sq_part1 = nu0sq_part2 = nu0sq_part3 = 0
  #
  # for(v in 1:V){
  #
  #   y_v = BOLD[,v]
  #   nu0sq_part1 = nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
  #   nu0sq_part2 = nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
  #   nu0sq_part3 = nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))
  # }
  #
  # nu0sq_hat = 1/(Q*V)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)

  nu0sq_hat <- theta$nu0_sq


  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  return(theta_new)
}



#' Computes posterior mean and precision matrix of s
#'
#' @param y_vec Vectorized, dimension-reduced fMRI data, grouped by locations. A vector of length QV.
#' @param D Sparse diagonal matrix of template standard deviations
#' @param Dinv_s0 The inverse of D times s0_vec
#' @param R_inv Estimate of inverse spatial correlation matrix (sparse)
#' @param theta List of current parameter estimates
#' @param P Permutation matrix for regrouping by locations (instead of by ICs.)
#' @param C_diag Diagonals of residual covariance of the first level model. A vector of length Q.
#'
#' @return A list containing the posterior mean \eqn{\mu} (mu) and precision \eqn{\Omega} (Omega) of s=(s1,...,sQ), along with the supporting vector m, where \eqn{\mu = \Omega^{-1}m}.
#'
#' @import Matrix
#' @importFrom INLA inla.qsolve
#'
compute_mu_s <- function(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag){

  ntime <- length(C_diag)
  Q <- ncol(theta$A)
  V <- nrow(P)/Q

  A <- theta$A
  nu0_sq <- theta$nu0_sq

  #set up B, C, and products thereof
  ones = Diagonal(V)
  B = kronecker(ones, A)
  nu0C_inv = kronecker(ones, diag(1/(C_diag*nu0_sq)))
  Pt_Bt_nu0C_inv = t(P) %*% t(B) %*% nu0C_inv

  #compute m (using current parameter estimates in theta)
  m1_vec <- D %*% Pt_Bt_nu0C_inv %*% y_vec
  m2_vec <- R_inv %*% Dinv_s0
  m_vec <- m1_vec + m2_vec

  # Compute Omega (using current parameter estimates in theta)
  Omega <- R_inv + D %*% Pt_Bt_nu0C_inv %*% B %*% P %*% D

  # Compute mu_s|y by solving for x in the system of equations Omega*x = m
  Omega_inv_m <- inla.qsolve(Q = Omega, B=m_vec, method='solve') # <-- slowest part (up to 1 hour for Q=16, V=5500), but with inla.setOption(smtp='pardiso') goes down to 20 seconds!
  mu_s <- D %*% Omega_inv_m

  return(list(mu=mu_s, m=m_vec, Omega=Omega, Omega_inv_m=Omega_inv_m))
}


#' Computes SPDE matrices (F, G and GFinvG) and prior precision matrix for S
#'
#' @param mesh Object of class "templateICA_mesh" containing the triangular mesh (see `help(make_templateICA_mesh)`)
#' @param kappa Current estimates of SPDE parameter kappa for each latent field
#' @param C1 Constant, equal to 1/(4*pi) for a 2-dimensional field with alpha=2
#'
#' @return A list containing R inverse and SPDE matrices
#'
#' @import Matrix
#' @importFrom stats var
#'
compute_R_inv <- function(mesh, kappa, C1=1/(4*pi)){

  Q <- length(kappa)

  #SPDE matrices, needed to construct R_q_inv
  spde = mesh$spde
  Fmat = spde$param.inla$M0
  Gmat = 1/2*(spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G

  if(var(kappa)==0) onekappa <- TRUE
  if(length(kappa)==1) onekappa <- TRUE
  if(onekappa) kappa <- kappa[1]

  #get inmesh and notinmesh indices
  Amat = mesh$A #n_loc x n_mesh
  V = ncol(mesh$A) #number of mesh locations
  inmesh = which(colSums(Amat) > 0)
  notinmesh = setdiff(1:V, inmesh)


  #set up R^{-1} (QVxQV) as a sparse block diagonal matrix
  if(onekappa){ #just compute Rinv once
    Qmat = C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
    Q11 = Qmat[inmesh,inmesh] # = Amat %*% Qmat %*% t(Amat)
    Q12 = Qmat[inmesh, notinmesh]
    Q21 = Qmat[notinmesh, inmesh]
    Q22 = Qmat[notinmesh,notinmesh]
    Q22_inv <- solve(Q22)
    R_q_inv = Q11 - (Q12 %*% Q22_inv %*% Q21)
    R_inv_list <- rep(list(R_q_inv), Q)
  } else { # compute Rinv block-wise
    R_inv_list = vector('list', Q)
    for(q in 1:Q){
      kappa_q = kappa[q]
      R_q_inv = C1 * (kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG)
      R_inv_list[[q]] = R_q_inv
    }
  }
  R_inv <- bdiag(R_inv_list)

  return(list(R_inv=R_inv, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG))
}

#' Creates permutation matrix P to regroup elements of s by locations instead of by ICs
#'
#' @param Q The number of template ICs
#' @param V The number of spatial locations
#'
#' @return P Permutation matrix size QVxQV
#'
#' @details If s=(s1,...,sQ) is grouped by ICs 1,...Q, then Ps=(s(1),...,s(V)) is grouped by locations 1,...,V
#'
#' @importFrom Matrix sparseMatrix
#'
make_Pmat <- function(Q, V){
  cols = 1:(Q*V)
  rows_P1 = seq(1, (Q-1)*V+1, by=V)
  offset = rep(0:(V-1), each=Q)
  rows = rep(rows_P1, V) + offset
  P = t(sparseMatrix(i = rows, j = cols))
}


#' Construct block diagonal matrix of many small blocks
#'
#' @description Code for function provided in examples of bdiag function from Matrix package
#'
#' @param lmat List of matrices
#'
#' @import Matrix
#' @importFrom methods new
#'
#' @return A sparse matrix obtained by combining the arguments into a block diagonal matrix
#'
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}


#' Helper function for SQUAREM for estimating parameters
#'
#' @param theta_vec Vector of initial parameter values
#' @param template_mean Passed to UpdateTheta function
#' @param template_var Passed to UpdateTheta function
#' @param mesh Passed to UpdateTheta function
#' @param BOLD Passed to UpdateTheta function
#' @param C_diag Passed to UpdateTheta function
#' @param s0_vec Passed to UpdateTheta function
#' @param D Passed to UpdateTheta function
#' @param Dinv_s0 Passed to UpdateTheta function
#' @param common_smoothness Passed to UpdateTheta function
#' @param verbose Passed to UpdateTheta function
#'
#' @return Vector of updated parameter values
#'
UpdateThetaSQUAREM <- function(theta_vec, template_mean, template_var, mesh, BOLD, C_diag, s0_vec, D, Dinv_s0, common_smoothness, verbose){

  Q = nrow(template_mean)

  #convert theta vector to list format
  A <- theta_vec[1:(Q^2)]
  nu0_sq <- theta_vec[(Q^2)+1]
  kappa <- theta_vec[(Q^2+1)+(1:Q)]
  theta <- list(A = matrix(A, nrow=Q, ncol=Q),
                nu0_sq = nu0_sq,
                kappa = kappa)

  #update theta parameters
  if(verbose) cat('~~~~~~~~~~~ UPDATING PARAMETER ESTIMATES ~~~~~~~~~~~ \n')
  theta_new = UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=common_smoothness, verbose=verbose)

  #convert theta_new list to vector format
  theta_new$A <- as.matrix(theta_new$A)
  theta_new_vec <- unlist(theta_new[1:3]) #everything but LL
  names(theta_new_vec)[1] <- sum(theta_new$LL)
  return(theta_new_vec)
}


#' Helper function for SQUAREM for extracting log likelihood
#'
#' @param theta_vec Vector of current parameter values
#' @param template_mean Not used, but squarem will return error without
#' @param template_var  Not used, but squarem will return error without
#' @param mesh  Not used, but squarem will return error without
#' @param BOLD  Not used, but squarem will return error without
#' @param C_diag  Not used, but squarem will return error without
#' @param s0_vec  Not used, but squarem will return error without
#' @param D  Not used, but squarem will return error without
#' @param Dinv_s0  Not used, but squarem will return error without
#' @param common_smoothness  Not used, but squarem will return error without
#' @param verbose  Not used, but squarem will return error without
#'
#' @return Negative log-likelihood given current values of parameters
#'
LL_SQUAREM <- function(theta_vec, template_mean, template_var, mesh, BOLD, C_diag, s0_vec, D, Dinv_s0, common_smoothness, verbose){

  LL <- names(theta_vec)[1]
  print(LL)
  LL <- as.numeric(LL)
  print(LL)
  return(-1*LL)

}

