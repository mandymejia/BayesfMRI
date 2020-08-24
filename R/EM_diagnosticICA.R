#' EM Algorithm for Diagnostic ICA Models
#'
#' @param template_mean (A list of G matrices, each VxL) template mean estimates for each group 1 to G
#' @param template_var (A list of G matrices, each VxL) template variance estimates for each group 1 to G
#' @param BOLD  (VxL matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at mixing matrix A (LxL)
#' @param C_diag (Lx1) diagonal elements of residual covariance after dimension reduction
#' @param maxiter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#' @param verbose If TRUE, display progress of algorithm
#'
#' @export
#' @return  A list with 5 elements: theta (list of final parameter estimates), group_probs (posterior probabilities of group membership), subICmean (estimates of subject-level ICs), subICvar (variance of subject-level ICs), and success (flag indicating convergence (\code{TRUE}) or not (\code{FALSE}))
#'
EM_diagnosticICA = function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001, verbose=TRUE){

  #NOTE: This function expects the VxQ orientation for the templates and BOLD, since they are not transposed in diagnosticICA()
  # The corresponding function for templateICA expects the QxV orientation, since they ARE transposed in templateICA()

  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of brain locations
  L <- ncol(template_mean[[1]]) #number of ICs
  G <- length(template_mean)
  if(L > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(L > ntime) stop('Cannot estimate more ICs than time points.')

  iter = 1
  theta = theta0
  success = 1
  for(g in 1:G) {
    num_smallvar <- sum(template_var[[g]] < 1e-6)
    if(num_smallvar>0){
      if(verbose) cat(paste0('Setting ',num_smallvar,' (',round(100*num_smallvar/length(template_var[[g]]),1),'%) very small variance values in group ',g,' template to ',1e-6,'.\n'))
      template_var[[g]][template_var[[g]] < 1e-6] = 1e-6 #to prevent problems when inverting covariance
    }
  }

  err = 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))

    t00 <- Sys.time()
    theta_new = UpdateTheta.diagnosticICA(template_mean, template_var, BOLD, theta, C_diag, return_MAP=FALSE, verbose=verbose)
    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    A_old = theta$A
    A_new = theta_new$A
    #2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
    err = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")
    change = format(err, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change,' for A\n'))

    ### Move to next iteration
    theta <- theta_new
    iter = iter + 1
    if(iter > maxiter){
      success = 0;
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }

  MAP = UpdateTheta.diagnosticICA(template_mean, template_var, BOLD, theta, C_diag, return_MAP=TRUE, verbose=verbose)

  result <- list(group_probs = MAP$group_probs,
                 subjICmean=MAP$ICmean,
                 subjICvar=MAP$ICvar,
                 theta_MLE=theta,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean = template_mean,
                 template_var = template_var)

  return(result)
}


#' Parameter Estimates in EM Algorithm for Diagnostic ICA Model
#'
#' @param template_mean (A list of G matrices, each VxL) template mean estimates for each group 1 to G
#' @param template_var (A list of G matrices, each VxL) template variance estimates for each group 1 to G
#' @param BOLD  (VxL matrix) dimension-reduced fMRI data
#' @param theta (list) current parameter estimates
#' @param C_diag (Lx1) diagonal elements of residual covariance after dimension reduction
#' @param return_MAP If TRUE, returns the posterior mean and variance of the latent fields and group membership instead of the parameter estimates
#' @param verbose If TRUE, display progress of algorithm
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields
#'
UpdateTheta.diagnosticICA = function(template_mean, template_var, BOLD, theta, C_diag, return_MAP=FALSE, verbose=TRUE){

  L <- ncol(BOLD)
  nvox <- nrow(BOLD)
  G <- length(template_mean)

  #initialize new objects
  theta_new = list(A = matrix(NA, L, L))
  A_part1 = A_part2 = matrix(0, L, L) #two parts of product for A-hat (construct each looping over voxels)

  A <- theta$A
  C <- diag(C_diag)
  C_inv <- diag(1/(C_diag))
  At_Cinv <- t(A) %*% C_inv
  At_Cinv_A <- At_Cinv %*% A

  ##########################################
  ### POSTERIOR PROBABILITIES OF z
  ##########################################

  if(verbose) cat('Computing posterior probabilities of z (group membership) \n')

  exp_part <- rep(NA, G)
  pr_zy <- rep(NA, G)
  for(g in 1:G){

    #exp_part1 <- L*nvox*log(2*pi)
    exp_part1 <- 0 #irrelevant as long as noninformative prior on z=g
    exp_part2 <- sum(log(template_var[[g]])) #sum over v,ell

    exp_part3 <- 0
    for(v in 1:nvox){

      y_v = BOLD[v,]
      s0_gv = template_mean[[g]][v,]

      #mean and cov of y(v)|z
      mu_yz_v <- A %*% s0_gv
      Sigma_yz_v <- A %*% diag(template_var[[g]][v,]) %*% t(A) + C

      exp_part3_v <- t(y_v - mu_yz_v) %*% solve(Sigma_yz_v) %*% (y_v - mu_yz_v)
      exp_part3 <- exp_part3 + exp_part3_v

    }

    exp_part[g] <- -0.5 * (exp_part1 + exp_part2 + exp_part3[1])
  }

  for(g in 1:G){
    pr_zy_inv_g <- sum(exp(exp_part-exp_part[g])) # ( (e^M1 + e^M2 + e^M3) / e^M1 ) = e^(M1-M1) + e^(M2-M1) + e^(M3-M1) = 1 + e^(M2-M1) + e^(M3-M1)  <-- If any e^(Mk-Mg) Inf, the inverse will be zero so p(z=g|y)=0
    pr_zy[g] <- 1/pr_zy_inv_g
  }

  #fix numerical issues with very small values (values very close to 1 are numerically equal to 1, while values very close to zero are not)
  if(any(pr_zy==1)){
    pr_zy[pr_zy!=1] <- 0
  }

  if(verbose) print(paste(pr_zy, collapse=','))

  # if(is.infinite(exp(max(exp_part)-min(exp_part)))) {
  #   #this is for two groups, need to generalize for G>2
  #   pr_zy[which.max(exp_part)] <- 1
  #   pr_zy[which.min(exp_part)] <- 0
  # } else {
  #   exp_part <- exp_part - mean(exp_part) #this minimizes the most extreme magnitude in exp_part (to avoid )
  #   prod_pr_yz <- exp(exp_part)
  #   #compute p(z=g|y), g=1,...,G
  #   denom <- sum((1/G)*prod_pr_yz) #can change 1/G for unequal prior probs
  #   pr_zy <- (1/G)*prod_pr_yz / denom
  # }

  ##########################################
  ### POSTERIOR MOMENTS OF s
  ##########################################

  if(verbose) cat('Computing posterior moments of s (IC maps) \n')

  mu_sy <- array(NA, dim=c(nvox,L))
  mu_mut_sy <- array(NA, dim=c(nvox,L,L))
  var_sy <- array(NA, dim=c(nvox,L)) #only needed if return_MAP=TRUE
  for(v in 1:nvox){

    mu_sy_v <- rep(0, L)
    mu_mut_sy_v <- Sigma_sy_v <- array(0, dim=c(L,L))
    for(g in 1:G){

      D_gv_inv <- diag(1/template_var[[g]][v,])
      s0_gv <- template_mean[[g]][v,]
      Sigma_syz_gv <- solve(At_Cinv_A + D_gv_inv) #posterior covariance of s|z=g
      mu_syz_gv <- Sigma_syz_gv %*% (At_Cinv %*% BOLD[v,] + D_gv_inv %*% s0_gv) #posterior mean of s|z=g
      mu_mut_syz_gv <- mu_syz_gv %*% t(mu_syz_gv) + Sigma_syz_gv #posterior second moment of s|z=g

      mu_sy_v <- mu_sy_v + pr_zy[g]*mu_syz_gv #weighted sum over g=1..G --> posterior mean of s
      mu_mut_sy_v <- mu_mut_sy_v + pr_zy[g]*mu_mut_syz_gv
      Sigma_sy_v <- Sigma_sy_v + pr_zy[g]*Sigma_syz_gv #weighted sum over g=1..G --> posterior second moment of s

    }

    mu_sy[v,] <- mu_sy_v
    mu_mut_sy[v,,] <- mu_mut_sy_v
    var_sy[v,] <- diag(Sigma_sy_v) #only needed if return_MAP=TRUE

  }

  #check whether mu_mut_sy = mu_sy %*% t(mu_sy) + Sigma_sy

  if(return_MAP){
    result <- list(group_probs = pr_zy,
                   ICmean = mu_sy,
                   ICvar = var_sy)
    return(result)
  }

  ##########################################
  ### MLE OF A
  ##########################################

  if(verbose) cat('Updating A \n')

  A_part1 <- A_part2 <- matrix(0, L, L)
  for(v in 1:nvox){

    y_v = BOLD[v,]

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 = A_part1 + y_v %*% t(mu_sy[v,]) #LxL
    A_part2 = A_part2 + mu_mut_sy[v,,] #LxL

  }

  A_hat = (A_part1 %*% solve(A_part2))


  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  return(theta_new)
}


