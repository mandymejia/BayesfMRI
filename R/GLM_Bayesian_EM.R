      ## EM Model. ---------------------------------------------------------------
      #if(do_EM) {
        #stop()
        # if (!requireNamespace("MatrixModels", quietly = TRUE)) {
        #   stop("EM requires the `MatrixModels` package. Please install it.", call. = FALSE)
        # }
        # if (verbose>0) cat('\tEstimating Bayesian model with EM.\n')
        # Psi_k <- spde$Amat
        # Psi <- Matrix::bdiag(rep(list(Psi_k),nK))
        # A <- Matrix::crossprod(model_data$X %*% Psi)
        # # Initial values for kappa and tau
        # kappa2 <- 4
        # phi <- 1 / (4*pi*kappa2*4)
        # # Using values based on the classical GLM
        # if (verbose>0) cat("\t\tFinding best guess initial values.\n")
        # beta_hat <- MatrixModels:::lm.fit.sparse(model_data$X, model_data$y)
        # res_y <- (model_data$y - model_data$X %*% beta_hat)@x
        # sigma2 <- stats::var(res_y)
        # beta_hat <- matrix(beta_hat, ncol = nK*nS)
        # rcpp_spde <- create_listRcpp(spde$spde)
        # if(nS > 1) {
        #   field_cols <- sapply(seq(nS), function(ss) seq(nK) + nK *(ss - 1))
        #   beta_hat <- apply(field_cols,1,function(x) beta_hat[,x])
        # }
        # n_threads <- parallel::detectCores()
        # n_threads <- min(n_threads,nK,n_threads)
        # cl <- parallel::makeCluster(n_threads)
        # kappa2_phi_rcpp <- parallel::parApply(
        #   cl = cl,
        #   beta_hat,
        #   2,
        #   .initialKP,
        #   theta = c(kappa2, phi),
        #   spde = rcpp_spde,
        #   n_sess = nS,
        #   tol = emTol,
        #   verbose = FALSE
        # )
        # parallel::stopCluster(cl)
        # if (verbose>0) cat("\t\tDone!\n")
        # theta <- c(t(kappa2_phi_rcpp), sigma2)
        # theta_init <- theta
        # Ns <- 50 # This is a level of approximation used for the Hutchinson trace estimator
        # if(verbose>0) cat("\t\tStarting EM algorithm.\n")
        # em_output <-
        #   .findTheta(
        #     theta = theta,
        #     spde = rcpp_spde,
        #     y = model_data$y,
        #     X = model_data$X,
        #     QK = make_Q(theta, rcpp_spde, nS),
        #     Psi = as(Psi, "dgCMatrix"),
        #     A = as(A, "dgCMatrix"),
        #     Ns = 50,
        #     tol = emTol,
        #     verbose = verbose>0
        #   )
        # if(verbose>0) cat("\t\tEM algorithm complete!\n")
        # kappa2_new <- phi_new <- sigma2_new <- mu_theta <- NULL
        # list2env(em_output, envir = environment())
        # Qk_new <- mapply(spde_Q_phi,kappa2 = kappa2_new, phi = phi_new,
        #                  MoreArgs = list(spde=rcpp_spde), SIMPLIFY = F)
        # Q_theta <- Matrix::bdiag(Qk_new)
        # if(nS > 1) Q_theta <- Matrix::bdiag(lapply(seq(nS), function(x) Q_theta))
        # Sig_inv <- Q_theta + A/sigma2_new
        # m <- Matrix::t(model_data$X%*%Psi)%*%model_data$y / sigma2_new
        # mu_theta <- Matrix::solve(Sig_inv, m)
        # # Prepare results
        # field_estimates <- matrix(NA, nrow = length(mask), ncol = nK*nS)
        # field_estimates[mask == 1,] <- matrix(mu_theta,nrow = nV, ncol = nK*nS)
        # colnames(field_estimates) <- rep(field_names, nS)
        # field_estimates <- lapply(seq(nS), function(ss) field_estimates[,(seq(nK) + nK * (ss - 1))])
        # names(field_estimates) <- session_names
        # avg_field_estimates <- NULL
        # if(combine_sessions) avg_field_estimates <- Reduce(`+`,field_estimates) / nS
        # theta_estimates <- c(sigma2_new,c(phi_new,kappa2_new))
        # names(theta_estimates) <- c("sigma2",paste0("phi_",seq(nK)),paste0("kappa2_",seq(nK)))
        # #extract stuff needed for group analysis
        # tau2_init <- 1 / (4*pi*theta_init[seq(nK)]*theta_init[(seq(nK) + nK)])
        # mu_init <- c(log(1/tail(theta_init,1)), c(rbind(log(sqrt(tau2_init)),log(sqrt(theta_init[seq(nK)])))))
        # tau2 <- 1 / (4*pi*kappa2_new*phi_new)
        # mu_theta <- c(log(1/sigma2_new),c(rbind(log(sqrt(tau2)),log(sqrt(kappa2_new)))))
        # if (verbose>0) cat("\t\tDone!\n")