
    # FIR Model ------------------------------------------------------------------
    # result_FIR <- vector('list', length=nS)
    # for (ss in seq(nS)) {
    #
    #   #check whether to proceed with FIR modeling
    #   FIR_ss <- data[[ss]]$design_FIR
    #   if(is.null(FIR_ss)) next()
    #   nFIR <- ncol(FIR_ss)
    #   print(paste0('Number of FIR regressors: ', nFIR))
    #   print(paste0('Number of volumes: ', nrow(FIR_ss)))
    #   if(nFIR > nrow(FIR_ss)){
    #     warning('More FIR regressors than volumes. Consider reducing FIR_nsec.')
    #     next()
    #   }
    #
    #   if (verbose>0) cat("\tFitting FIR model.\n")
    #
    #   #set up vectorized data and big sparse design matrix
    #   #if(!do_pw) data_ss <- sparse_and_PW(data[[ss]]$BOLD, data[[ss]]$design)
    #   #if(do_pw) data_ss <- data[[ss]] #data has already been "organized" (big sparse design) in prewhitening step above
    #
    #   #y_ss <- matrix(data_ss$BOLD, nrow=nT) #a vector (grouped by location)
    #   y_ss <- data[[ss]]$BOLD #[TO DO] implement prewhitening case (may not be needed if we are not doing inference)
    #   X_ss <- cbind(FIR_ss, 1, data[[ss]]$nuisance) #need the intercept since FIR bases are not centered
    #
    #   #fit model
    #   beta_hat_s <- matrix(NA, nV_all, nFIR)
    #   XTX_inv <- try(Matrix::solve(Matrix::crossprod(X_ss)))
    #   if (inherits(XTX_inv, "try-error")) {
    #     stop("There is some numerical instability in your design matrix (due to very large or very small values). Scaling the design matrix is suggested.")
    #   }
    #   coef_s <- as.matrix(XTX_inv %*% t(X_ss) %*% y_ss) #a vector of (estimates for location 1, estimates for location 2, ...)
    #   beta_hat_s[mask==TRUE,] <- t(coef_s[1:nFIR,]) #drop the intercept and nuisance, transpose to V x nFIR
    #
    #   colnames(beta_hat_s) <- colnames(FIR_ss)
    #   result_FIR[[ss]] <- beta_hat_s
    # }
    # names(result_FIR) <- session_names
    #