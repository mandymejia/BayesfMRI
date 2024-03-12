#' GLM multi
#' 
#' @param y,X,X2 BOLD, design, nuisance
#' @param spatial Spatial info
#' @return Results for GLM multi
#' @keywords internal 
GLM_multi <- function(y, X, X2) {
  # Step 1: Identify no-signal locations. Compare null vs. canonical model using out-of-sample prediction error.

  nT <- nrow(y)

  #ingredients for prediction
  nT2 <- round(nT/2)
  inds1 <- 1:nT2
  inds2 <- setdiff(1:nT, inds1)
  y_list <- list(y[inds1,], y[inds2,])
  X2_list <- list(cbind(1, X2[inds1,]), cbind(1, X2[inds2,])) #for nuisance regression of the "held out" data
  which_can <- 22 #[TO DO] this must be provided as an argument
  X1_list <- list(X[inds1,,which_can], X[inds2,,which_can]) #canonical HRF task regressors

  RSS_OS <- matrix(0, nrow=nV_all, ncol=2)
  RSS_OS[mask==FALSE,] <- NA
  for(pred in 1:2){

    train <- pred
    test <- setdiff(1:2, train)

    # (i) estimate task coefficients based on training set (for canonical model only, not necessary for null model)
    y_train <- y_list[[train]]
    X2_train <- X2_list[[train]]
    X1_train <- X1_list[[train]]
    X_train_can <- cbind(X1_train, X2_train)
    print(head(X_train_can))
    save(X_train_can, file='~/Dropbox/RESEARCH/HRF-Adaptation/tmp.RData')
    XtX_inv <- try(Matrix::solve(Matrix::crossprod(X_train_can))) #this includes nuisance regressors
    if (inherits(XtX_inv, "try-error")) {
      warning(paste0("Numerical instability in design matrix"))
    }
    coefs_can <- (XtX_inv %*% t(X_train_can) %*% y_train)[1:nK,] #save task coefficients only (KxV)

    # (ii) do nuisance regression on test set
    y_test <- y_list[[test]]
    X2_test <- X2_list[[test]]
    X1_test <- X1_list[[test]]

    Imat <- diag(1, nrow(X2_test))
    Hmat <- X2_test %*% try(Matrix::solve(Matrix::crossprod(X2_test))) %*% t(X2_test)
    y_test_nuis <- (Imat - Hmat) %*% y_test #regress out nuisance from y
    X1_test_nuis <- (Imat - Hmat) %*% X1_test #regress out nuisance from task regressors

    # (iii) apply coefficients to generate prediction errors
    resid_list <- list(canonical = y_test_nuis - X1_test_nuis %*% coefs_can, null = y_test_nuis)
    RSS_OS_pred <- sapply(resid_list, function(x) colSums(x^2)) #Vx2
    RSS_OS[mask==TRUE,] <- RSS_OS[mask==TRUE,] + RSS_OS_pred #sum over both directions of prediction
  }

  noHRF <- (RSS_OS[,2] < RSS_OS[,1]) #for which locations is the null model RSS less than the canonical error RSS

  #loop over models
  nP <- dim(X)[3]
  RSS <- matrix(NA, nrow=nV_all, ncol=nP) #keep track of residual sum of squares (proxy for R^2 or AIC)

  if(verbose > 0) cat('\tFitting models: Model ')
  for(pp in 1:nP){

    if(verbose > 0) cat(paste0(pp,'\t'))

    X_sp <- cbind(X[,,pp], rep(1, nT), X2) #combine design with intercept and nuisance

    #1. Compute residual SD

    XtX_inv_pp <- try(Matrix::solve(Matrix::crossprod(X_sp)))
    if (inherits(XtX_inv_pp, "try-error")) {
      warning(paste0("Numerical instability in design matrix for model ",pp))
    }
    coef_pp <- XtX_inv_pp %*% t(X_sp) %*% y
    resid_pp <- y - X_sp %*% coef_pp #TxV matrix
    RSS[mask==TRUE,pp] <- sqrt(colSums(resid_pp^2)/(nT - ncol(X_sp)))

    #determine best model (minimum residual squared error)
    bestmodel <- apply(RSS, 1, function(x){
      wm <- which.min(x)
      varx <- var(x, na.rm=TRUE)
      if(is.na(varx)) varx <- 0
      if(varx==0) wm <- NA
      wm
    })

    result_multiple <- list(bestmodel = bestmodel,
                                  noHRF = noHRF,
                                  RSS = RSS,
                                  RSS_OS = RSS_OS)
  }
  result_multiple
}