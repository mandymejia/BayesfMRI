#' GLM multi
#'
#' @param y,X,X2 BOLD, design, nuisance
#' @param Xc (Optional) canonical design matrix
#' @param verbose verbose?
#' @return Results for GLM multi
#' @importFrom stats as.formula var pf pchisq
#' @keywords internal
GLM_multi <- function(y, X, X2, Xc=NULL, verbose=TRUE) {
  # Step 1: Identify no-signal locations. Compare null vs. canonical model using out-of-sample prediction error.
  nT <- nrow(y)
  nV_D <- ncol(y)
  nK <- ncol(X)

  # F = (RSS0 - RSS1)/RSS1 * (n - p1)/(p1 - p0), where
  #   model 0 = null model & model 1 = canonical HRF model.
  #   So p1 - p0 = nK (# of task regressors).

  # For each voxel:
  # a) Fit canonical HRF model (model 1), record RSS1
  # b) Fit null model (model 0), record RSS0
  # c) Compute F-statistic
  # d) Compute p-value, since F ~ F(nK, n-p1) ( Basically a chi-sq, since n -> Inf.  nK*F ~ X2(nK))

  Fstat <- pvalF <- rep(0, nV_D)

  if (!is.null(Xc)) {
    X1 <- cbind(Xc, rep(1, nT), X2) #combine design with intercept and nuisance
    X0 <- cbind(rep(1, nT), X2) #no task regressors in null model

    #compute RSS (steps (a) and (b))
    XtX_inv1 <- try(Matrix::solve(Matrix::crossprod(X1)))
    XtX_inv0 <- try(Matrix::solve(Matrix::crossprod(X0)))
    if (inherits(XtX_inv1, "try-error")) { warning(paste0("Numerical instability in design matrix for canonical model ")) }
    if (inherits(XtX_inv0, "try-error")) { warning(paste0("Numerical instability in design matrix for null model ")) }
    coef1 <- XtX_inv1 %*% t(X1) %*% y
    coef0 <- XtX_inv0 %*% t(X0) %*% y
    resid1 <- y - X1 %*% coef1 #TxV matrix
    resid0 <- y - X0 %*% coef0 #TxV matrix
    RSS1 <- colSums(resid1^2)
    RSS0 <- colSums(resid0^2)

    #compute F-statistic (step c)
    DOF1 <- nT - ncol(X1) #will be over-estimated, but shouldn't matter much because it's essentially Inf
    Fstat <- (RSS0 - RSS1)/RSS1 * DOF1/nK
    pvalF <- 1 - pf(Fstat, df1 = nK, df2 = DOF1)
    #pval_Chi2 <- 1 - pchisq(nK*Fstat, df = nK)
  }

  # #ingredients for prediction
  # nT2 <- round(nT/2)
  # inds1 <- 1:nT2
  # inds2 <- setdiff(1:nT, inds1)
  # y_list <- list(y[inds1,], y[inds2,])
  # X2_list <- list(cbind(1, X2[inds1,]), cbind(1, X2[inds2,])) #for nuisance regression of the "held out" data
  # X1_list <- list(design_can[inds1,], design_can[inds2,]) #canonical HRF task regressors
  #
  # RSS_OS <- matrix(0, nrow=nV_D, ncol=2)
  # RSS_OS[mask==FALSE,] <- NA
  # for(pred in 1:2){
  #
  #   train <- pred
  #   test <- setdiff(1:2, train)
  #
  #   # (i) estimate task coefficients based on training set (for canonical model only, not necessary for null model)
  #   y_train <- y_list[[train]]
  #   X2_train <- X2_list[[train]]
  #   X1_train <- X1_list[[train]]
  #   X_train_can <- cbind(X1_train, X2_train)
  #   nK <- ncol(X1_train)
  #
  #   XtX_inv <- try(Matrix::solve(Matrix::crossprod(X_train_can))) #this includes nuisance regressors
  #   if (inherits(XtX_inv, "try-error")) {
  #     warning(paste0("Numerical instability in design matrix"))
  #   }
  #   coefs_can <- (XtX_inv %*% t(X_train_can) %*% y_train)[1:nK,] #save task coefficients only (KxV)
  #
  #   # (ii) do nuisance regression on test set
  #   y_test <- y_list[[test]]
  #   X2_test <- X2_list[[test]]
  #   X1_test <- X1_list[[test]]
  #
  #   Imat <- diag(1, nrow(X2_test))
  #   Hmat <- X2_test %*% try(Matrix::solve(Matrix::crossprod(X2_test))) %*% t(X2_test)
  #   y_test_nuis <- (Imat - Hmat) %*% y_test #regress out nuisance from y
  #   X1_test_nuis <- (Imat - Hmat) %*% X1_test #regress out nuisance from task regressors
  #
  #   # (iii) apply coefficients to generate prediction errors
  #   resid_list <- list(canonical = y_test_nuis - X1_test_nuis %*% coefs_can, null = y_test_nuis)
  #   RSS_OS_pred <- sapply(resid_list, function(x) colSums(x^2)) #Vx2
  #   RSS_OS <- RSS_OS + RSS_OS_pred #sum over both directions of prediction
  # }
  #noHRF <- (RSS_OS[,2] < RSS_OS[,1]) #for which locations is the null model RSS less than the canonical error RSS

  #loop over models
  nP <- dim(X)[3]
  RSS <- matrix(NA, nrow=nV_D, ncol=nP) #keep track of residual sum of squares (proxy for R^2 or AIC)

  if(verbose > 0) { cat('\tFitting models: Model ') }
  for(pp in 1:nP){
    if(verbose > 0) { cat(paste0(pp,'\t')) }

    X_sp <- cbind(X[,,pp], rep(1, nT), X2) #combine design with intercept and nuisance

    #1. Compute residual SD

    XtX_inv_pp <- try(Matrix::solve(Matrix::crossprod(X_sp)))
    if (inherits(XtX_inv_pp, "try-error")) {
      warning(paste0("Numerical instability in design matrix for model ",pp))
    }
    coef_pp <- XtX_inv_pp %*% t(X_sp) %*% y
    resid_pp <- y - X_sp %*% coef_pp #TxV matrix
    RSS[,pp] <- sqrt(colSums(resid_pp^2)/(nT - ncol(X_sp)))
  }
  if (verbose > 0) { cat("\n") }

  #determine best model (minimum residual squared error)
  bestmodel <- apply(RSS, 1, function(x){
    wm <- which.min(x)
    varx <- var(x, na.rm=TRUE)
    if(is.na(varx)) varx <- 0
    if(varx==0) wm <- NA
    wm
  })

  list(
    bestmodel = bestmodel,
    Fstat = Fstat,
    pvalF = pvalF
  )
}
