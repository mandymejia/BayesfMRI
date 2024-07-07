checkX_cor <- function(X){
  suppressWarnings(corX <- cor(X))
  diag(corX) <- NA
  corX
}

checkX_VIF <- function(X){
  nT <- nrow(X)
  nK <- ncol(X)
  
  y <- rnorm(nT) #add fake y variable, has no influence
  Xnames <- paste0("X", seq(nK))
  df <- as.data.frame(cbind(X, y)); names(df) <- c(Xnames,"y")
  f <- as.formula(paste0('y ~ ',paste(Xnames, collapse = " + ")))
  suppressWarnings(try(car::vif(try(lm(f, data = df)))))
}
