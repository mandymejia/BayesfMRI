#y is the TxV data matrix containing the fMRI timeseries
#X is the TxK design matrix with K task-related columns 

organize_data <- function(y, X){
	
	ntime <- nrow(y)
	nvox <- ncol(y)

	#check orientation, send warning message and transpose if necessary
	if(ntime > nvox){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		y <- t(y)
		ntime <- nrow(y)
		nvox <- ncol(y)
	}

	y <- as.vector(y) #makes a vector (y_1,...,y_V), where y_v is the timeseries for data location v
	ix <- 1:(ntime*nvox)
	iy <- rep(1:nvox, each = ntime)

	K <- ncol(X)
	for(k in 1:K){
		X_k <- sparseMatrix(ix, iy, x=rep(X[,k], nvox)) 
		if(k==1) A <- X_k else A <- cbind(A, X_k)
	}

	result <- list(y=y, A=A)
	return(result)
}