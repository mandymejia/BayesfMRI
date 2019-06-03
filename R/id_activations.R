#'
#'
#' @param object An object of class ‘"inla"’, a result of a call to inla
#' @param name Name of latent field on which to identify activations
#' @param mask Logical vector used to map beta estimates back to whole-brain field
#' @param mesh SPDE triangular mesh.  Only required if area.limit is specified.
#' @param session_names Names of sessions included in INLA model that resulted in object
#' @param threshold Activation threshold (e.g. 1%)
#' @param alpha Significance level (e.g. 0.05)
#' @param area.limit Below this value, activations will be considered spurious.  If NULL, no limit.
#'
#' @return Activation
#'
#' @export
#' @importFrom excursions excursions.inla excursions.inla.no.spurious
#'
#' @examples
id_activations <- function(object, name, mask, mesh=NULL, session_names, threshold, alpha, area.limit=NULL){

	nvox2 <- sum(mask)
	n_sess <- length(session_names)

	if(length(object$summary.random[[1]]$mean) != nvox2*n_sess) stop('Length of estimate vectors must equal number of sessions X number of voxels')

	mask2 <- mask
	mask2[mask==FALSE] <- NA #put NAs outside of the mask
	act <- rep(list(mask2), n_sess)
	names(act) <- session_names
	for(v in 1:n_sess){
		sess_name <- session_names[v]
		inds_v <- (1:nvox2) + (v-1)*nvox2 #indices of beta vector corresponding to session v
		if(is.null(area.limit)){
			res.exc <- excursions.inla(object, name=name, ind=inds_v, u=threshold, type='>', method='QC', alpha=alpha, F.limit=0.1)
		} else {
			res.exc <- excursions.inla.no.spurious(object, mesh=mesh, name=name, ind=inds_v, u=threshold, type='>', method='QC', alpha=alpha, area.limit = area.limit, use.continuous=FALSE, verbose=FALSE)
		}
		act_v <- res.exc$E[inds_v]
		act_v[is.na(act_v)] <- 0
		act[[v]][mask==TRUE] <- act_v
	}
	return(act)
}

