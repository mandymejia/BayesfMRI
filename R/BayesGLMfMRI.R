#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh.
#' @param mesh A `inla.mesh` object.  Must be provided if and only if `vertices` and `faces` are not.
#' @param mask A vector of 0s and 1s of length V, where locations with value 0 will be excluded from analysis.
#' @param scale If TRUE, scale timeseries data so estimates represent percent signal change.  Else, do not scale.
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @examples \dontrun{}
BayesGLMfMRI <- function(data, vertices, faces, mesh, mask, scale=TRUE){

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V0 = full number of data locations
  # V = masked number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  #need to check that sessions are consistent in terms of V, K?

  #INLA:::inla.dynload.workaround() #avoid error on creating mesh

  #check that only mesh OR vertices+faces supplied
  has_mesh <- !missing(mesh)
  has_verts_faces <- !missing(vertices) & !missing(faces)
  has_howmany <- has_mesh + has_verts_faces
  if(has_howmany != 1) stop('Must supply EITHER mesh OR vertices and faces.')

  #maybe also allow the user to supply a mesh object
  #if mesh and mask both supplied, will need to deal with that...

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists, but it is not')

  V <- ncol(data[[1]]$BOLD)
  K <- ncol(data[[1]]$design)
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  if(missing(mask)) mask <- rep(1, V)
  if(missing(mesh)) mesh <- make_mesh(vertices, faces, mask)

  spde <- inla.spde2.matern(mesh)
  #areas <- compute_vertex_areas(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL

  for(s in 1:n_sess){

      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD
      BOLD_s <- BOLD_s[,mask==1]

      #scale data to represent % signal change
      BOLD_s <- scale_timeseries(t(BOLD_s))

      #regress nuisance parameters from BOLD data and design matrix
      if('nuisance' %in% names(data[[s]])){
        design_s <- data[[s]]$design
        nuisance_s <- data[[s]]$nuisance
        y_reg <- nuisance_regress(BOLD_s, nuisance_s)
        X_reg <- nuisance_regress(design_s, nuisance_s)
      } else {
        y_reg <- BOLD_s
        X_reg <- data[[s]]$design
      }

      #set up data and design matrix
      data_org <- organize_data(y_reg, X_reg)
      y_vec <- data_org$y
      X_list <- list(data_org$A)
      names(X_list) <- session_names[s]

      y_all <- c(y_all, y_vec)
      X_all_list <- c(X_all_list, X_list)
  }

  #construct betas and repls objects
  replicates_list <- organize_replicates(n_sess=n_sess, n_task=K, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  formula <- make_formula(beta_names = names(betas), repl_names = names(repls), model_name = 'spde', hyper_initial = c(-2,2))
  model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

  #estimate model using INLA
  INLA_result <- estimate_model(formula=formula, data=model_data, A=model_data$X, prec_initial=1)

  #extract useful stuff from INLA model result
  beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names) #posterior means of latent task field
  theta_posteriors <- get_posterior_densities(object=INLA_result, spde) #hyperparameter posterior densities

  #identify areas of activation if activation threshold(s) specified by user

  #construct object to be returned
  result <- list(model=INLA_result, mesh=mesh, sessions=session_names, beta_estimates=beta_estimates, theta_posteriors=theta_posteriors)
  return(result)


#   # ID AREAS OF ACTIVATION
#
#   binarize <- function(x, p){ return(x > p) }
#
#   #6 hours with cluster-wise, 2 hours with voxel-wise only
#   t0 <- Sys.time()
#   for(u in 1:U){
#
#     thr_u <- thresholds[u]
#     print(paste0('Threshold: ', thr_u))
#
#     #voxel-wise joint PPM
#     time_u <- system.time(excur_u <- id_activations(object=result, name='bbeta1', mask=mask2_sh, session_names=session_names, threshold=thr_u, alpha=0.05))
#     print(time_u)
#     row_htu <- (comptime_excur_s$hemisphere == h) & (comptime_excur_s$task_type == task_type) & (comptime_excur_s$threshold == thr_u)
#     comptime_excur_s$comptime[row_htu] <- time_u[3] #elapsed time (sec)
#     save(excur_u, file=paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'.Rdata'), sep='_'))
#
#     for(a in c(95,99)){
#       print(paste0('Significance: ', a))
#       #voxel-wise joint PPM
#       excur_ua <- lapply(excur_u, binarize, p=a/100)
#
#       # #cluster-wise joint PPM
#       # time_ua <- system.time(excur_clust_ua <- id_activations(object=result, name='bbeta1', mask=mask2_sh, mesh=mesh_sh, session_names=session_names, threshold=thr_u, alpha=(100-a)/100, area.limit=10))
#       # row_htua <- (comptime_excur_clust_s$hemisphere == h) & (comptime_excur_clust_s$task_type == task_type) & (comptime_excur_clust_s$threshold == thr_u) & (comptime_excur_clust_s$alpha == a)
#       # comptime_excur_clust_s$comptime[row_htua] <- time_ua[3] #elapsed time (sec)
#       # excur_clust_ua <- lapply(excur_clust_ua, binarize, p=a/100)
#
#       #save activation sets for each session
#       for(v in 1:n_sess){
#         print(v)
#         sess_name <- session_names[v]
#         fname_v <- paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'_',a,'.csv'), sep='_')
#         write.csv(excur_ua[[v]], fname_v, row.names=FALSE)
#         # fname_v <- paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'_',a,'_clust10.csv'), sep='_')
#         # write.csv(excur_clust_ua[[v]], fname_v, row.names=FALSE)
#
#         #compute size of active area
#         inds_v <- which(excur_ua[[v]])
#         area_v <- sum(areas_full_sh[inds_v])
#         row_v <- which(paste(active_areas_sht$visit, active_areas_sht$task, sep='_')==sess_name & active_areas_sht$threshold==thr_u & active_areas_sht$significance==a)
#         active_areas_sht$area[row_v] <- area_v
#       }
#     }
#   }
#   Sys.time() - t0
#
#   active_areas_s <- rbind(active_areas_s, active_areas_sht)
#   fname <- paste0('../../Results/',s,'_activeareas.Rdata')
#   save(active_areas_s, file=fname)
#
}
