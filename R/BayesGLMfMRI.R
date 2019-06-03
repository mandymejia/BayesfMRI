#' Applies spatial Bayesian GLM to task fMRI data
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param vertices A Vx3 matrix of vertex locations of the triangular mesh in Euclidean space.
#' @param faces A Wx3 matrix, where each row contains the vertex indices for a given face or triangle in the triangular mesh.
#' @param mask A vector of 0s and 1s of length V, where locations with value 0 will be excluded from analysis.
#'
#' @return A list containing...
#' @export
#'
#' @examples
BayesGLMfMRI <- function(data, vertices, faces, mask=NULL){

  #check whether data is a list
  #check whether each element of data is a session (use is.session)
  # V0 = full number of data locations
  # V = masked number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  INLA:::inla.dynload.workaround() #avoid error on creating mesh

  mesh <- make_mesh(vertices, faces, mask)
  spde <- inla.mesh.create(mesh)
  areas <- compute_vertex_areas(mesh)

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL

  session_names <- names(data)
  n_sess <- length(session_names)
  for(s in 1:n_sess){

      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD
      BOLD_s <- BOLD_s[,mask==1]

      #scale data to represent % signal change
      BOLD_s <- scale_timeseries(BOLD_s)

      #regress nuisance parameters from BOLD data and design matrix
      design_s <- data[[s]]$design
      nuisance_s <- data[[s]]$nuisance
      y_reg <- nuisance_regress(BOLD_s, nuisance_s)
      X_reg <- nuisance_regress(design_s, nuisance_s)

      #set up data and design matrix
      data_org <- organize_data(y_reg, X_reg)
      y_vec <- data_org$y
      X_list <- list(data_org$A)
      names(X_list) <- session_names[s]

      y_all <- c(y_all, y_vec)
      X_all_list <- c(X_all_list, X_list)
  }

  #construct betas and repls objects
  replicates_list <- organize_replicates(n_sess=n_sess, nx=nx, mesh=mesh)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  formula <- make_formula(beta_names = names(betas), repl_names = names(repls), model_name = 'spde', hyper_initial = c(-2,2))
  model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

  #estimate model using INLA
  result <- estimate_model(formula=formula, data=model_data, A=model_data$X, prec_initial=1)

  row_ht <- (comptime_s$hemisphere == h) & (comptime_s$task_type == task_type)
  comptime_s$comptime[row_ht] <- mod_time
  comptime_s$units[row_ht] <- attributes(mod_time)$units

  # EXTRACT AND SAVE ESTIMATES

  session_names <- names(X_all_list)
  betas_sh <- extract_estimates(object=result, mask=mask2_sh, session_names=session_names)
  for(v in 1:n_sess){
    sess_name <- session_names[v]
    fname_v <- paste(paste0('../../Results/beta_estimates/',s), sess_name, paste0(h,'.csv'), sep='_')
    write.csv(betas_sh[[v]], fname_v, row.names=FALSE)
  }

  # ID AREAS OF ACTIVATION

  binarize <- function(x, p){ return(x > p) }

  #6 hours with cluster-wise, 2 hours with voxel-wise only
  t0 <- Sys.time()
  for(u in 1:U){

    thr_u <- thresholds[u]
    print(paste0('Threshold: ', thr_u))

    #voxel-wise joint PPM
    time_u <- system.time(excur_u <- id_activations(object=result, name='bbeta1', mask=mask2_sh, session_names=session_names, threshold=thr_u, alpha=0.05))
    print(time_u)
    row_htu <- (comptime_excur_s$hemisphere == h) & (comptime_excur_s$task_type == task_type) & (comptime_excur_s$threshold == thr_u)
    comptime_excur_s$comptime[row_htu] <- time_u[3] #elapsed time (sec)
    save(excur_u, file=paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'.Rdata'), sep='_'))

    for(a in c(95,99)){
      print(paste0('Significance: ', a))
      #voxel-wise joint PPM
      excur_ua <- lapply(excur_u, binarize, p=a/100)

      # #cluster-wise joint PPM
      # time_ua <- system.time(excur_clust_ua <- id_activations(object=result, name='bbeta1', mask=mask2_sh, mesh=mesh_sh, session_names=session_names, threshold=thr_u, alpha=(100-a)/100, area.limit=10))
      # row_htua <- (comptime_excur_clust_s$hemisphere == h) & (comptime_excur_clust_s$task_type == task_type) & (comptime_excur_clust_s$threshold == thr_u) & (comptime_excur_clust_s$alpha == a)
      # comptime_excur_clust_s$comptime[row_htua] <- time_ua[3] #elapsed time (sec)
      # excur_clust_ua <- lapply(excur_clust_ua, binarize, p=a/100)

      #save activation sets for each session
      for(v in 1:n_sess){
        print(v)
        sess_name <- session_names[v]
        fname_v <- paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'_',a,'.csv'), sep='_')
        write.csv(excur_ua[[v]], fname_v, row.names=FALSE)
        # fname_v <- paste(paste0('../../Results/activations/',s), sess_name, paste0(h,'_thr',thr_u,'_',a,'_clust10.csv'), sep='_')
        # write.csv(excur_clust_ua[[v]], fname_v, row.names=FALSE)

        #compute size of active area
        inds_v <- which(excur_ua[[v]])
        area_v <- sum(areas_full_sh[inds_v])
        row_v <- which(paste(active_areas_sht$visit, active_areas_sht$task, sep='_')==sess_name & active_areas_sht$threshold==thr_u & active_areas_sht$significance==a)
        active_areas_sht$area[row_v] <- area_v
      }
    }
  }
  Sys.time() - t0

  active_areas_s <- rbind(active_areas_s, active_areas_sht)
  fname <- paste0('../../Results/',s,'_activeareas.Rdata')
  save(active_areas_s, file=fname)

  ### EXTRACT MARGINAL POSTERIORS OF HYPERPARAMETERS

  theta_posteriors_sht <- get_posterior_densities(result, spde_sh, names(betas))



}
