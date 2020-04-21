#' Applies spatial Bayesian GLM to task fMRI data for 3D subcortical volumes
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param locations Vx3 matrix of x,y,z coordinates of each voxel
#' @param labels Vector of length V of region labels
#' @param groups_df Data frame indicating the name and model group of each region.  See Details.
#' @param scale If TRUE, scale timeseries data so estimates represent percent signal change.  If FALSE, just center the data and design to exclude the baseline field.
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is TRUE. Required for running \code{id_activations} on \code{BayesGLM} model object.
#' @param outfile File name where results will be written (for use by \code{BayesGLM_grp}).
#' @param GLM If TRUE, classical GLM estimates will also be returned
#' @param num.threads Maximum number of threads the inla-program will use for model estimation
#' @param verbose Boolean indicating if the inla-program should run in a verbose mode (default FALSE).
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @details The subcortical data is separated into regions, whose sizes range from approximately 100 voxels to approximately 9000 voxels.  Smaller regions are grouped together to improve model fit.
#' The \code{groups_df} argument specifies which regions are grouped together.  This argument should be a data frame with R rows (the number of regions) and three columns: label, region, and group.
#' The label column is the numerical identifier of each region; the region column contains the region names, and the group column contains the model group assignments (e.g. 1,2,3). Regions to be excluded
#' from analysis are indicated by NA in the group assignment.
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#' @examples \dontrun{}
BayesGLM_vol3D <- function(data, locations, labels, groups_df, scale=TRUE, return_INLA_result=FALSE, outfile = NULL, GLM = TRUE, num.threads = 6, verbose=FALSE){

  # Check to see that the INLA package is installed
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("This function requires the INLA package (see www.r-inla.org/download)")


  # Check to see if PARDISO is installed
  if(!exists("inla.pardiso.check", mode = "function")){
    warning("Please update to the latest version of INLA for full functionality and PARDISO compatibility (see www.r-inla.org/download)")
  }else{
    if(inla.pardiso.check() == "FAILURE: PARDISO IS NOT INSTALLED OR NOT WORKING"){
      warning("Consider enabling PARDISO for faster computation (see inla.pardiso())")}
    #inla.pardiso()
  }

  #check that all elements of the data list are valid sessions and have the same number of locations and tasks
  session_names <- names(data)
  n_sess <- length(session_names)

  if(!is.list(data)) stop('I expect data to be a list, but it is not')
  data_classes <- sapply(data, 'class')
  if(! all.equal(unique(data_classes),'list')) stop('I expect data to be a list of lists (sessions), but it is not')

  V <- ncol(data[[1]]$BOLD) #number of data locations
  K <- ncol(data[[1]]$design) #number of tasks
  for(s in 1:n_sess){
    if(! is.session(data[[s]])) stop('I expect each element of data to be a session object, but at least one is not (see `is.session`).')
    if(ncol(data[[s]]$BOLD) != V) stop('All sessions must have the same number of data locations, but they do not.')
    if(ncol(data[[s]]$design) != K) stop('All sessions must have the same number of tasks (columns of the design matrix), but they do not.')
  }

  if(is.null(outfile)){
    warning('No value supplied for outfile, which is required for group modeling (see help(BayesGLM_group)).')
  }

  if(nrow(locations) != V) stop('The locations argument should have V rows, the number of data locations in BOLD.')
  if(ncol(locations) != 3) stop('The locations argument should have 3 columns indicating the x,y,z coordinate of each data location.')

  ### Fit classical GLM

  if(GLM){
    # estimate model using GLM
    GLM_result <- classicalGLM(data)
  } else {
    GLM_result <- NULL
  }

  ### Run SPDE object for each group of regions

  regions_set <- unique(labels)
  num_regions <- length(regions_set)
  if(class(groups_df) != 'data.frame') stop('The groups_df argument should be a data frame with the following columns: label, region, group.')
  if(!all.equal(sort(names(groups_df)), sort(c('label','region','group')))) stop('The groups_df argument should be a data frame with the following columns: label, region, group.')
  if(nrow(groups_df) != num_regions) stop('The groups_df argument should contain one row for each region represented in the labels vector.')


  ### Fit model for each group of regions

  groups_df <- groups_df[!is.na(groups_df$group),]
  group_set <- unique(groups_df$group)

  beta_estimates_all <- matrix(NA, nrow=V, ncol=K)
  beta_estimates_all <- rep(list(beta_estimates_all), n_sess)
  INLA_result_all <- vector('list', length=length(group_set))
  spde_all <- vector('list', length=length(group_set))
  names(INLA_result_all) <- paste0('model_group_',group_set)
  names(spde_all) <- paste0('model_group_',group_set)

  for(grp in group_set){
    label_set_grp <- groups_df$label[groups_df$group==grp]
    name_set_grp <- groups_df$region[groups_df$group==grp]
    inds_grp <- labels %in% label_set_grp
    locs_grp <- locations[inds_grp,]
    labels_grp <- labels[inds_grp]

    paste0('Estimating Model 1 (', paste(name_set_grp, collapse = ', '), ')')

    spde_grp <- create_spde_vol3D(locs=locs_grp, labs=labels_grp, lab_set=label_set_grp)
    #plot(spde_grp)
    spde <- spde_grp$spde

    #collect data and design matrices
    y_all <- c()
    X_all_list <- NULL

    for(s in 1:n_sess){

      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD[,inds_grp]

      #scale data to represent % signal change (or just center if scale=FALSE)
      BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale, transpose = FALSE)
      design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline

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
      data_org <- organize_data(y_reg, X_reg, transpose = FALSE)
      y_vec <- data_org$y
      X_list <- list(data_org$A)
      names(X_list) <- session_names[s]

      y_all <- c(y_all, y_vec)
      X_all_list <- c(X_all_list, X_list)
    }


    #construct betas and repls objects
    replicates_list <- organize_replicates(n_sess=n_sess, n_task=K, mesh = spde_grp)
    betas <- replicates_list$betas
    repls <- replicates_list$repls

    beta_names <- names(betas)
    repl_names <- names(repls)
    n_beta <- length(names(betas))
    hyper_initial <- c(-2,2)
    hyper_initial <- rep(list(hyper_initial), n_beta)
    hyper_vec <- paste0(', hyper=list(theta=list(initial=', hyper_initial, '))')

    formula_vec <- paste0('f(',beta_names, ', model = spde, replicate = ', repl_names, hyper_vec, ')')
    formula_vec <- c('y ~ -1', formula_vec)
    formula_str <- paste(formula_vec, collapse=' + ')
    formula <- as.formula(formula_str, env = globalenv())

    model_data <- make_data_list(y=y_all, X=X_all_list, betas=betas, repls=repls)

    save(formula, model_data, spde, file='sample_data_inla5.Rdata')

    # estimate model using INLA
    t0 <- Sys.time()
    INLA_result_grp <- estimate_model(formula=formula, data=model_data, A=model_data$X, spde=spde, prec_initial=1, num.threads = num.threads, verbose=verbose)
    print(Sys.time() - t0)

    #extract beta estimates and project back to data locations for current group
    beta_estimates_grp <- extract_estimates(object=INLA_result_grp, session_names=session_names) #posterior means of latent task field
    for(s in 1:n_sess){
      beta_estimates_all[[s]][inds_grp,] <- as.matrix(spde_grp$Amat %*% beta_estimates_grp[[s]])
    }

    #extract theta posteriors
    theta_posteriors_grp <- get_posterior_densities_vol3D(object=INLA_result_grp, spde) #hyperparameter posterior densities
    theta_posteriors_grp$model_group <- grp
    if(grp == group_set[1]) theta_posteriors_all <- theta_posteriors_grp else theta_posteriors_all <- rbind(theta_posteriors_all, theta_posteriors_grp)

    #save spde and INLA objects
    spde_all[[which(group_set == grp)]] <- spde_grp
    if(return_INLA_result) INLA_result_all[[which(group_set == grp)]] <- INLA_result_grp
  }

  # for testing/debugging
  #betas_df <- data.frame(classical=as.vector(GLM_result), Bayesian=as.vector(betas), task=rep(1:6, each=224), region=rep(c(rep(3,117),rep(4,107)),6))
  #ggplot(betas_df, aes(x=classical, y=Bayesian)) + facet_grid(region ~ task) + geom_point() + geom_smooth() + geom_abline(intercept=0, slope=1, col='red')

  #construct object to be returned
  if(return_INLA_result){
    result <- list(INLA_result = INLA_result_all,
                   spde_obj = spde_all,
                   mesh = list(n = spde_obj$spde$n.spde),
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   theta_posteriors = theta_posteriors,
                   call = match.call(),
                   GLM_result = GLM_result)
  } else {
    result <- list(INLA_result = NULL,
                   spde_obj = spde_obj,
                   mesh = list(n = spde_obj$spde$n.spde),
                   session_names = session_names,
                   beta_names = beta_names,
                   beta_estimates = beta_estimates,
                   theta_posteriors = theta_posteriors,
                   call = match.call(),
                   GLM_result = GLM_result)
  }


  class(result) <- "BayesGLM"

  if(!is.null(outfile)){
    if(!dir.exists(dirname(outfile))){dir.create(dirname(outfile))}
    ext <- strsplit(outfile, split=".", fixed = TRUE)[[1]]
    ext <- ext[length(ext)] #get last part of file name
    if(ext != "RDS"){
      outfile <- sub(ext, "RDS", outfile)
    }
    message('File saved at: ', outfile)
    save(result, file = outfile)
  }

  return(result)

}
