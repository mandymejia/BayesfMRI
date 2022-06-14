#' BayesGLM for 3D volume
#'
#' Applies spatial Bayesian GLM to task fMRI data for 3D subcortical volumes
#'
#' The subcortical data is separated into regions, whose sizes range
#'  from approximately 100 voxels to approximately 9000 voxels.  Smaller regions
#'  are grouped together to improve model fit.
#'  The \code{groups_df} argument specifies which regions are grouped together.
#'  This argument should be a data frame with R rows (the number of regions) and
#'  three columns: label, region, and group.
#'  The label column is the numerical identifier of each region; the region
#'  column contains the region names, and the group column contains the model
#'  group assignments (e.g. 1,2,3). Regions to be excluded
#'  from analysis are indicated by NA in the group assignment.
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#'  List element names represent session names.
#' @param locations Vx3 matrix of x,y,z coordinates of each voxel
#' @param labels Vector of length V of region labels
#' @param groups_df Data frame indicating the name and model group of each region.  See Details.
#' @param scale Option for scaling the BOLD response.
#'
#' 	If \code{"auto"} (default), will use mean scaling except if demeaned data
#' 	is detected, in which case sd scaling will be used instead.
#'
#' 	\code{"mean"} scaling will scale the data to percent local signal change.
#'
#' 	\code{"sd"} scaling will scale the data by local standard deviation.
#'
#' 	\code{"none"} will only center the data, not scale it.
#' @inheritParams return_INLA_result_Param_FALSE
#' @param outfile File name where results will be written (for use by \code{BayesGLM_grp}).
#' @param GLM If TRUE, classical GLM estimates will also be returned
#' @inheritParams num.threads_Param
#' @inheritParams verbose_Param_inla
#'
#' @return A list containing...
#' @importFrom stats as.formula
#'
#' @export
BayesGLM_vol3D <- function(data, locations, labels, groups_df, scale=c("auto", "mean", "sd", "none"),
  return_INLA_result=FALSE, outfile = NULL, GLM = TRUE, num.threads = 6, verbose=FALSE){

  check_INLA(FALSE)

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
    warning('No value supplied for outfile, which is required for group modeling (see `help(BayesGLM2)`).')
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
  if(!inherits(groups_df, 'data.frame')) stop('The groups_df argument should be a data frame with the following columns: label, region, group.')
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

    cat(paste0('Estimating Model 1 (', paste(name_set_grp, collapse = ', '), ')'),"\n")

    spde_grp <- create_spde_vol3D(locs=locs_grp, labs=labels_grp, lab_set=label_set_grp)
    #plot(spde_grp)
    spde <- spde_grp$spde

    #collect data and design matrices
    y_all <- c()
    X_all_list <- NULL

    for(s in 1:n_sess){

      #extract and mask BOLD data for current session
      BOLD_s <- data[[s]]$BOLD[,inds_grp]

      #scale data to represent % signal change (or just center if scale=="none")
      BOLD_s <- scale_timeseries(t(BOLD_s), scale=scale, transpose = FALSE)
      design_s <- scale(data[[s]]$design, scale=FALSE) #center design matrix to eliminate baseline

      #regress nuisance parameters from BOLD data and design matrix
      if('nuisance' %in% names(data[[s]])){
        design_s <- data[[s]]$design
        nuisance_s <- data[[s]]$nuisance
        y_reg <- nuisance_regression(BOLD_s, nuisance_s)
        X_reg <- nuisance_regression(design_s, nuisance_s)
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
    replicates_list <- organize_replicates(n_sess=n_sess, beta_names = beta_names, mesh = spde_grp)
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

#' Extracts posterior density estimates for hyperparameters
#'
#' @inheritSection INLA_Description INLA Requirement
#'
#' @param object An object of class \code{"inla"}, a result of a call to
#'  \code{inla()}
#' @param spde The model used for the latent fields in the \code{inla()} call,
#'  an object of class \code{"inla.spde"}
#'
#' @return Long-form data frame containing posterior densities for the
#'  hyperparameters associated with each latent field
#'
#' @export
get_posterior_densities_vol3D <- function(object, spde){

  check_INLA(FALSE)

  hyper_names <- names(object$marginals.hyperpar)

  for(h in hyper_names){
    df.h <- INLA::inla.extract.el(object$marginals.hyperpar, h)
    df.h <- as.data.frame(df.h)
    names(df.h) <- c('x','y')
    df.h$hyper_name <- h
    df.h$beta_name <- ifelse(grepl('bbeta',h), #if bbeta appears in name
                             gsub('.+bbeta','bbeta',h), #rename as bbeta*
                             NA) #NA if this is the precision hyperparameter
    df.h$theta_name <- ifelse(grepl('Theta',h), #if bbeta appears in name
                              gsub(' for.+','',h), #rename as Theta*
                              NA) #NA if this is the precision hyperparameter

    if(h==hyper_names[1]) df <- df.h else df <- rbind(df, df.h)
  }
  return(df)
}



### More updated function to apply BayesGLM to a single slice.  Designed mainly for simulations.

#'  BayesGLM for 2D slice
#'
#'  Spatial Bayesian GLM for fMRI task activation on 2d slice volumetric data
#'
#' @param BOLD A list of sessions, each with a three-dimensional array in which
#'   the first two dimensions correspond to the size of the fMRI slice in space
#'   and the last dimension corresponds to time
#' @param binary_mask (optional) a binary brain slice image used to mask
#'   the BOLD data and make a more efficient network mesh for the
#'   neighborhood definitions
#' @param design,onsets,TR Either provide \code{design}, or provide both \code{onsets} and \code{TR}.
#'
#'   \code{design} is a \eqn{T x K} task design matrix (or list of such
#'   matrices, for multiple-session modeling) with column names representing
#'   tasks. Each column represents the expected BOLD response due to each task,
#'   a convolution of the hemodynamic response function (HRF) and the task
#'   stimulus. Note that the scale of the regressors will affect the scale and
#'   interpretation of the beta coefficients, so imposing a proper scale (e.g.,
#'   set maximum to 1) is recommended.
#'
#'   \code{onsets} is a matrix of onsets (first column) and durations (second column)
#'   for each task in seconds, organized as a list where each element of the
#'   list corresponds to one task. Names of list should be task names. (Or for
#'   multi-session modeling, a list of such lists.)
#'
#'   \code{TR} is the temporal resolution of the data in seconds.
#' @param nuisance (Optional) A TxJ matrix of nuisance signals (or list of such
#'   matrices, for multiple-session modeling).
#' @param dHRF Set to \code{1} to add the temporal derivative of each column
#'  in the design matrix, \code{2} to add the second derivatives too, or
#'  \code{0} to not add any columns. Default: \code{1}.
#' @param hpf,DCT Add DCT bases to \code{nuisance} to apply a temporal 
#'  high-pass filter to the data? Only one of these arguments should be provided.
#'  \code{hpf} should be the filter frequency; if it is provided, \code{TR}
#'  must be provided too. The number of DCT bases to include will be computed
#'  to yield a filter with as close a frequency to \code{hpf} as possible.
#'  Alternatively, \code{DCT} can be provided to directly specify the number 
#'  of DCT bases to include.
#'  
#'  Default: \code{DCT=4} (use four DCT bases for high-pass filtering; for
#'  typical \code{TR} this amounts to lower filter frequency than the
#'  approximately .01 Hz used in most studies.)
#' @inheritParams scale_BOLD_Param
#' @inheritParams scale_design_Param
#' @inheritParams num.threads_Param
#' @param GLM_method Either 'Bayesian' for spatial Bayesian GLM only,
#'   'classical' for the classical GLM only, or 'both' to return both classical
#'   and Bayesian estimates of task activation.
#' @param session_names (Optional) A vector of names corresponding to each
#'   session.
#' @inheritParams return_INLA_result_Param_TRUE
#' @param outfile (Optional) File name (without extension) of output file for
#'   BayesGLM result to use in Bayesian group modeling.
#' @inheritParams verbose_Param_inla
#' @inheritParams contrasts_Param_inla
#' @inheritParams avg_sessions_Param
#' @param trim_INLA (logical) should the \code{INLA_result} objects within the
#'   result be trimmed to only what is necessary to use `id_activations()`? Default: `TRUE`.
#'
#' @importFrom utils head
#'
#' @return An object of class \code{"BayesGLM"}, a list containing...
#'
#' @export
BayesGLM_slice <- function(
  BOLD,
  design = NULL,
  onsets=NULL,
  TR=NULL,
  nuisance=NULL,
  dHRF=c(0, 1, 2),
  hpf=NULL,
  DCT=if(is.null(hpf)) {4} else {NULL},
  binary_mask = NULL,
  scale_BOLD = c("auto", "mean", "sd", "none"),
  scale_design = TRUE,
  num.threads = 4,
  GLM_method = 'both',
  session_names = NULL,
  return_INLA_result = TRUE,
  outfile = NULL,
  verbose = FALSE,
  contrasts = NULL,
  avg_sessions = TRUE,
  trim_INLA = TRUE) {

  do_Bayesian <- (GLM_method %in% c('both','Bayesian'))
  do_classical <- (GLM_method %in% c('both','classical'))

  # Check nuisance arguments.
  dHRF <- as.numeric(
    match.arg(as.character(dHRF), as.character(c(0, 1, 2)))
  )
  if (!is.null(DCT)) { 
    stopifnot(is.numeric(DCT) && length(DCT)==1 && DCT>=0 && DCT==round(DCT)) 
    if (DCT==0) { DCT <- NULL }
  }
  if (!is.null(hpf)) { 
    stopifnot(is.numeric(hpf) && length(hpf)==1 && hpf>=0)
    if (hpf==0) { hpf <- NULL }
  }

  check_INLA(require_PARDISO=do_Bayesian)

  image_dims <- head(dim(BOLD[[1]]),-1)
  if (is.null(binary_mask))
    binary_mask <- matrix(1, nrow = image_dims[1], ncol = image_dims[2])

  mesh <- make_slice_mesh(binary_mask)

  # Name sessions and check compatibility of multi-session arguments
  n_sess <- length(BOLD)
  if(n_sess == 1 & avg_sessions) avg_sessions <- FALSE
  if(n_sess==1){
    if(is.null(session_names)) session_names <- 'single_session'
  } else {
    if(is.null(session_names)) session_names <- paste0('session', 1:n_sess)
  }
  if(length(session_names) != n_sess)
    stop('If session_names is provided, it must be of the same length as BOLD')

  cat('\n SETTING UP DATA \n')

  if(is.null(design)) {
    make_design <- TRUE
    design <- vector('list', length=n_sess)
  } else {
    make_design <- FALSE
  }

  for(ss in 1:n_sess){
    if(make_design){
      cat(paste0('    Constructing design matrix for session ', ss, '\n'))
      design[[ss]] <- make_HRFs(onsets[[ss]], TR=TR, duration=ntime, deriv=dHRF)
    }
  }

  ### Check that design matrix names consistent across sessions
  if(n_sess > 1){
    tmp <- sapply(design, colnames)
    tmp <- apply(tmp, 1, function(x) length(unique(x)))
    if(max(tmp) > 1)
      stop('task names must match across sessions for multi-session modeling')
  }

  cat('\n RUNNING MODEL \n')

  classicalGLM <- NULL
  BayesGLM <- NULL

  ### FORMAT DESIGN MATRIX
  for(ss in 1:n_sess){
    if(scale_design){
      design[[ss]] <- scale_design_mat(design[[ss]])
    } else {
      design[[ss]] <- scale(design[[ss]], scale=FALSE) #center design matrix
      # to eliminate baseline
    }
  }

  ### ADD ADDITIONAL NUISANCE REGRESSORS
  DCTs <- vector("numeric", n_sess)
  for (ss in 1:n_sess) {
    # DCT highpass filter
    if (!is.null(hpf) || !is.null(DCT)) {
      # Get the num. of bases for this session.
      if (!is.null(hpf)) {
        DCTs[ss] <- round(dct_convert(ntime[ss], TR, f=hpf))
      } else {
        DCTs[ss] <- DCT
      }
      # Generate the bases and add them.
      DCTb_ss <- dct_bases(ntime[ss], DCTs[ss])
      if (DCTs[ss] > 0) {
        if (!is.null(nuisance)) {
          nuisance[[ss]] <- cbind(nuisance[[ss]], DCTb_ss)
        } else {
          nuisance[[ss]] <- DCTb_ss
        }
      }
    }
  }

  scale_design <- F # This is done to prevent double-scaling in BayesGLM

  #set up session list
  # mat_BOLD <- sapply(BOLD, function(y_t) {
  #   # Remove any NA voxels and output the response as a matrix
  #   y <- apply(y_t,3, identity)
  #   y_exclude <- apply(y,1, function(yv) any(is.na(yv)))
  #   y <- y[!y_exclude,]
  #   y <- t(y)
  # }, simplify = F)
  session_data <- vector('list', n_sess)
  names(session_data) <- session_names
  for(ss in 1:n_sess){
    sess <- list(BOLD = BOLD[[ss]], design=design[[ss]])
    if(!is.null(nuisance)) sess$nuisance <- nuisance[[ss]]
    session_data[[ss]] <- sess
  }

  ### FIT GLM(s)

  if(do_classical) classicalGLM_out <- classicalGLM(session_data,
                                                    scale_BOLD=scale_BOLD,
                                                    scale_design = scale_design)
  if(do_Bayesian) {
    BayesGLM_out <- BayesGLM(session_data,
                             mesh = mesh,
                             scale_BOLD=scale_BOLD,
                             scale_design = scale_design,
                             num.threads = num.threads,
                             return_INLA_result = return_INLA_result,
                             outfile = outfile,
                             verbose = verbose,
                             avg_sessions = avg_sessions,
                             trim_INLA = trim_INLA)

    # Create a conversion matrix
    in_binary_mask <- which(binary_mask == 1, arr.ind = T)
    in_binary_mask <- in_binary_mask[,2:1]
    convert_mat_A <- INLA::inla.spde.make.A(mesh = mesh, loc = in_binary_mask)
    # Extract the point estimates
    point_estimates <- sapply(session_names, function(sn){
      as.matrix(convert_mat_A %*% BayesGLM_out$beta_estimates[[sn]])
    }, simplify = F)
    if(avg_sessions)
      avg_point_estimates <- BayesGLM_out$avg_beta_estimates
  }

  classical_slice <- Bayes_slice <- vector('list', n_sess)
  names(classical_slice) <- names(Bayes_slice) <- session_names
  for(ss in 1:n_sess){
    num_tasks <- ncol(design[[ss]])
    if(do_classical){
      classical_slice[[ss]] <- sapply(seq(num_tasks), function(tn) {
        image_coef <- binary_mask
        image_coef[image_coef == 1] <- classicalGLM_out[[ss]][,tn]
        image_coef[binary_mask == 0] <- NA
        return(image_coef)
      },simplify = F)
    }
    if(do_Bayesian){
      mat_coefs <- point_estimates[[ss]]
      Bayes_slice[[ss]] <- sapply(seq(num_tasks), function(tn) {
        image_coef <- binary_mask
        image_coef[image_coef == 1] <- mat_coefs[,tn]
        image_coef[binary_mask == 0] <- NA
        return(image_coef)
      },simplify = F)
      if(n_sess > 1 & avg_sessions) {
        Bayes_slice$avg_over_sessions <- sapply(seq(num_tasks), function(tn) {
          image_coef <- binary_mask
          image_coef[image_coef == 1] <- avg_point_estimates[,tn]
          image_coef[binary_mask == 0] <- NA
          return(image_coef)
        },simplify = F)
      }
    }
  }

  if (do_Bayesian) {
    beta_names <- BayesGLM_out$beta_names
  } else {
    beta_names <- NULL
    BayesGLM_out <- NULL
  }

  if(!do_classical)
    classicalGLM_out <- NULL

  result <- list(session_names = session_names,
                 beta_names = beta_names,
                 betas_Bayesian = Bayes_slice,
                 betas_classical = classical_slice,
                 GLMs_Bayesian = BayesGLM_out,
                 GLMs_classical = classicalGLM_out,
                 design = design,
                 mask = binary_mask)
  class(result) <- "BayesGLM_slice"
  return(result)
}


