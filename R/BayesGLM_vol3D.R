#' Applies spatial Bayesian GLM to task fMRI data for 3D subcortical volumes
#'
#' @param data A list of sessions, where each session is a list with elements
#' BOLD, design and nuisance.  See \code{?create.session} and \code{?is.session} for more details.
#' List element names represent session names.
#' @param spde_obj Spatial partial differential equation for subcortical region of interest. See \code{?create_spde_vol3D} for more details.
#' @param scale If TRUE, scale timeseries data so estimates represent percent signal change.  Else, do not scale.
#' @param return_INLA_result If TRUE, object returned will include the INLA model object (can be large).  Default is TRUE. Required for running \code{id_activations} on \code{BayesGLM} model object.
#' @param outfile File name where results will be written (for use by \code{BayesGLM_grp}).
#' @param GLM If TRUE, classical GLM estimates will also be returned
#'
#' @return A list containing...
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
#' @examples \dontrun{}
BayesGLM_vol3D <- function(data, spde_obj, scale=TRUE, return_INLA_result=TRUE, outfile = NULL, GLM = TRUE, num.threads = 6){

  ### HERE

  ### Take in as arguments locations, labels, and label names/groupings
  ### Run create_spde_vol3D for each grouping

  #check whether data is a list OR a session (for single-session analysis)
  #check whether each element of data is a session (use is.session)
  # V = number of data locations
  # T = length of time series for each session (vector)
  # K = number of unique tasks in all sessions

  #need to check that sessions are consistent in terms of V, K?

  #INLA:::inla.dynload.workaround() #avoid error on creating mesh

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
    warning('No value supplied for outfile, which is required for post-hoc group modeling.')
  }

  # Create SPDE
  if(is.null(spde_obj)){
    stop("SPDE object is required. See ?create_spde_vol3D")
  }
  spde <- spde_obj$spde

  #collect data and design matrices
  y_all <- c()
  X_all_list <- NULL

  for(s in 1:n_sess){

    #extract and mask BOLD data for current session
    BOLD_s <- data[[s]]$BOLD

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
  replicates_list <- organize_replicates(n_sess=n_sess, n_task=K, mesh = spde_obj)
  betas <- replicates_list$betas
  repls <- replicates_list$repls

  #organize the formula and data objects
  #formula <- make_formula(beta_names = names(betas), repl_names = names(repls), hyper_initial = c(-2,2))
  #formula <- as.formula(formula)

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

  # estimate model using INLA
  INLA_result <- estimate_model(formula=formula, data=model_data, A=model_data$X, spde=spde, prec_initial=1, num.threads = num.threads)
  # INLA_result <- readRDS('/Users/sarahryan/Desktop/StuffWithMandy/myINLA_result.RDS')

  #extract useful stuff from INLA model result
  beta_estimates <- extract_estimates(object=INLA_result, session_names=session_names) #posterior means of latent task field
  theta_posteriors <- get_posterior_densities_vol3D(object=INLA_result, spde) #hyperparameter posterior densities

  if(GLM){
    # estimate model using GLM
    GLM_result <- solve(t(design_s) %*% design_s) %*% t(design_s) %*% BOLD_s
    GLM_result <- t(GLM_result)
  }else{
    GLM_result <- NULL
  }

  #construct object to be returned
  if(return_INLA_result){
    result <- list(INLA_result = INLA_result,
                   spde_obj = spde_obj,
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
    ext <- ext[length(ext)]
    if(ext != "RDS"){
      outfile <- sub(ext, "RDS", outfile)
    }
    message('File saved at: ', outfile)
    save(result, file = outfile)
  }

  return(result)

}
