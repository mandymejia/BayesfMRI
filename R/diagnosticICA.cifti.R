#' Run diagnostic ICA based on CIFTI-format BOLD data and CIFTI-based template
#'
#' @param cifti_fname File path of CIFTI-format timeseries data (ending in .dtseries.nii).
#' @param templates LIST or VECTOR of templates, each the result of call to \code{estimate_template.cifti()} (one template for each group).
#' Either a list of objects of class \code{template.cifti} OR a vector of file path prefixes (the part before "_mean.dscalar.nii"
#' and "_var.dscalar.nii") of cifti files written out by \code{estimate_template.cifti()}.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).  For spatial
#'  modeling (performed if `surfL_fname` and/or `surfR_fname` provided), subcortical
#'  locations will be ignored.
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial
#' standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If TRUE, display progress of algorithm
#' @param write_dir Where should any output files be written? \code{NULL} (default) will write them to the current working directory.
#'
#' @return A list containing the subject IC estimates (class 'xifti'), the subject IC variance estimates (class 'xifti'), and the result of the model call to \code{diagnosticICA} (class 'dICA')
#' @export
#'
diagnosticICA.cifti <- function(cifti_fname,
                              templates,
                              brainstructures=c('left','right'),
                              scale=TRUE,
                              maxQ=NULL,
                              maxiter=100,
                              epsilon=0.001,
                              verbose=TRUE,
                              write_dir=NULL){

  if (is.null(write_dir)) { write_dir <- getwd() }

  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  do_left <- do_right <- do_sub <- FALSE
  if('left' %in% brainstructures) do_left <- TRUE
  if('right' %in% brainstructures) do_right <- TRUE
  if('subcortical' %in% brainstructures) do_sub <- TRUE

  G <- length(templates)
  if(G>5) stop(paste0('Length of templates is ',G,' which is a large number of groups. Check that your templates argument is correct.'))
  cat(paste0('Number of group templates provided: ',G,'\n'))

  # GET TEMPLATE MEAN AND VARIANCE FOR EACH GROUP
  template_class <- sapply(templates, class)
  if(length(unique(template_class))>1) stop('All elements of templates argument must be of the same class')
  template_mean_mat <- template_var_mat <- vector('list', length=G)

  for(g in 1:G){

    #Obtain template_mean_g and template_var_g
    if(template_class[1]=='template.cifti'){
        template_mean_g <- templates[[g]]$template_mean #class xifti
        template_var_g <- templates[[g]]$template_var #class xifti
    } else if(template_class[1]=='character'){
      fname_mean <- paste0(templates[[g]],'_mean.dscalar.nii')
      fname_var <- paste0(templates[[g]],'_var.dscalar.nii')
      if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
      if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
      template_mean_g <- read_cifti(fname_mean, brainstructures=brainstructures)
      template_var_g <- read_cifti(fname_mean, brainstructures=brainstructures)
    } else {
      stop('template argument must be an object of class template.cifti or file path prefix to result of estimate_template.cifti() (same as out_fname argument passed to estimate_template.cifti().')
    }

    #Extract data matrix from template_mean_g and template_var_g
    template_mean_mat[[g]] <- rbind(template_mean_g$data$cortex_left, template_mean_g$data$cortex_right, template_mean_g$data$subcort)
    template_var_mat[[g]] <- rbind(template_var_g$data$cortex_left, template_var_g$data$cortex_right, template_var_g$data$subcort)
  }

  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(cifti_fname)) stop(paste0('The BOLD timeseries file ',cifti_fname,' does not exist.'))
  BOLD <- read_cifti(cifti_fname, brainstructures = brainstructures)
  BOLD_mat <- rbind(BOLD$data$cortex_left, BOLD$data$cortex_right, BOLD$data$subcort)

  #set up xifti objects for IC mean and variance estimates
  clear_data <- function(x){
    if(!is.null(x$data$cortex_left)) x$data$cortex_left <- matrix(0, nrow(x$data$cortex_left), 1)
    if(!is.null(x$data$cortex_right)) x$data$cortex_right <- matrix(0, nrow(x$data$cortex_right), 1)
    if(!is.null(x$data$subcort)) x$data$subcort <- matrix(0, nrow(x$data$subcort), 1)
    return(x)
  }
  subjICmean_xifti <- subjICvar_xifti <- clear_data(BOLD)


  # CALL DIAGNOSTIC ICA FUNCTION
  result_mod <- diagnosticICA(template_mean = template_mean_mat,
                            template_var = template_var_mat,
                            BOLD = BOLD_mat,
                            scale = scale,
                            maxQ=maxQ,
                            maxiter=maxiter,
                            epsilon=epsilon,
                            verbose=verbose)

  # MAP ESTIMATES AND VARIANCE TO XIFTI FORMAT
  n_left <- n_right <- n_sub <- 0
  if(do_left) {
    n_left <- nrow(subjICmean_xifti$data$cortex_left)
    subjICmean_xifti$data$cortex_left <- result_mod$subjICmean[1:n_left,]
    subjICvar_xifti$data$cortex_left <- result_mod$subjICvar[1:n_left,]
  }
  if(do_right) {
    n_right <- nrow(subjICmean_xifti$data$cortex_right)
    subjICmean_xifti$data$cortex_right <- result_mod$subjICmean[n_left+(1:n_right),]
    subjICvar_xifti$data$cortex_right <- result_mod$subjICvar[n_left+(1:n_right),]
  }
  if(do_sub) {
    n_sub <- nrow(subjICmean_xifti$data$subcort)
    subjICmean_xifti$data$subcort <- result_mod$subjICmean[n_left + n_right + (1:n_sub),]
    subjICvar_xifti$data$subcort <- result_mod$subjICvar[n_left + n_right + (1:n_sub),]
  }

  # RETURN XIFTI RESULTS AND MODEL RESULT
  list(subjICmean_xifti = subjICmean_xifti,
       subjICvar_xifti = subjICvar_xifti,
       model_result = result_mod)

}

