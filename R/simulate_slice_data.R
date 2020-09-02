#' Create simulated slice data for the BayesGLM_surface function
#'
#' @param num_sessions a number
#' @param num_tasks a number
#' @param active_centers a matrix, with the number of rows corresponding to the
#'   number of activation regions. There should be two columns, corresponding
#'   to the x- and y-coordinates for the centers.
#' @param active_size a vector of length \code{nrow(active_centers)}
#' @param beta_weights a matrix, with the number of rows equal to
#'   \code{num_tasks} and the number of columns equal to the number of rows in
#'   \code{active_centers}
#' @param vary_active logical, if there are multiple sessions, should the
#'   activation centers vary between sessions?
#' @param num_time length of the time series that is generated
#' @param binary_template (optional) a binary brain slice image
#'
#' @return A list in which the first element, \code{data}, is a list of
#'   sessions, each with components \code{BOLD} (a matrix) and \code{design}
#'   (a matrix). The second list element is named \code{betas}, which is a list
#'   of sessions, each containing a list of the true values of the coefficients
#'   used to generate the data to be analyzed.
#' @export
#' @importFrom neuRosim specifydesign specifyregion
#' @import stats
#'
#' @examples
#' data <- simulate_slice_data()
simulate_slice_data <-
  function(num_sessions = 1,
           num_tasks = 2,
           active_centers = NULL,
           active_size = NULL,
           beta_weights = NULL,
           vary_active = T,
           num_time = 200,
           binary_template = NULL) {
  epoch_length <- floor(num_time / 5)
  onset_start <- seq(0,floor(epoch_length / 2), length.out = num_tasks)
  tasks <- sapply(seq(num_tasks), function(tt) {
    task_n <- neuRosim::specifydesign(
      onsets = seq(onset_start[tt], num_time, by = epoch_length),
      durations = 1,
      totaltime = num_time,
      TR = 1,
      effectsize = 1.3,
      conv = "double-gamma",
      param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
    )
    return(task_n)
  }, simplify = T)
  # Create the responses
  if(is.null(binary_template)) {
    data("binary_template", package = "BayesfMRI",envir = environment())
  }
  if(is.null(active_centers)){
    bin_dims <- dim(binary_template)
    active_centers <- rbind(
      round(c(0.8,0.5) * bin_dims),
      round(c(0.25,0.5) * bin_dims),
      round(c(0.5,0.3) * bin_dims)
    )
  }
  if(is.null(active_size)) {
    active_size <- sample(2:4, size = nrow(active_centers))
  }
  if(is.null(beta_weights)) {
    largest_margin <- max(num_tasks, nrow(active_centers))
    base_weights <- stats::toeplitz(rev(seq(largest_margin)) / largest_margin)
    beta_weights <- base_weights[seq(num_tasks),seq(nrow(active_centers))]
  }
  y_i <- sapply(seq(num_sessions), function(i) {
    # Make the active regions
    activation_regions <- mapply(function(cs,ss) {
      if(vary_active) {
        region_out <-
          neuRosim::specifyregion(
            dim = c(45, 54),
            coord = cs + c(sample(seq(-3,3),size = 1),sample(seq(-3,3),size = 1)),
            radius = ss,
            form = "sphere",
            fading = runif(1,0.1,0.5))
        region_out <- region_out * rgamma(1,20, 20)
      } else {
        region_out <-
          neuRosim::specifyregion(
            dim = c(45, 54),
            coord = cs,
            radius = ss,
            form = "sphere",
            fading = 0.3)
      }
      return(region_out)
    },cs = split(active_centers, row(active_centers)),ss = active_size,
    SIMPLIFY = F)
    # Make the coefficients using the active regions
    beta_coefficients <- sapply(seq(num_tasks), function(j) {
      beta_components <- mapply(`*`, x = beta_weights[j,],
                                y = activation_regions, SIMPLIFY = F)
      beta <- Reduce(`+`, beta_components)
      beta[binary_template == 0] <- NA
      return(beta)
    }, simplify = F)
    # Make the mean values for the BOLD data
    task_means <- mapply(`%o%`,X = beta_coefficients,
                         Y = split(tasks, col(tasks)), SIMPLIFY = F)
    y_means <- Reduce(`+`, task_means)
    # Use the mean responses to create the simulated response data
    y_t <- apply(y_means, seq(length(dim(y_means)) - 1), function(bv) {
      if(is.na(bv[1])) {
        return(rep(NA,length(bv)))
      } else {
        out <- 250 + bv + arima.sim(list(ar = 0.3), n = 200, sd = 2)
        return(out)
      }
    })
    # Remove any NA voxels and output the response as a matrix
    y <- apply(y_t,1, identity)
    y_exclude <- apply(y,1, function(yv) any(is.na(yv)))
    y <- y[!y_exclude,]
    y <- t(y)
    # Return everything
    return(list(session = list(BOLD=y,design = tasks), betas = beta_coefficients))
  }, simplify = F)
  names(y_i) <- paste("session",seq(num_sessions), sep ="_")
  data <- sapply(y_i, `[[`, i = "session", simplify = F)
  betas <- sapply(y_i, `[[`, i = "betas", simplify = F)
  return(list(data = data, betas = betas))
}
