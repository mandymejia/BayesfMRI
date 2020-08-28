#' Create simulated slice data for the BayesGLM_surface function
#'
#' @param sessions a number
#' @param num_tasks a number
#' @param active_centers a matrix
#' @param active_size a vector of length \code{nrow(active_centers)}
#' @param vary_active logical, if there are multiple sessions, should the activation centers vary between sessions?
#' @param num_time length of the time series that is generated
#' @param binary_template (optional) a binary brain slice image
#'
#' @return a list..
#' @export
#' @importFrom neuRosim specifydesign specifyregion
#'
#' @examples
#' data <- simulate_slice_data()
simulate_slice_data <-
  function(num_sessions = 1,
           num_tasks = 2,
           active_centers = NULL,
           active_size = NULL,
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
  if(is.null(binary_template)) data("binary_template")
  if(is.null(active_centers)){
    bin_dims <- dim(binary_template)
    active_centers <- rbind(
      round(c(0.8,0.5) * bin_dims),
      round(c(0.25,0.5) * bin_dims),
      round(c(0.5,0.3) * bin_dims)
    )
  }
  if(is.null(active_size)) active_size <- sample(2:4, size = nrow(active_centers))
  y_i <- sapply(seq(num_sessions), function(i) {
    activation_regions <- mapply(function(cs,ss) {
      region_out <-
        specifyregion(
          dim = c(45, 54),
          coord = cs + c(sample(seq(-3,3),size = 1),sample(seq(-3,3),size = 1)),
          radius = ss,
          form = "sphere",
          fading = runif(1,0.1,0.5))
      return(region_out * rgamma(1,20, 20))
    },cs = split(active_centers, row(active_centers)),ss = sizes, SIMPLIFY = F)
    # Create beta_1
    beta_1 <- ifelse(binary_template == 0, NA, 0)
    beta_1 <- beta_1 + Reduce(`+`,sapply(activation_regions[-3],`*`,
                                         y = rgamma(1,20,20), simplify = F))
    # Create beta_2
    beta_2 <- ifelse(binary_template == 0, NA, 0)
    beta_2 <- beta_2 + Reduce(`+`,sapply(activation_regions[-1],`*`,
                                         y = rgamma(1,20,20), simplify = F))
    # Create the mean response
    task_1_means <- beta_1 %o% c(t1)
    task_2_means <- beta_2 %o% c(t2)
    y_means <- task_1_means + task_2_means
    # Use the mean responses to create the simulated response data
    y_t <- apply(y_means, 1:2, function(bv) {
      if(is.na(bv[1])) {
        return(rep(NA,length(bv)))
      } else {
        out <- 250 + bv + arima.sim(list(ar = 0.3), n = 200, sd = 2)
        return(out)
      }
    })
    y <- apply(y_t,1, identity)
    y_exclude <- apply(y,1, function(yv) any(is.na(yv)))
    y <- y[!y_exclude,]
    y <- t(y)
    return(list(y=y,betas = list(bbeta_1 = beta_1, bbeta_2 = beta_2)))
  }, simplify = F)

}
