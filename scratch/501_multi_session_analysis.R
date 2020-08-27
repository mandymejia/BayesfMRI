# Library and seed
set.seed(47401)
library(BayesfMRI)
# Load the binary template image
data(binary_template)
# Create the covariates
library(neuRosim)
# Task 1
t1 <-
  specifydesign(
    onsets = seq(0, 200, by = 40),
    durations = 1,
    totaltime = 200,
    TR = 1,
    effectsize = 1.3,
    conv = "double-gamma",
    param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
  )
# Task 2
t2 <-
  specifydesign(
    onsets = seq(20, 200, by = 40),
    durations = 1,
    totaltime = 200,
    TR = 1,
    effectsize = 1.3,
    conv = "double-gamma",
    param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
  )

# Set the values for the activations
num_session <- 2
centers <- matrix(c(36,28, 12,28, 23,16), nrow = 3, byrow = T)
sizes <- 2:4
mean_smoothness <- 0.2

# Create the responses
y_i <- sapply(seq(num_session), function(i) {
  activation_regions <- mapply(function(cs,ss) {
    region_out <-
      specifyregion(
        dim = c(45, 54),
        coord = cs + c(sample(seq(-3,3),size = 1),sample(seq(-3,3),size = 1)),
        radius = ss,
        form = "sphere",
        fading = runif(1,0.1,0.5))
    return(region_out * rgamma(1,20, 20))
    },cs = split(centers, row(centers)),ss = sizes, SIMPLIFY = F) 
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
  return(y)
}, simplify = F)
# Create the sessions
data <- sapply(y_i, function(yi) {
  out <- list(BOLD = yi,
              design = cbind(t1,t2))
  return(out)
},simplify = F)
names(data) <- paste("session",1:2,sep="_")
# Create the mesh object for INLA
library(INLA)
in_mask <- which(binary_template == 1, arr.ind = T)
in_mask <- in_mask[,2:1]
boundary <- inla.nonconvex.hull(in_mask, resolution = 100)
mesh <- inla.mesh.2d(loc = in_mask, boundary = boundary, max.edge = c(2,4))
# Run the model in INLA
multi_session_result <- BayesGLM_surface(data=data, mesh=mesh, verbose=F)
# Save the result
saveRDS(multi_session_result, "~/github/BayesfMRI/vignettes/501_multi_session_BayesGLM_results.rds")
