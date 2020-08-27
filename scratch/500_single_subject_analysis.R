set.seed(47401)
library(BayesfMRI)
# Load the binary template image
data(binary_template)
# Create the activation profiles
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
t1 <- (t1 - mean(t1)) / max(t1)
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
t2 <- (t2 - mean(t2)) / max(t2)
# Region 1
r1 <-
  specifyregion(
    dim = c(45, 54),
    coord = c(36, 28),
    radius = 2,
    form = "sphere",
    fading = 0.2
  )
# Region 2
r2 <-
  specifyregion(
    dim = c(45, 54),
    coord = c(12, 28), # Center of the activation
    radius = 3, # How big is the activation regions
    form = "sphere", # The activation are will be circular
    fading = 0.2 # How fast does the activation region fade (0 is no fade, 1 is fast fade)
  )
# Region 3
r3 <-
  specifyregion(
    dim = c(45, 54),
    coord = c(23, 16),
    radius = 4,
    form = "sphere",
    fading = 0.2
  )
# Create beta_1
beta_1 <- ifelse(binary_template == 0, NA, 0)
beta_1 <- beta_1 + r1 + r2
# Create beta_2
beta_2 <- ifelse(binary_template == 0, NA, 0)
beta_2 <- beta_2 + 0.8*(r2 + r3)
# Create the response data
task_1_means <- beta_1 %o% c(t1)
task_2_means <- beta_2 %o% c(t2)
y_means <- task_1_means + task_2_means
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
# Get data ready for INLA
library(INLA)
in_mask <- which(binary_template == 1, arr.ind = T)
in_mask <- in_mask[,2:1]
boundary <- inla.nonconvex.hull(in_mask, resolution = 100)
mesh <- inla.mesh.2d(loc = in_mask, boundary = boundary, max.edge = c(2,4))
# Make a session list (for INLA)
session <- list(
  BOLD = y,
  design = cbind(t1,t2)
)
# Make the data list (for INLA)
data <- list(
  single_session = session
)
# Make a conversion matrix to transform back to an image
convert_mat <- inla.spde.make.A(mesh = mesh, loc = in_mask)
# Run the classical GLM
single_subject_classical <- classicalGLM(data = data, scale_BOLD = T, scale_design = T)
# Register PARDISO
inla.setOption(pardiso.license="~/licenses/pardiso.lic")
inla.pardiso.check()
# Run the Bayes
single_subject_result <- BayesGLM_surface(data = data, mesh = mesh, verbose = F)
# Save the results
saveRDS(single_subject_classical, "500_single_subject_classical_results.rds")
saveRDS(single_subject_result, "500_single_subject_BayesGLM_results.rds")
