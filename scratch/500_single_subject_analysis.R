library(BayesfMRI)
set.seed(47401)
simulated_data <-
  simulate_slice_data(
    num_sessions = 1,
    num_tasks = 2,
    active_centers = matrix(c(36, 28, 12, 28, 23, 16), 3, 2, byrow = T),
    active_size = 2:4,
    beta_weights = matrix(c(1, 1, 0, 0, 0.8, 0.8), 2, 3, byrow = T),
    vary_active = FALSE,
    num_time = 200,
    binary_template = NULL
  )
# Run the classical GLM
single_subject_classical <- classicalGLM(data = simulated_data$data,
                                         scale_BOLD = F, scale_design = F)
# Register PARDISO
inla.setOption(pardiso.license="~/licenses/pardiso.lic")
inla.pardiso.check()
# Run the Bayes
data("binary_template")
mesh <- make_2d_mesh(binary_template)
single_subject_result <- BayesGLM(data = simulated_data$data, mesh = mesh, scale_BOLD = F,
                                  scale_design = F, verbose = F)
# Save the results
saveRDS(single_subject_classical, "../scratch/500_single_subject_classical_results.rds")
saveRDS(single_subject_result, "../scratch/500_single_subject_BayesGLM_results.rds")
