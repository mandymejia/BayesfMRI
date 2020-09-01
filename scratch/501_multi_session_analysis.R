library(BayesfMRI)
# Setting the seed for reproducibility
set.seed(47401)
multi_data <-
  simulate_slice_data(
    num_sessions = 2,
    num_tasks = 2,
    active_centers = matrix(c(36, 28, 12, 28, 23, 16), 3, 2, byrow = T),
    active_size = 2:4,
    beta_weights = matrix(c(1, 1, 0, 0, 0.8, 0.8), 2, 3, byrow = T),
    vary_active = T,
    num_time = 200,
    binary_template = NULL
  )
# Create the mesh object for INLA
data("binary_template")
mesh <- make_2d_mesh(binary_template)
# Run the model in INLA
multi_session_result <- BayesGLM(data=multi_data$data, mesh=mesh, scale_BOLD = FALSE, scale_design = FALSE, verbose=F)
saveRDS(multi_session_result, "../scratch/501_multi_session_BayesGLM_results.rds")# Save the result
