test_that("BayesGLM_cifti", {
  library(brainSim)
  library(INLA)
  inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
  sim <- brainSim::simulate_cifti_multiple("/Applications/workbench",
                                           hemisphere = "left",ntime = 100,
                                           resamp_res = 1000)
  res <- BayesGLM_cifti(cifti_fname = list(sim$simulated_cifti$subj1_sess1_run1),
                        surfL_fname = sim$simulated_cifti$subj1_sess1_run1$surf$cortex_left,
                        surfR_fname = NULL,
                        brainstructures = "left",
                        design = sim$design,onsets = NULL, TR = 1,
                        scale_BOLD = T, scale_design = T,
                        GLM_method = "Bayesian", ar_order = 0, ar_smooth = 0,
                        nuisance = NULL, nuisance_include = c("drift","dHRF"),
                        session_names = NULL,resamp_res = NULL,num.threads = 4,
                        verbose = F,outfile = NULL,return_INLA_result = T,
                        avg_sessions = TRUE,trim_INLA = T,tol = 0.001)
})
