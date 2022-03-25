# This is a test script for the BayesfMRI package. If this runs successfully,
# that means that the package is working as intended.

# Original ----------------------------
wb_path <- "/Applications/workbench"

# For Damon ---------------------------
# wb_path <- "~/Desktop/workbench"

library(brainSim)
library(ciftiTools)
ciftiTools.setOption("wb_path", wb_path)
library(BayesfMRI)
save_path <- "~/Desktop"
set.seed(47408)
num_errs <- 0
for(num_runs in 1:2){
  cat("SIMULATING DATA FOR",num_runs,"RUN(S)\n")
  sim_data <-
    simulate_cifti_multiple(
      wb_path = wb_path,
      brainstructure = "both",
      n_subjects = 1,
      n_sessions = 1,
      n_runs = num_runs,
      ntasks = 2,
      ntime = 100,
      resamp_res = 1000,
      max_amplitude = 2,
      onsets = NULL,
      durations = NULL,
      TR = 1,
      ar_error = c(0.35, 0.15)
  )
  # True values for comparing performance
  true_coefs <- sapply(c("left","right"), function(h){
    out <- sapply(sim_data$coef_cifti, function(sess_coeffs) {
      c(sess_coeffs$data[[paste0("cortex_",h)]])
    }, simplify = F)
    # If there are multiple sessions, and we're interested in the average
    out <- Reduce(cbind,out)
    if("matrix" %in% class(out)) out <- apply(out,1,mean)
    return(out)
  })
  for(hem in c("left","right")) {
    cat("....",toupper(hem),"HEMISPHERE\n")
    for(PW in rev(c("not prewhitened","prewhitened"))) {
      # stop break for testing (top) ----
      # stop()
      cat(".... DATA",toupper(PW),"\n")
      if (PW == "prewhitened") {
        ar_order = 6
        ar_smooth = 6
      } else {
        ar_order <- ar_smooth <- 0
      }
      if(num_runs == 1){
        des <- sim_data$design
      } else {
        des <- rep(list(sim_data$design), num_runs)
      }
      result_obj <-
        try(
          BayesGLM_cifti(
            cifti_fname = sim_data$simulated_cifti,
            surfL_fname = sim_data$simulated_cifti[[1]]$surf$cortex_left,
            surfR_fname = sim_data$simulated_cifti[[1]]$surf$cortex_right,
            brainstructures = hem,
            design = des,
            onsets = NULL,
            TR = 1,
            nuisance = NULL,
            nuisance_include = c("drift", "dHRF"),
            scale_BOLD = TRUE,
            scale_design = TRUE,
            Bayes = TRUE, # This is the construction in BayesfMRI 1.9
            ar_order = ar_order,
            ar_smooth = ar_smooth,
            resamp_res = NULL,
            num.threads = 6,
            verbose = FALSE,
            outfile = NULL,
            return_INLA_result = TRUE,
            avg_sessions = TRUE,
            session_names = NULL,
            trim_INLA = TRUE
          )
        )
      # stop break for testing (bottom) ----
      # stop()
      if(class(result_obj) =="try-error") {
        cat(num_runs,"Run data",hem,"hemisphere","data",PW,"FAILED with message:\n",result_obj,"/n")
        num_errs <- num_errs + 1
      } else {
        if("avg" %in% names(result_obj$betas_classical)) {
          class_ests <- c(result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]])
          bayes_ests <- c(result_obj$betas_Bayesian$avg$data[[paste0("cortex_",hem)]])
        } else {
          class_ests <-
            unlist(sapply(result_obj$betas_classical, function(x) {
              out <- c(x$data[[paste0("cortex_", hem)]])
              return(out)
            }, simplify = F))
          bayes_ests <-
            unlist(sapply(result_obj$betas_Bayesian, function(x) {
              out <- c(x$data[[paste0("cortex_", hem)]])
              return(out)
            }, simplify = F))
        }
        class_corr <- cor(class_ests, true_coefs[[hem]])
        bayes_corr <- cor(bayes_ests, true_coefs[[hem]])
        cat(num_runs,"Run data",hem,"hemisphere","data",PW,"SUCCEEDED\n")
        cat("....Correlation w/ true:",class_corr,"(Classical)",
            bayes_corr,"(Bayesian)\n")
      }
    }
  }
}
cat("FINAL NUMBER OF ERRORS:",num_errs)
