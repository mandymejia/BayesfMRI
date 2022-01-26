# This is a test script for the BayesfMRI package. If this runs successfully,
# that means that the package is working as intended.
# Cortical surface ----
library(brainSim)
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
library(INLA)
inla.setOption(pardiso.license = "~/licenses/pardiso.lic")
library(BayesfMRI)
save_path <- "~/Desktop"
set.seed(47408)
num_errs <- 0
for(num_runs in 1:2){
  cat("SIMULATING DATA FOR",num_runs,"RUN(S)\n")
  sim_data <-
    simulate_cifti_multiple(
      wb_path = "/Applications/workbench/",
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

  # Writing this cifti out is necessary in some versions of BayesfMRI
  cifti_fname <- sapply(seq(num_runs), function(rn) {
    if(rn == 1) {
      surfL_file = file.path(save_path,"surfL.gii")
      surfR_file = file.path(save_path,"surfR.gii")
    } else {
      surfL_file = NULL
      surfR_file = NULL
    }
    cf_nm <- file.path(save_path, paste0("run", rn, ".dtseries.nii"))
    write_cifti(
      xifti = sim_data$simulated_cifti[[rn]],
      cifti_fname = cf_nm,
      surfL_fname = surfL_file,
      surfR_fname = surfR_file,
      verbose = F
    )
    return(cf_nm)
  })
  for(hem in c("left","right")) {
    cat("....",toupper(hem),"HEMISPHERE\n")
    for(PW in c("not prewhitened","prewhitened")) {
      # stop break for testing ----
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
            cifti_fname = cifti_fname,
            surfL_fname = file.path(save_path, "surfL.gii"),
            surfR_fname = file.path(save_path, "surfR.gii"),
            brainstructures = hem,
            design = des,
            onsets = NULL,
            TR = 1,
            nuisance = NULL,
            nuisance_include = c("drift", "dHRF"),
            scale_BOLD = TRUE,
            scale_design = TRUE,
            # Bayes = TRUE, # This is the construction in BayesfMRI 1.9
            GLM_method = "all", # This is the construction for BayesfMRI 1.8.EM
            ar_order = ar_order,
            ar_smooth = ar_smooth,
            resamp_res = NULL,
            num.threads = 6,
            verbose = FALSE,
            outfile = NULL,
            return_INLA_result = TRUE,
            avg_sessions = TRUE,
            session_names = NULL,
            groups_df = NULL, # This is for BayesfMRI 1.8.EM
            tol = 1e-1, # This is for BayesfMRI 1.8.EM
            trim_INLA = TRUE
          )
        )
      if(class(result_obj) =="try-error") {
        cat(num_runs,"Run data",hem,"hemisphere","data",PW,"FAILED with message:\n",result_obj,"/n")
        num_errs <- num_errs + 1
      } else {
        if("avg" %in% names(result_obj$betas_classical)) {
          class_ests <- c(result_obj$betas_classical$avg$data[[paste0("cortex_",hem)]])
          bayes_ests <- c(result_obj$betas_Bayesian$avg$data[[paste0("cortex_",hem)]])
          em_ests <- c(result_obj$betas_EM$avg$data[[paste0("cortex_",hem)]])
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
          em_ests <-
            unlist(sapply(result_obj$betas_EM, function(x) {
              out <- c(x$data[[paste0("cortex_", hem)]])
              return(out)
            }, simplify = F))
        }
        class_corr <- cor(class_ests, true_coefs[[hem]])
        bayes_corr <- cor(bayes_ests, true_coefs[[hem]])
        em_corr <- cor(em_ests, true_coefs[[hem]])
        cat(num_runs,"Run data",hem,"hemisphere","data",PW,"SUCCEEDED\n")
        cat("....Correlation w/ true:",class_corr,"(Classical)",
            bayes_corr,"(Bayesian)", em_corr,"(EM)\n")
      }
    }
  }
}
cat("FINAL NUMBER OF ERRORS:",num_errs)

# Subcortex ----
# We don't have simulation set up for subcortical volumes (yet), so we use
# real data
library(ciftiTools)
ciftiTools.setOption("wb_path","/Applications/workbench")
library(BayesfMRI)
save_path <- "~/Desktop"
set.seed(47408)
num_errs <- 0
main_dir <- "/Volumes/GoogleDrive/My Drive/danspen/HCP_Motor_Task_Dan"
data_dir <- file.path(main_dir,"visit1_data") # Macbook
result_dir <- "~/Desktop"
TR = 0.72 #temporal resolution of data
thetas <- NULL # No starting values for precision parameters
regions <- c('Accumbens-l','Accumbens-r', #3,4 -- 200 voxels - BASAL GANGLIA --> MODEL 1
             'Amygdala-l','Amygdala-r',   #5,6 -- 600 voxels  IMPORTANT --> MODEL 2
             'Brain Stem',                #7 -- 3,402 voxels  EXCLUDE
             'Caudate-L','Caudate-R',     #8,9 -- 1400 voxels  BASAL GANGLIA --> MODEL 1
             'Cerebellum-L','Cerebellum-R', #10,11 -- 9000 voxels each --> MODELS 3 AND 4
             'Diencephalon-L','Diencephalon-R', #12,13 -- 1400 voxels  --> MODEL 2
             'Hippocampus-L','Hippocampus-R', #14,15 -- 1500 voxels  IMPORTANT --> MODEL 2
             'Pallidum-L','Pallidum-R', #16,17 -- 500 voxels  BASAL GANGLIA --> MODEL 1
             'Putamen-L','Putamen-R', #18,19 -- 2000 voxels  BASAL GANGLIA --> MODEL 1
             'Thalamus-L','Thalamus-R' #20,21 -- 2500 voxels IMPORTANT --> MODEL 2
)
groups <- c(1,1,rep(NA,17))
groups_df <- data.frame(label=3:21, region=regions, group=groups)
subject <- "103818"; visit <- 1; h <- 3
dir_s <- file.path(data_dir, subject, 'MNINonLinear', 'fsaverage_LR32k')
dir1_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_LR')
dir2_s <- file.path(data_dir,subject, 'MNINonLinear', 'Results', 'tfMRI_MOTOR_RL')
fname1_ts <- file.path(dir1_s,'tfMRI_MOTOR_LR_Atlas.dtseries.nii')
fname2_ts <- file.path(dir2_s,'tfMRI_MOTOR_RL_Atlas.dtseries.nii')
hem <- 'subcortical'
cols_h <- 5
tasks_h <- 'rh'
names_tasks_h <- "right_hand"
K <- 1
#Set up list of onsets
onsets1 <- onsets2 <- vector('list', length=K)
names(onsets1) <- names(onsets2) <- names_tasks_h
for(t in tasks_h){
  ind_t <- which(tasks_h==t)
  fname1_t <- file.path(dir1_s,'EVs',paste0(t,'.txt'))
  onsets1[[ind_t]] <- read.table(fname1_t, header=FALSE)
  fname2_t <- file.path(dir2_s,'EVs',paste0(t,'.txt'))
  onsets2[[ind_t]] <- read.table(fname2_t, header=FALSE)
}
#Set up nuisance regressors
motion1 <- as.matrix(read.table(file.path(dir1_s,'Movement_Regressors.txt'), header=FALSE))
motion2 <- as.matrix(read.table(file.path(dir2_s,'Movement_Regressors.txt'), header=FALSE))
use_fnames <- c(fname1_ts,fname2_ts)
start_time <- proc.time()[3]
result_svh <- BayesGLM_cifti(cifti_fname = use_fnames, # Multi-session
                             surfL_fname = NULL,
                             surfR_fname = NULL,
                             brainstructures = hem,
                             design = NULL,
                             onsets = list(onsets1, onsets2), # Multi-session
                             TR = TR,
                             nuisance = list(motion1, motion2), # Multi-session
                             nuisance_include = c('drift','dHRF'),
                             scale_BOLD = TRUE,
                             scale_design = TRUE,
                             GLM_method = 'EM',
                             ar_order = 0,
                             ar_smooth = 0,
                             session_names = c('LR','RL'), # Multiple sessions
                             resamp_res = NULL, # Don't forget to change this
                             num.threads = 6,
                             verbose = TRUE,
                             outfile = NULL,
                             return_INLA_result = T,
                             avg_sessions = T,
                             trim_INLA = T,
                             groups_df = groups_df,
                             tol = 1e-1)
total_time <- proc.time()[3] - start_time
result_svh$total_time <- total_time
