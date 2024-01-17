# [Build --> Install and Restart]

# Setup ------------------------------------------------------------------------
# [Edit these]
doINLA <- TRUE
saveResults <- TRUE
overwriteResults <- TRUE
resamp_res <- 9000
my_wb <- "~/Desktop/workbench" # path to your Connectome Workbench

dir_results <- "tests/results_notInPackage"
thisResultName <- gsub(
  ".", "_",
  as.character(packageVersion("BayesfMRI")[[1]]), fixed=TRUE
)
dir_resultsThis <- file.path(dir_results, thisResultName, "sim")
if (!overwriteResults && dir.exists(dir_resultsThis)) { stop("Results exist already.") }
if (!dir.exists(dir_resultsThis)) { dir.create(dir_resultsThis) }

library(testthat)
if (doINLA) { library(INLA) }
library(fMRItools)
library(ciftiTools)
ciftiTools.setOption('wb_path', my_wb)
library(BayesfMRI)

# Simulate data ----------------------------------------------------------------
bsim <- readRDS(file.path(dir_results, "simData/bsim.rds"))
bsim2 <- readRDS(file.path(dir_results, "simData/bsim_wSpatialVar.rds"))

# BayesGLM ---------------------------------------------------------------------
bfmri_fname <- file.path(dir_resultsThis, "bfmri_1sess.rds")
if (!file.exists(bfmri_fname)) {
  bfmri <- lapply(bsim$simulated_cifti[seq(6)], BayesGLM_cifti,
    brainstructures = "both",
    design = bsim$design,
    Bayes = TRUE,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 3000,
    verbose = FALSE,
    return_INLA = "trimmed"
  )
  saveRDS(bfmri, bfmri_fname)
}

bfmri_fname <- file.path(dir_resultsThis, "bfmri_2sess.rds")
if (!file.exists(bfmri_fname)) {
  bfmri <- lapply(
    list(
      bsim$simulated_cifti[c(1,4)],
      bsim$simulated_cifti[c(2,5)],
      bsim$simulated_cifti[c(3,6)]
    ),
    BayesGLM_cifti,
    brainstructures = "left",
    design = list(bsim$design, bsim$design),
    Bayes = TRUE,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 3000,
    verbose = FALSE
    return_INLA = "full",
    meanTol = 1,
    varTol = 1,
    scale_design = FALSE,
    scale_BOLD = "none"
  )
  saveRDS(bfmri, bfmri_fname)
}

bfmri_fname <- file.path(dir_resultsThis, "bfmri_2sess_scrambledDesign.rds")
if (!file.exists(bfmri_fname)) {
  bd2 <- bsim$design[,c(2,1,3)]
  colnames(bd2) <- colnames(bd2)[c(2,1,3)]
  bfmri <- lapply(
    list(
      bsim$simulated_cifti[c(1,4)],
      bsim$simulated_cifti[c(2,5)],
      bsim$simulated_cifti[c(3,6)]
    ),
    BayesGLM_cifti,
    brainstructures = "left",
    design = list(bsim$design, bd2),
    Bayes = TRUE,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 1500,
    DCT=0,
    verbose = FALSE,
    return_INLA = "minimal",
    session_names = c("normal", "scrambled")
  )
  saveRDS(bfmri, bfmri_fname)
}

bfmri_fname <- file.path(dir_resultsThis, "bfmri_2sess_spatialVar.rds")
if (!file.exists(bfmri_fname)) {
  bfmri <- lapply(
    list(
      bsim2$simulated_cifti[c(1,4)],
      bsim2$simulated_cifti[c(2,5)],
      bsim2$simulated_cifti[c(3,6)]
    ),
    BayesGLM_cifti,
    brainstructures = "left",
    design = list(bsim2$design, bsim2$design),
    Bayes = TRUE,
    ar_order = 6,
    ar_smooth = 5,
    resamp_res = 3000,
    verbose = FALSE,
    return_INLA = "trimmed"
  )
  saveRDS(bfmri, bfmri_fname)
}

# BayesGLM2 --------------------------------------------------------------------
nS <- 1; nN <- 6

### bfmri_1sess ----
bglm_fname <- file.path(dir_resultsThis, "bglm2_1sess.rds")
if (!file.exists(bglm_fname)) {
  contrasts <- list(
    field1 = rep(rep(c(1, 0, 0), nS), nN),
    field2 = rep(rep(c(0, 1, 0), nS), nN),
    field3 = rep(rep(c(0, 0, 1), nS), nN),
    f12diff = rep(rep(c(1, -1, 0), nS), nN),                 # first field, minus the second: w/ != contrast, expect union
    field1_m7 = rep(rep(c(1, 0, 0), nS), nN)                 # first field, gonna mult by 7: same as `field1`
  )
  contrasts <- lapply(contrasts, function(q){ q / sum(abs(q)) })
  contrasts$field1_m7 <- contrasts$field1_m7 * 7
  b2 <- list()
  b2$a <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_1sess.rds")),
    nsamp_theta=5, nsamp_beta=10
  )
  b2$b <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_1sess.rds")),
    contrasts=contrasts, excursion_type = ">", quantiles=.7,
    nsamp_theta=10, nsamp_beta=20
  )
  b2$c <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_1sess.rds")),
    contrasts=contrasts[rev(seq(length(contrasts)))],
    excursion_type = c(">", "!=", "<", "<", "<"), quantiles=.001,
    alpha=.01,
    nsamp_theta=3, nsamp_beta=5
  )
  saveRDS(b2, bglm_fname)
}
b2 <- readRDS(bglm_fname)
if (FALSE) {
  plot(b2$a)
  plot(b2$b)
  plot(b2$c)
}
plot(
  b2$c, what="act",
  fname=file.path(dir_resultsThis, paste0("bglm2_1sess.png")),
  together="idx", together_ncol=length(contrasts)
)

##### umm -----
b2$d <- BayesGLM2(
  readRDS(file.path(dir_resultsThis, "bfmri_1sess.rds")),
  contrasts=contrasts, excursion_type = ">", quantiles=.7,
  nsamp_theta=10, nsamp_beta=20, gamma=1, alpha=.05
)

nS <- 2; nN <- 3
### bfmri_2sess ----
contrasts <- list(
  field1 = rep(rep(c(1, 0, 0), nS), nN),
  field2 = rep(rep(c(0, 1, 0), nS), nN),
  field3 = rep(rep(c(0, 0, 1), nS), nN),
  f12diff = rep(rep(c(1, -1, 0), nS), nN),                 # first field, minus the second: w/ != contrast, expect union
  fl1ss2 = rep(c(c(0, 0, 0), rep(c(1, 0, 0), nS-1)), nN),  # first field, second session only: same as `field1`
  fl1sd = rep(c(c(1, 0, 0), rep(c(-1, 0, 0), nS-1)), nN),  # first field, first sesh minus second: empty
  field1_m7 = rep(rep(c(1, 0, 0), nS), nN)                 # first field, gonna mult by 7: same as `field1`
)
contrasts <- lapply(contrasts, function(q){ q / sum(abs(q)) })
contrasts$field1_m7 <- contrasts$field1_m7 * 7
bglm_fname <- file.path(dir_resultsThis, "bglm2_2sess.rds")
if (!file.exists(bglm_fname)) {
  b2 <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_2sess.rds")),
    contrasts <- contrasts, excursion_type = "!=", quantiles=.2,
    nsamp_theta=30, nsamp_beta=70
  )
  saveRDS(b2, bglm_fname)
}
b2 <- readRDS(bglm_fname)
plot(
  b2, what="act",
  fname=file.path(dir_resultsThis, paste0("bglm2_2sess.png")),
  together="idx", together_ncol=length(contrasts)
)

### bfmri_2sess_scrambledDesign ----
contrasts <- list(
  field1 = rep(rep(c(1, 0, 0), nS), nN),                   # :fields 1 and 2 (b/c scrambled)
  fl1ss2 = rep(c(c(0, 0, 0), c(1, 0, 0)), nN)              # :second field (b/c scrambled)
)
contrasts <- lapply(contrasts, function(q){ q / sum(abs(q)) })
bglm_fname <- file.path(dir_resultsThis, "bglm2_2sess_scrambledDesign.rds")
if (!file.exists(bglm_fname)) {
  b2 <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_2sess_scrambledDesign.rds")),
    contrasts=contrasts, excursion_type=">",
    nsamp_theta=30, nsamp_beta=70
  )
  saveRDS(b2, bglm_fname)
}
b2 <- readRDS(bglm_fname)
plot(
  b2, what="act",
  fname=file.path(dir_resultsThis, paste0("bglm2_2sess_scrambledDesign.png")),
  together="idx", together_ncol=length(contrasts)
)

### bfmri_2sess_spatialVar ----
contrasts <- list(
  field1 = rep(rep(c(1, 0, 0), nS), nN),
  field2 = rep(rep(c(0, 1, 0), nS), nN),
  field3 = rep(rep(c(0, 0, 1), nS), nN),
  f12diff = rep(rep(c(1, -1, 0), nS), nN),                 # first field, minus the second
  fl1ss1 = rep(c(c(1, 0, 0), rep(c(0, 0, 0), nS-1)), nN),  # first field, first session only
  fl1ss2 = rep(c(c(0, 0, 0), rep(c(1, 0, 0), nS-1)), nN),  # first field, second session only
  fl3sb2 = c(rep(0, 6), rep(c(0, 0, 1), 2), rep(0, 6)),    # third field, second subject only
  field2_sdiff = c(c(0,1,0, 0,-1,0), rep(0, 6*(nN-1)))     # second field, first sesh minus second, subj1
)
contrasts <- lapply(contrasts, function(q){ q / sum(abs(q)) })
bglm_fname <- file.path(dir_resultsThis, "bglm2_2sess_spatialVar.rds")
if (!file.exists(bglm_fname)) {
  b2 <- BayesGLM2(
    readRDS(file.path(dir_resultsThis, "bfmri_2sess_spatialVar.rds")),
    contrasts=contrasts, excursion_type=">",
    nsamp_theta=30, nsamp_beta=70
  )
  saveRDS(b2, bglm_fname)
}
b2 <- readRDS(bglm_fname)
plot(
  b2, what="act",
  fname=file.path(dir_resultsThis, paste0("bglm2_2sess_spatialVar.png")),
  together="idx", together_ncol=length(contrasts)
)
