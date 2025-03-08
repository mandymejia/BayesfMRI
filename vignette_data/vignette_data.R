library(ciftiTools)
ciftiTools.setOption('wb_path','~/Applications') #where your Connectome Workbench installation is located
library(BayesfMRI)

# Data:
# HCP subject 100307
# session 1 = tfMRI_EMOTION_LR
# session 2 = tfMRI_EMOTION_RL

setwd('~/Documents/Github/BayesfMRI/vignette_data/')

sessions <- c('tfMRI_EMOTION_LR','tfMRI_EMOTION_RL')
fnames0 <- paste0(sessions, '_Atlas_MSMAll.dtseries.nii')
fname_BOLD <- 'BOLD_10k.dtseries.nii'

subjects <- c('100206','100307','100408','100610','101006')
nS <- length(subjects)

# ---------------------------------------------------------------
# Resample BOLD data

for(ii in 1:nS){
  print(subj_ii <- subjects[ii])

  for(jj in 1:2){

    sess_jj <- sessions[jj]

    #read in BOLD
    fname0_jj <- file.path('original', subj_ii, sess_jj, fnames0[jj])
    BOLD <- read_cifti(fname0_jj, resamp_res = 10000, brainstructures = "all")

    #write resampled BOLD
    fname_jj <- file.path('original', subj_ii, sess_jj, fname_BOLD)
    write_cifti(BOLD, cifti_fname = fname_jj)
  }
}

# ---------------------------------------------------------------
# Run BayesGLM for each session

for(ii in 1:nS){

  print(subj_ii <- subjects[ii])

  for(jj in 1:2){

    print(sess_jj <- sessions[jj])

    #read in resampled BOLD
    fname_jj <- file.path('original', subj_ii, sess_jj, fname_BOLD)
    BOLD <- read_cifti(fname_jj, brainstructures = "all")
    TR <- BOLD$meta$cifti$time_step
    nT <- ncol(BOLD)

    #task design and motion regressors
    fname_motion <- file.path("original",subj_ii, sess_jj, 'Movement_Regressors.txt')
    fname_events <- file.path("original",subj_ii, sess_jj, c('EVs/fear.txt', 'EVs/neut.txt'))
    motion <- as.matrix(read.table(fname_motion, header=FALSE))
    events <- lapply(fname_events, read.table, header=FALSE)
    names(events) <- c('fear', 'neut')
    design <- make_design(EVs=events, nTime = nT, TR = TR)

    # #to avoid collinearity
    # if(ii==1 & jj==1) motion <- motion[,-6]
    #
    # motion_df <- as.data.frame(motion)
    # motion_df$time <- 1:nrow(motion)
    # motion_df <- reshape2::melt(motion_df, id.vars='time')
    # ggplot(motion_df, aes(x=time, y=value, group=variable, color=variable)) + geom_line()
    #
    #run BayesGLM
    bglm <- BayesGLM(BOLD=BOLD,
                     design=design$design,
                     brainstructures="all",
                     subROI=c('Amygdala-L','Amygdala-R','Hippocampus-L','Hippocampus-R'),
                     TR=TR,
                     nuisance=motion,
                     scale_BOLD='mean',
                     hpf=.01,
                     nbhd_order=1,
                     ar_order=3,
                     ar_smooth=0,
                     Bayes=TRUE,
                     verbose=1 ,
                     meanTol=1)
    saveRDS(bglm, file = file.path("glms_full", paste0('subj',ii,'_sess',jj,'.rds')))
    #bglm <- readRDS(file = file.path("glms_full", paste0('subj',ii,'_sess',jj,'.rds')))

    act <- activations(bglm, Bayes = TRUE, gamma = 0, alpha = 0.05, verbose = 0)
    saveRDS(act, file = file.path("act", paste0('subj',ii,'_sess',jj,'.rds')))

    act0 <- activations(bglm, Bayes = FALSE, correction = "FWER", gamma = 0, alpha = 0.05, verbose = 0)
    act00 <- activations(bglm, Bayes = FALSE, correction = "FDR", gamma = 0, alpha = 0.05, verbose = 0)
    saveRDS(act0, file = file.path("act", paste0('subj',ii,'_sess',jj,'_0.rds')))
    saveRDS(act00, file = file.path("act", paste0('subj',ii,'_sess',jj,'_00.rds')))

  }
}

# Run BayesGLM2 (precook due to size of glms)

fnames <- list.files("glms_full")

# Second-level analysis 1: Fear & Neutral (mean over subjects)
bglm2 <- BayesGLM2(file.path("glms_full", fnames), excursion_type = '>', num_cores = 20)
saveRDS(bglm2, file='bglm2.rds')

# Second-level analysis 2: Fear - Neutral (mean over subjects)
bglm2b <- BayesGLM2(file.path("glms_full", fnames), excursion_type = '>', num_cores = 20,
                    contrasts = list(fear_vs_neut = rep(c(1, -1), 10)/10))
saveRDS(bglm2b, file='bglm2b.rds')






