# R script
#
# Tutorial 0 - Main output matrices
# micapipe v0.2.0
# R version 4.3.1
#
# Created by RRC in September 2021 (the second year of the pandemic)
# Updated by CJ in March 2024

# Set the environment
require('RColorBrewer')      # version 1.1-2
require('viridis')           # version 0.5.1
require('gifti')             # version 0.8.0

# Set the working directory to your subject's directory
setwd("out/micapipe_v0.2.0/sub-HC001/ses-01")

# This variable will be different for each subject
subjectID <- 'sub-HC001_ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID

# Here we define the atlas
atlas <- 'schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS


# ------------------------------------------------------------------------ #
#### Structural connectomes ####

### Full structural connectome
# Set the path to the the structural cortical connectome
cnt_sc_cor <- paste0('dwi/connectomes/', subjectID, '_space-dwi_atlas-', atlas, '_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii')

# Load the cortical connectome
mtx_sc <- readgii(cnt_sc_cor)$data$shape

# Fill the lower triangle of the matrix
mtx_sc[lower.tri(mtx_sc)] <- t(mtx_sc)[lower.tri(mtx_sc)]

# Plot the log matrix
image(log(mtx_sc), axes=FALSE, main=paste0("SC ", atlas), col=brewer.pal(9, "Purples") )


###  Full structural connectome edge lengths
# Set the path to the the structural cortical connectome
cnt_sc_EL <- paste0('dwi/connectomes/', subjectID, '_space-dwi_atlas-', atlas, '_desc-iFOD2-40M-SIFT2_full-edgeLengths.txt')

# Load the cortical connectome
mtx_scEL <- as.matrix(read.csv(cnt_sc_EL, sep=" ", header=FALSE,))

# Fill the lower triangle of the matrix
mtx_scEL[lower.tri(mtx_scEL)] <- t(mtx_scEL)[lower.tri(mtx_scEL)]

# Plot the log matrix
image(log(mtx_scEL), axes=FALSE, main=paste0("SC ", atlas), col=brewer.pal(9, "Purples"))


# ------------------------------------------------------------------------ #
#### Functional connectomes ####

# Set the path to the the functional cortical connectome
acq_func <- 'desc-se_task-rest_acq-AP_bold'
cnt_fs <- paste0('func/', acq_func, '/surf/', subjectID, '_surf-fsLR-32k_atlas-', atlas, '_desc-FC.shape.gii')

# Load the cortical connectome
mtx_fs <- readgii(cnt_fs)$data$shape

# Fill the lower triangle of the matrix
mtx_fs[lower.tri(mtx_fs)] <- t(mtx_fs)[lower.tri(mtx_fs)]

# Plot the matrix
image(mtx_fs, axes=FALSE, main=paste0("FC ", atlas), col=brewer.pal(9, "Reds"))

### Time series
# Set the path to the the time series file
cnt_time <- paste0('func/', acq_func, '/surf/', subjectID, '_surf-fsLR-32k_desc-timeseries_clean.shape.gii')

# Load the time series
mtx_time <- readgii(cnt_time)$data$shape

# Plot as a matrix
image(mtx_time, axes=FALSE, main=paste0("Time series ", atlas), col=plasma(64))

# ------------------------------------------------------------------------ #
#### MPC connectomes ####

# Set the path to the the MPC cortical connectome
cnt_mpc <- paste0('mpc/acq-T1map/', subjectID,'_atlas-', atlas, '_desc-MPC.shape.gii')

# Load the cortical connectome
mtx_mpc <- readgii(cnt_mpc)$data$shape

# Fill the lower triangle of the matrix
mtx_mpc[lower.tri(mtx_mpc)] <- t(mtx_mpc)[lower.tri(mtx_mpc)]

# Plot the matrix
image(mtx_mpc, axes=FALSE, main=paste0("MPC ", atlas), col=brewer.pal(9, "Greens"))

### Intensity profiles (Profile x ROI)
# Set the path to the the time series file
cnt_int <- paste0('mpc/acq-T1map/', subjectID,'_atlas-', atlas, '_desc-intensity_profiles.shape.gii')

# Load the time series
mtx_int <- readgii(cnt_int)$data$shape

# Plot as a matrix
image(t(mtx_int), axes=FALSE, main=paste0("Intensity profiles ", atlas), col=brewer.pal(9, "Greens"))

# ------------------------------------------------------------------------ #
#### Geodesic distance connectomes ####

# Set the path to the the geodesic distance connectome
cnt_gd <- paste0('dist/', subjectID, '_atlas-', atlas, '_GD.shape.gii')

# Load the cortical connectome
mtx_gd <- readgii(cnt_gd)$data$shape

# Plot the matrix
image(t(mtx_gd), axes=FALSE, main=paste0("GD ", atlas), col=brewer.pal(9, "Blues"))
