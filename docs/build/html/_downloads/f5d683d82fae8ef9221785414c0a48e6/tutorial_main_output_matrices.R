# R script
#
# Tutorial 0 - Main output matrices 
# micapipe v0.1.1
# R version 3.6.3
#
# Created by RRC on September 2021 (the second year of the pademic)

# Set the environment
require('RColorBrewer')      # version 1.1-2
require('viridis')           # version 0.5.1

# Set the working directory to the out directory
setwd("~/tmp/micaConn/micapipe_tutorials") # <<<<<<<<<<<< CHANGE THIS PATH

# This variable will be different for each subject
subjectID <- 'sub-HC001_ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
subjectDir <- 'micapipe/sub-HC001/ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's DIRECTORY

# Here we define the atlas 
atlas <- 'schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS


# ------------------------------------------------------------------------ # 
#### Structural connectomes ####

### Full structural connectome
# Set the path to the the structural cortical connectome
cnt_sc_cor <- paste0(subjectDir, '/dwi/connectomes/', subjectID, '_space-dwi_atlas-', atlas, '_desc-iFOD2-40M-SIFT2_full-connectome.txt')

# Load the cortical connectome
mtx_sc <- as.matrix(read.csv(cnt_sc_cor, sep=" ", header=FALSE))

# Fill the lower triangle of the matrix
mtx_sc[lower.tri(mtx_sc)] <- t(mtx_sc)[lower.tri(mtx_sc)]

# Plot the log matrix
image(log(mtx_sc), axes=FALSE, main=paste0("SC ", atlas), col=brewer.pal(9, "Purples") )


###  Full structural connectome edge lengths
# Set the path to the the structural cortical connectome
cnt_sc_EL <- paste0(subjectDir, '/dwi/connectomes/', subjectID, '_space-dwi_atlas-', atlas, '_desc-iFOD2-40M-SIFT2_full-edgeLengths.txt')

# Load the cortical connectome
mtx_scEL <- as.matrix(read.csv(cnt_sc_EL, sep=" ", header=FALSE,))

# Fill the lower triangle of the matrix
mtx_scEL[lower.tri(mtx_scEL)] <- t(mtx_scEL)[lower.tri(mtx_scEL)]

# Plot the log matrix
image(log(mtx_scEL), axes=FALSE, main=paste0("SC ", atlas), col=brewer.pal(9, "Purples"))


# ------------------------------------------------------------------------ # 
#### Functional connectomes ####

# Set the path to the the functional cortical connectome
cnt_fs <- paste0(subjectDir, '/func/surfaces/', subjectID, '_rsfmri_space-fsnative_atlas-', atlas, '_desc-FC.txt')

# Load the cortical connectome
mtx_fs <- as.matrix(read.csv(cnt_fs, sep=" ", header=FALSE))

# Fill the lower triangle of the matrix
mtx_fs[lower.tri(mtx_fs)] <- t(mtx_fs)[lower.tri(mtx_fs)]

# Plot the matrix
image(mtx_fs, axes=FALSE, main=paste0("FC ", atlas), col=brewer.pal(9, "Reds"))

### Time series
# Set the path to the the time series file
cnt_time <- paste0(subjectDir, '/func/surfaces/', subjectID, '_rsfmri_space-fsnative_atlas-', atlas, '_desc-timeseries.txt')

# Load the time series
mtx_time <- as.matrix(read.csv(cnt_time, sep=" ", header=FALSE))

# Plot as a matrix
image(mtx_time, axes=FALSE, main=paste0("Time series ", atlas), col=plasma(64))

# ------------------------------------------------------------------------ # 
#### MPC connectomes ####

# Set the path to the the MPC cortical connectome
cnt_mpc <- paste0(subjectDir, '/anat/surfaces/micro_profiles/', subjectID, '_space-fsnative_atlas-', atlas, '_desc-MPC.txt')

# Load the cortical connectome
mtx_mpc <- as.matrix(read.csv(cnt_mpc, sep=" ", header=FALSE))

# Fill the lower triangle of the matrix
mtx_mpc[lower.tri(mtx_mpc)] <- t(mtx_mpc)[lower.tri(mtx_mpc)]

# Plot the matrix
image(mtx_mpc, axes=FALSE, main=paste0("MPC ", atlas), col=brewer.pal(9, "Greens"))


### Intensity profiles (Profile x ROI)
# Set the path to the the time series file
cnt_int <- paste0(subjectDir, '/anat/surfaces/micro_profiles/', subjectID, '_space-fsnative_atlas-', atlas, '_desc-intensity_profiles.txt')

# Load the time series
mtx_int <- as.matrix(read.csv(cnt_int, sep=" ", header=FALSE))

# Plot as a matrix
image(mtx_int, axes=FALSE, main=paste0("Intensity profiles ", atlas), col=brewer.pal(9, "Greens"))

# ------------------------------------------------------------------------ # 
#### Geodesic distance connectomes ####

# Set the path to the the geodesic distance connectome
cnt_gd <- paste0(subjectDir, '/anat/surfaces/geo_dist/', subjectID, '_space-fsnative_atlas-', atlas, '_GD.txt')

# Load the cortical connectome
mtx_gd <- as.matrix(read.csv(cnt_gd, sep=" ", header=FALSE))

# Plot the matrix
image(mtx_gd, axes=FALSE, main=paste0("GD ", atlas), col=brewer.pal(9, "Blues"))
