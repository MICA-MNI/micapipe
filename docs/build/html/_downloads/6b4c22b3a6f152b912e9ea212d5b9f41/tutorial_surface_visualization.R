# R script
#
# Tutorial 1 - Surface Visualization 
# micapipe v0.1.1
#
# Created by RRC on September 2021 (the second year of the pademic)

# Set the working directory to your subjec's directory
setwd("~/out/micapipe/sub-HC001/ses-01") # <<<<<<<<<<<< CHANGE THIS PATH

# Load required packages
library("RColorBrewer")
library("viridis")

# This variable will be different for each subject
subjectID <- 'sub-HC001_ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID

# Here we define the atlas 
atlas <- 'schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS


# ------------------------------------------------------------------------ # 
#### Morphology ####
# Thickness
# Thickness: Native surface



# Thickness: fsaverage5



# Thickness: conte 69



# Curvature
# Curvature: Native surface



# Curvature: fsaverage5



# Curvature: conte 69



# ------------------------------------------------------------------------ # 
#### Smoothed ####
# Thickness fsaverage5 fwhm=10mm



# Thickness conte69 fwhm=10mm



# Curvature fsaverage5 fwhm=10mm



# Curvature conte69 fwhm=10mm



# ------------------------------------------------------------------------ # 
#### Conte 69 ####
# Conte 69: Pial surface



# Conte 69: Middle surface



# Conte 69: White surface



# ------------------------------------------------------------------------ # 
#### Native sphere ####



# ------------------------------------------------------------------------ # 
####  rsfMRI on surface ####



# ------------------------------------------------------------------------ # 
####  MPC native surface ####



# ------------------------------------------------------------------------ # 
####  Schaefer-400 labels ####

