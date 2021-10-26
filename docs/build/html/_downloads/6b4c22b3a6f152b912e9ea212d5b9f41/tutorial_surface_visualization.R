# R script
#
# Tutorial 1 - Surface Visualization 
# micapipe v0.1.1
# R version 3.6.3
#
# Created by RRC on September 2021 (the second year of the pademic)

# ------------------------------------------------------------------------ # 
####  Set the environment #### 
require('RColorBrewer')      # version 1.1-2
require('viridis')           # version 0.5.1
require('fsbrain')           # version 0.4.2
require('freesurferformats') # version 0.1.14
require('rgl')               # version 0.100.54

# Set the working directory to the out directory
setwd('~/tmp/micaConn/micapipe_tutorials') # <<<<<<<<<<<< CHANGE THIS PATH

# This variable will be different for each subject
subjectID <- 'sub-HC001_ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
subjectDir <- 'micapipe/sub-HC001/ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's DIRECTORY

# Set paths and variables
dir_conte <- paste0(subjectDir, '/anat/surfaces/conte69/')
dir_morph <- paste0(subjectDir, '/anat/surfaces/morphology/')
dir_mpc <- paste0(subjectDir, '/anat/surfaces/micro_profiles/')

# Helper function
plot_surface <-function(brainMesh, legend='', view_angles=c('sd_lateral_lh', 'sd_medial_lh', 'sd_medial_rh', 'sd_lateral_rh'), img_only=FALSE, horizontal=TRUE) {
  try(img <- vis.export.from.coloredmeshes(brainMesh, colorbar_legend = legend, grid_like = FALSE, view_angles = view_angles, img_only = img_only, horizontal=horizontal))
  while (rgl.cur() > 0) { rgl.close() }; file.remove(list.files(path = getwd(), pattern = 'fsbrain'))
  return(img)
}

# ------------------------------------------------------------------------ # 
#### Morphology thickness ####
####  Thickness: Inflated native surface
# Set the path to the surface
th.lh <- paste0(dir_morph, subjectID, '_space-fsnative_desc-lh_thickness.mgh')
th.rh <- paste0(dir_morph, subjectID, '_space-fsnative_desc-rh_thickness.mgh')

# Plot the surface
th_nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=th.lh, morph_data_rh=th.rh, surface='inflated', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_nat, 'Thickness [mm]')


####  Thickness: fsaverage5
# fsaverage5 thickness
# Set the path to the surface
th.lh.fs5 <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-lh_thickness.mgh')
th.rh.fs5 <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-rh_thickness.mgh')

# Plot the surface
th_fs5 <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=th.lh.fs5, morph_data_rh=th.rh.fs5, surface='inflated', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_fs5, 'Thickness [mm]')

####  Thickness: conte 69
# Set the path to the surface
th.lh.c69 <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-lh_thickness.mgh')
th.rh.c69 <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-rh_thickness.mgh')

# Plot the surface
th_c69 <- vis.data.on.subject('freesurfer/', 'conte69', morph_data_lh=th.lh.c69, morph_data_rh=th.rh.c69, surface='conte69.gii', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_c69, 'Thickness [mm]')

####  Morphology - Curvature ####
# Colormap
RdYlGn <- colorRampPalette(brewer.pal(11,'RdYlGn'))

####  Curvature: Native surface
# Set the path to the surface
cv.lh <- paste0(dir_morph, subjectID, '_space-fsnative_desc-lh_curvature.mgh')
cv.rh <- paste0(dir_morph, subjectID, '_space-fsnative_desc-rh_curvature.mgh')

# Plot the surface
cv_nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=cv.lh, morph_data_rh=cv.rh, surface='inflated', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_nat, 'Curvature [1/mm]')


####  Curvature: fsaverage5
# Set the path to the surface
cv.lh.fs5 <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-lh_curvature.mgh')
cv.rh.fs5 <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-rh_curvature.mgh')

# Plot the surface
cv_fs5 <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=cv.lh.fs5, morph_data_rh=cv.rh.fs5, surface='inflated', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_fs5, 'Curvature [1/mm]')


#### Curvature: conte 69
# Set the path to the surface
cv.lh.c69 <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-lh_curvature.mgh')
cv.rh.c69 <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-rh_curvature.mgh')

# Plot the surface
cv_c69 <- vis.data.on.subject('freesurfer', 'conte69', morph_data_lh=cv.lh.c69, morph_data_rh=cv.rh.c69, surface='conte69.gii', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_c69, 'Curvature [1/mm]')


# ------------------------------------------------------------------------ # 
#### Smoothed ####
# Thickness fsaverage5 fwhm=10mm
# Set the path to the surface
th.lh.fs5.10mm <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-lh_thickness_10mm.mgh')
th.rh.fs5.10mm <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-rh_thickness_10mm.mgh')

# Plot the surface
th_fs5.10mm <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=th.lh.fs5.10mm, morph_data_rh=th.rh.fs5.10mm, surface='pial', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_fs5.10mm, 'Thickness [mm]')


### Thickness conte69 fwhm=10mm
# Set the path to the surface
th.lh.c69.10mm <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-lh_thickness_10mm.mgh')
th.rh.c69.10mm <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-rh_thickness_10mm.mgh')

# Plot the surface
th_c69.10mm <- vis.data.on.subject('freesurfer/', 'conte69', morph_data_lh=th.lh.c69.10mm, morph_data_rh=th.rh.c69.10mm, surface='conte69.gii', draw_colorbar = TRUE, 
                                   views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_c69.10mm, 'Thickness [mm]')


### Curvature fsaverage5 fwhm=10mm
# Set the path to the surface
cv.lh.fs5.10mm <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-lh_curvature_10mm.mgh')
cv.rh.fs5.10mm <- paste0(dir_morph, subjectID, '_space-fsaverage5_desc-rh_curvature_10mm.mgh')

# Plot the surface
cv_fs5.10mm <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=cv.lh.fs5.10mm, morph_data_rh=cv.rh.fs5.10mm, surface='pial', draw_colorbar = TRUE, 
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_fs5.10mm, 'Curvature [1/mm]')


# Curvature conte69 fwhm=10mm
# Set the path to the surface
cv.lh.c69.10mm <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-lh_curvature_10mm.mgh')
cv.rh.c69.10mm <- paste0(dir_morph, subjectID, '_space-conte69-32k_desc-rh_curvature_10mm.mgh')

# Plot the surface
cv_c69.10mm <- vis.data.on.subject('freesurfer', 'conte69', morph_data_lh=cv.lh.c69.10mm, morph_data_rh=cv.rh.c69.10mm, surface='conte69.gii', draw_colorbar = TRUE, 
                                   views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_c69.10mm, 'Curvature [1/mm]')


# ------------------------------------------------------------------------ # 
#### Conte 69 ####
### Conte 69: Pial surface ###
# Colormap
grays <- colorRampPalette(c('gray65', 'gray65', 'gray65'))

# Set the path to the surface
c69.pial.lh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-lh_pial.surf.gii') )
c69.pial.rh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-rh_pial.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(c69.pial.lh, morph_data = rnorm(nrow(c69.pial.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(c69.pial.rh, morph_data = rnorm(nrow(c69.pial.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
c69.pial <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(c69.pial, 'conte69 pial')

### Conte 69: Middle surface  ### 
# Set the path to the surface
c69.mid.lh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-lh_midthickness.surf.gii') )
c69.mid.rh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-rh_midthickness.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(c69.mid.lh, morph_data = rnorm(nrow(c69.mid.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(c69.mid.rh, morph_data = rnorm(nrow(c69.mid.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
c69.mid <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(c69.mid, 'conte69 mid')

### Conte 69: White surface  ### 
# Set the path to the surface
c69.wm.lh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-lh_white.surf.gii') )
c69.wm.rh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_space-conte69-32k_desc-rh_white.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(c69.wm.lh, morph_data = rnorm(nrow(c69.wm.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(c69.wm.rh, morph_data = rnorm(nrow(c69.wm.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
c69.wm <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(c69.wm, 'conte69 white')

# ------------------------------------------------------------------------ # 
#### Native sphere ####
# Colormap
grays <- colorRampPalette(c('white', 'gray65','black'))

# Set the path to the surface
sph.lh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_lh_sphereReg.surf.gii'))
sph.rh <- read.fs.surface(filepath = paste0(dir_conte, subjectID,'_rh_sphereReg.surf.gii'))

# Set the color limits
lf= limit_fun(-0.2, 0.2)

# Create the coloredmeshes
cml = coloredmesh.from.preloaded.data(sph.lh, morph_data = lf(read.fs.mgh(cv.lh)), hemi = 'lh', makecmap_options = list('colFn'=grays))
cmr = coloredmesh.from.preloaded.data(sph.rh, morph_data = lf(read.fs.mgh(cv.rh)), hemi = 'rh', makecmap_options = list('colFn'=grays))
sph.nat <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), rglactions = list('no_vis'=T))

# Plot the surface
plot_surface(sph.nat, 'Native sphere curvature [1/mm]')

# ------------------------------------------------------------------------ # 
#### Microstructural profile covariance ####
###  MPC native surface
# Create a mask
mask.lh <- ifelse(read.fs.morph(th.lh)<0.5,0,1)
mask.rh <- ifelse(read.fs.morph(th.rh)<0.5,0,1)

# Set the path to the surface
mpc.lh <- paste0(dir_mpc, subjectID, '_space-fsnative_desc-lh_MPC-10.mgh')
mpc.rh <- paste0(dir_mpc, subjectID, '_space-fsnative_desc-rh_MPC-10.mgh')

# Load the data
mpc <- list(lh=read.fs.morph(mpc.lh)*mask.lh,  rh=read.fs.morph(mpc.rh)*mask.rh )

# Set color range based on MPC distribution
Qt <- round(quantile(c(mpc$lh[mpc$lh!=0], mpc$rh[mpc$rh!=0]), probs = c(0.05,0.95)),1)

# Plot the surface
mpc.nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=mpc$lh, morph_data_rh=mpc$rh, surface='pial', draw_colorbar = TRUE,
                               views=NULL, rglactions = list('trans_fun'=limit_fun(Qt[1],Qt[2]), 'no_vis'=T),  makecmap_options = list('colFn'=viridis))
plot_surface(mpc.nat, 'MPC-10')


###  MPC fsaverage5 surface
# Create a mask
mask.lh.fs5 <- ifelse(read.fs.morph(th.lh.fs5)<0.5,0,1)
mask.rh.fs5 <- ifelse(read.fs.morph(th.rh.fs5)<0.5,0,1)

# Set the path to the surface
mpc.lh.fs5 <- paste0(dir_mpc, subjectID, '_space-fsaverage5_desc-lh_MPC-10.mgh')
mpc.rh.fs5 <- paste0(dir_mpc, subjectID, '_space-fsaverage5_desc-rh_MPC-10.mgh')

# Load the data
mpc.fs5 <- list(lh=read.fs.morph(mpc.lh.fs5)*mask.lh.fs5,  rh=read.fs.morph(mpc.rh.fs5)*mask.rh.fs5 )

# Plot the surface
mpc.fs5 <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=mpc.fs5$lh, morph_data_rh=mpc.fs5$rh, surface='pial', draw_colorbar = TRUE,
                               views=NULL, rglactions = list('trans_fun'=limit_fun(Qt[1],Qt[2]), 'no_vis'=T),  makecmap_options = list('colFn'=viridis))
plot_surface(mpc.fs5, 'MPC-10')


###  MPC conte69 surface
# Create a mask
mask.lh.c69 <- ifelse(read.fs.morph(th.lh.c69)<0.5,0,1)
mask.rh.c69 <- ifelse(read.fs.morph(th.rh.c69)<0.5,0,1)

# Set the path to the surface
mpc.lh.c69 <- paste0(dir_mpc, subjectID, '_space-conte69-32k_desc-lh_MPC-10.mgh')
mpc.rh.c69 <- paste0(dir_mpc, subjectID, '_space-conte69-32k_desc-rh_MPC-10.mgh')

# Load the data
mpc.c69 <- list(lh=read.fs.morph(mpc.lh.c69)*mask.lh.c69,  rh=read.fs.morph(mpc.rh.c69)*mask.rh.c69 )

# Plot the surface
mpc.c69 <- vis.data.on.subject('freesurfer/', 'conte69', morph_data_lh=mpc.c69$lh, morph_data_rh=mpc.c69$rh, surface='conte69.gii', draw_colorbar = TRUE,
                               views=NULL, rglactions = list('trans_fun'=limit_fun(Qt[1],Qt[2]), 'no_vis'=T),  makecmap_options = list('colFn'=viridis))
plot_surface(mpc.c69, 'MPC-10')


# ------------------------------------------------------------------------ # 
####  Schaefer-400 labels ####
# Plot the surface
schaefer.400 <- vis.subject.annot('freesurfer/', subjectID, 'schaefer-400_mics', 'both', surface='pial',
                           views=NULL, rglactions = list('no_vis'=T))
plot_surface(schaefer.400, 'Schaefer-400')

####  Economo labels ####
# Plot the surface
economo <- vis.subject.annot('freesurfer/', subjectID, 'economo_mics', 'both', surface='pial',
                           views=NULL, rglactions = list('no_vis'=T))
plot_surface(economo, 'economo', img_only=TRUE)
