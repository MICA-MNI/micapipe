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
dir_surf <- paste0(subjectDir, '/surf/')
dir_maps <- paste0(subjectDir, '/maps/')

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
th.lh <- paste0(dir_maps, subjectID, '_space-fsnative_desc-lh_thickness.mgh')
th.rh <- paste0(dir_maps, subjectID, '_space-fsnative_desc-rh_thickness.mgh')

# Plot the surface
th_nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=th.lh, morph_data_rh=th.rh, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_nat, 'Thickness [mm]')


####  Thickness: fsaverage5
# fsaverage5 thickness
# Set the path to the surface
th.lh.fs5 <- paste0(dir_maps, subjectID, '_space-fsaverage5_desc-lh_thickness.mgh')
th.rh.fs5 <- paste0(dir_maps, subjectID, '_space-fsaverage5_desc-rh_thickness.mgh')

# Plot the surface
th_fs5 <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=th.lh.fs5, morph_data_rh=th.rh.fs5, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_fs5, 'Thickness [mm]')

####  Thickness: fsLR-32k
# Set the path to the surface
th.lh.f32k <- paste0(dir_maps, subjectID, 'hemi-L_surf-fsLR-32k_thickness.mgh')
th.rh.f32k <- paste0(dir_maps, subjectID, 'hemi-R_surf-fsLR-32k_thickness.mgh')

# Plot the surface
th_f32k <- vis.data.on.subject('freesurfer/', 'conte69', morph_data_lh=th.lh.f32k, morph_data_rh=th.rh.f32k, surface='conte69.gii', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))
plot_surface(th_f32k, 'Thickness [mm]')

####  Morphology - Curvature ####
# Colormap
RdYlGn <- colorRampPalette(brewer.pal(11,'RdYlGn'))

####  Curvature: Native surface
# Set the path to the surface
cv.lh <- paste0(dir_maps, subjectID, '_space-fsnative_desc-lh_curvature.mgh')
cv.rh <- paste0(dir_maps, subjectID, '_space-fsnative_desc-rh_curvature.mgh')

# Plot the surface
cv_nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=cv.lh, morph_data_rh=cv.rh, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_nat, 'Curvature [1/mm]')


####  Curvature: fsaverage5
# Set the path to the surface
cv.lh.fs5 <- paste0(dir_maps, subjectID, '_space-fsaverage5_desc-lh_curvature.mgh')
cv.rh.fs5 <- paste0(dir_maps, subjectID, '_space-fsaverage5_desc-rh_curvature.mgh')

# Plot the surface
cv_fs5 <- vis.data.on.subject('freesurfer/', 'fsaverage5', morph_data_lh=cv.lh.fs5, morph_data_rh=cv.rh.fs5, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_fs5, 'Curvature [1/mm]')


#### Curvature: fsLR-32k
# Set the path to the surface
cv.lh.f32k <- paste0(dir_maps, subjectID, 'hemi-L_surf-fsLR-32k_curvature.mgh')
cv.rh.f32k <- paste0(dir_maps, subjectID, 'hemi-R_surf-fsLR-32k_curvature.mgh')

# Plot the surface
cv_f32k <- vis.data.on.subject('freesurfer', 'conte69', morph_data_lh=cv.lh.f32k, morph_data_rh=cv.rh.f32k, surface='conte69.gii', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))
plot_surface(cv_f32k, 'Curvature [1/mm]')

# ------------------------------------------------------------------------ #
#### fsLR-32k ####
### fsLR-32k: Pial surface ###
# Colormap
grays <- colorRampPalette(c('gray65', 'gray65', 'gray65'))

# Set the path to the surface
f32k.pial.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-L_surf-fsLR-32k_pial.surf.gii') )
f32k.pial.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-R_surf-fsLR-32k_pial.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(f32k.pial.lh, morph_data = rnorm(nrow(f32k.pial.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(f32k.pial.rh, morph_data = rnorm(nrow(f32k.pial.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.pial <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.pial, 'conte69 pial')

### fsLR-32k: Middle surface  ###
# Set the path to the surface
f32k.mid.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-L_surf-fsLR-32k_midthickness.surf.gii') )
f32k.mid.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-R_surf-fsLR-32k_midthickness.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(f32k.mid.lh, morph_data = rnorm(nrow(f32k.mid.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(f32k.mid.rh, morph_data = rnorm(nrow(f32k.mid.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.mid <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.mid, 'conte69 mid')

### fsLR-32k: White surface  ###
# Set the path to the surface
f32k.wm.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-L_surf-fsLR-32k_white.surf.gii') )
f32k.wm.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'hemi-R_surf-fsLR-32k_white.surf.gii') )

# Plot the surface
cml = coloredmesh.from.preloaded.data(f32k.wm.lh, morph_data = rnorm(nrow(f32k.wm.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(f32k.wm.rh, morph_data = rnorm(nrow(f32k.wm.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.wm <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.wm, 'conte69 white')

# ------------------------------------------------------------------------ #
#### Native sphere ####
# Colormap
grays <- colorRampPalette(c('white', 'gray65','black'))

# Set the path to the surface
sph.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_lh_sphereReg.surf.gii'))
sph.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_rh_sphereReg.surf.gii'))

# Set the color limits
lf= limit_fun(-0.2, 0.2)

# Create the coloredmeshes
cml = coloredmesh.from.preloaded.data(sph.lh, morph_data = lf(read.fs.mgh(cv.lh)), hemi = 'lh', makecmap_options = list('colFn'=grays))
cmr = coloredmesh.from.preloaded.data(sph.rh, morph_data = lf(read.fs.mgh(cv.rh)), hemi = 'rh', makecmap_options = list('colFn'=grays))
sph.nat <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), rglactions = list('no_vis'=T))

# Plot the surface
plot_surface(sph.nat, 'Native sphere curvature [1/mm]')

# ------------------------------------------------------------------------ #
#### Superficial White Matter (SWM) in fsnative surface ####
###  SWM 1,2,3mm
for (mm in 1:3) {
  # Set the path to the surface
  f32k.swm.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_surf-fsnative_label-swm',mm,'.0mm.surf.gii') )
  f32k.swm.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_surf-fsnative_label-swm',mm,'.0mm.surf.gii') )

  # Plot the surface
  cml = coloredmesh.from.preloaded.data(f32k.swm.lh, morph_data = rep(0, nrow(f32k.swm.lh$vertices)), makecmap_options = list('colFn'=grays) )
  cmr = coloredmesh.from.preloaded.data(f32k.swm.rh, morph_data = rep(0, nrow(f32k.swm.rh$vertices)), makecmap_options = list('colFn'=grays) )
  brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
             rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=F))
}


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
