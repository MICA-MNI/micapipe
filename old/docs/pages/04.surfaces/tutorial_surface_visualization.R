# R script
#
# Tutorial 1 - Surface Visualization
# micapipe v0.2.0
# R version 4.3.1
#
# Created by RRC in September 2021 (the second year of the pandemic)
# Updated by CJ in March 2024

# ------------------------------------------------------------------------ #
####  Set the environment ####
require('RColorBrewer')      # version 1.1.3
require('viridis')           # version 0.6.5
require('fsbrain')           # version 0.5.5
require('freesurferformats') # version 0.1.18
require('rgl')               # version 1.2.8
require('gifti')             # version 0.8.0

# Define working directory and subject-specific information. CHANGE THESE PATHS as appropriate.
setwd('~/derivatives_RtD/')       # working directory
subjectID <- 'sub-HC001_ses-01'                  # subject ID
subjectDir <- 'micapipe_v0.2.0/sub-HC001/ses-01' # subject directory path

# Set paths for surface and morphometry data.
dir_surf <- paste0(subjectDir, '/surf/')
dir_maps <- paste0(subjectDir, '/maps/')

# Helper function to plot brain surfaces.
plot_surface <-function(brainMesh, legend='', view_angles=c('sd_lateral_lh', 'sd_medial_lh', 'sd_medial_rh', 'sd_lateral_rh'), img_only=FALSE) {
  try(img <- vis.export.from.coloredmeshes(brainMesh, colorbar_legend = legend, grid_like = FALSE, view_angles = view_angles, img_only = img_only, horizontal=TRUE))
  while (rgl.cur() > 0) { rgl.close() }; file.remove(list.files(path = getwd(), pattern = 'fsbrain'))
  return(img)
}

# Define color maps for brain visualization.
RdYlGn <- colorRampPalette(brewer.pal(11,'RdYlGn'))
bw <- colorRampPalette(c('black','gray65', 'white'))
grays <- colorRampPalette(c('gray65', 'gray65', 'gray65'))

# Set color limits for data visualization
lft = limit_fun(1.5, 4)    # thickness color scale
lfc = limit_fun(-0.2, 0.2) # curvature color scale

# ------------------------------------------------------------------------ #
#### Morphology thickness ####

##  Thickness: Inflated native surface
# Define paths for left and right hemisphere thickness data using FreeSurfer's native format (fsnative).
th.lh <- paste0(dir_maps, subjectID, "_hemi-L_surf-fsnative_label-thickness.func.gii")
th.rh <- paste0(dir_maps, subjectID, "_hemi-R_surf-fsnative_label-thickness.func.gii")

# Use 'vis.data.on.subject' to visualize thickness data on the inflated cortical surface, specific to FreeSurfer format.
th_nat <- vis.data.on.subject('freesurfer', subjectID, morph_data_lh=th.lh, morph_data_rh=th.rh, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T),  makecmap_options = list('colFn'=inferno))

# Display the native surface with mapped thickness data.
plot_surface(th_nat, 'Thickness [mm]')

##  Thickness: fsaverage5
# This section provides a pathway for visualizing thickness without requiring FreeSurfer-formatted surfaces, utilizing the fsaverage5 template.
# Load fsaverage5 surface paths for left and right hemispheres.
fs5.lh <-  read.fs.surface(filepath = 'micapipe/surfaces/fsaverage5/surf/lh.inflated')
fs5.rh <-  read.fs.surface(filepath = 'micapipe/surfaces/fsaverage5/surf/rh.inflated')

# Define paths for thickness morphometric data on fsaverage5 for each hemisphere.
fs5.lh.th <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsaverage5_label-thickness.func.gii')
fs5.rh.th <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsaverage5_label-thickness.func.gii')

# Create colored meshes for both hemispheres based on thickness data.
cml.fs5.th <- coloredmesh.from.preloaded.data(fs5.lh, morph_data = lft(unlist(readgii(fs5.lh.th)$data)), hemi = 'lh', makecmap_options = list('colFn'=inferno))
cmr.fs5.th <- coloredmesh.from.preloaded.data(fs5.rh, morph_data = lft(unlist(readgii(fs5.rh.th)$data)), hemi = 'rh', makecmap_options = list('colFn'=inferno))
th_fs5 <- brainviews(views = 't4', coloredmeshes=list('lh'=cml.fs5.th, 'rh'=cmr.fs5.th), rglactions = list('no_vis'=T))

# Display the fsaverage5 surface with thickness data.
plot_surface(th_fs5, 'Thickness [mm]')

##  Thickness: fsLR-32k
# Load fsLR-32k surface paths for left and right hemispheres.
f32k.lh <-  read.fs.surface(filepath = 'micapipe/surfaces/fsLR-32k.L.surf.gii', format = "gii")
f32k.rh <-  read.fs.surface(filepath = 'micapipe/surfaces/fsLR-32k.R.surf.gii', format = "gii")

# Define paths for thickness morphometry data on fsLR-32k for each hemisphere.
f32k.lh.th <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsLR-32k_label-thickness.func.gii')
f32k.rh.th <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsLR-32k_label-thickness.func.gii')

# Create colored meshes based on thickness data.
cml.f32k.th <- coloredmesh.from.preloaded.data(f32k.lh, morph_data = lft(unlist(readgii(f32k.lh.th)$data)), hemi = 'lh', makecmap_options = list('colFn'=inferno))
cmr.f32k.th <- coloredmesh.from.preloaded.data(f32k.rh, morph_data = lft(unlist(readgii(f32k.rh.th)$data)), hemi = 'rh', makecmap_options = list('colFn'=inferno))
th_f32k <- brainviews(views = 't4', coloredmeshes=list('lh'=cml.f32k.th, 'rh'=cmr.f32k.th), rglactions = list('no_vis'=T))

# Display the fsLR-32k surface with thickness data.
plot_surface(th_f32k, 'Thickness [mm]')

# ------------------------------------------------------------------------ #
####  Morphology - Curvature ####

##  Curvature: Native surface
# Define paths for left and right hemisphere curvature data in FreeSurfer's native format (fsnative).
cv.lh <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsnative_label-curv.func.gii')
cv.rh <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsnative_label-curv.func.gii')

# Use 'vis.data.on.subject' to visualize curvature on the subject's inflated surface.
cv_nat <- vis.data.on.subject('freesurfer/', subjectID, morph_data_lh=cv.lh, morph_data_rh=cv.rh, surface='inflated', draw_colorbar = TRUE,
                              views=NULL, rglactions = list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T),  makecmap_options = list('colFn'=RdYlGn))

# Display the native surface with mapped curvature data.
plot_surface(cv_nat, 'Curvature [1/mm]')

##  Curvature: fsaverage5
# Load fsaverage5 surface paths for left and right hemispheres.
fs5.lh <- read.fs.surface(filepath = 'micapipe/surfaces/fsaverage5/surf/lh.inflated')
fs5.rh <- read.fs.surface(filepath = 'micapipe/surfaces/fsaverage5/surf/rh.inflated')

# Define paths to the morphometry data on fsaverage5.
fs5.lh.cv <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsaverage5_label-curv.func.gii')
fs5.rh.cv <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsaverage5_label-curv.func.gii')

# Generate colored meshes using morphometry data.
cml_fs5 <- coloredmesh.from.preloaded.data(fs5.lh, morph_data = lfc(unlist(readgii(fs5.lh.cv)$data)), hemi = 'lh', makecmap_options = list('colFn'=RdYlGn))
cmr_fs5 <- coloredmesh.from.preloaded.data(fs5.rh, morph_data = lfc(unlist(readgii(fs5.rh.cv)$data)), hemi = 'rh', makecmap_options = list('colFn'=RdYlGn))
cv_fs5 <- brainviews(views = 't4', coloredmeshes=list('lh'=cml_fs5, 'rh'=cmr_fs5), rglactions = list('no_vis'=T))

# Display the fsaverage5 surface with morphometry data.
plot_surface(cv_fs5, 'Curvature [1/mm]')

## Curvature: fsLR-32k
# Load fsLF-32k surface paths for left and right hemispheres.
f32k.lh <- read.fs.surface(filepath = 'micapipe/surfaces/fsLR-32k.L.inflated.surf.gii')
f32k.rh <- read.fs.surface(filepath = 'micapipe/surfaces/fsLR-32k.R.inflated.surf.gii')

# Define paths to the morphometry data on fsLF32k.
f32k.lh.cv <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsLR-32k_label-curv.func.gii')
f32k.rh.cv <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsLR-32k_label-curv.func.gii')

# Generate colored meshes using morphometry data.
cml_f32k <- coloredmesh.from.preloaded.data(f32k.lh, morph_data = lfc(unlist(readgii(f32k.lh.cv)$data)), hemi = 'lh', makecmap_options = list('colFn'=RdYlGn))
cmr_f32k <- coloredmesh.from.preloaded.data(f32k.rh, morph_data = lfc(unlist(readgii(f32k.rh.cv)$data)), hemi = 'rh', makecmap_options = list('colFn'=RdYlGn))
cv.f32k <- brainviews(views = 't4', coloredmeshes=list('lh'=cml_f32k, 'rh'=cmr_f32k), rglactions = list('no_vis'=T))

# Display the fsLR-32k surface with morphometry data.
plot_surface(cv.f32k, 'Curvature [1/mm]')

## Curvature: fsLR-5k
# Load fsLF-5k surface paths for left and right hemispheres.
f5k.lh <- read.fs.surface(filepath = 'micapipe/surfaces/fsLR-5k.L.inflated.surf.gii')
f5k.rh <- read.fs.surface(filepath = 'micapipe/surfaces/fsLR-5k.R.inflated.surf.gii')

# Define paths to the morphometry data on fsLR-5k.
f5k.lh.cv <- paste0(dir_maps, subjectID, '_hemi-L_surf-fsLR-5k_label-curv.func.gii')
f5k.rh.cv <- paste0(dir_maps, subjectID, '_hemi-R_surf-fsLR-5k_label-curv.func.gii')

# Generate colored meshes using morphometry data.
cml_f5k <- coloredmesh.from.preloaded.data(f5k.lh, morph_data = lfc(unlist(readgii(f5k.lh.cv)$data)), hemi = 'lh', makecmap_options = list('colFn'=RdYlGn))
cmr_f5k <- coloredmesh.from.preloaded.data(f5k.rh, morph_data = lfc(unlist(readgii(f5k.rh.cv)$data)), hemi = 'rh', makecmap_options = list('colFn'=RdYlGn))
cv.f5k <- brainviews(views = 't4', coloredmeshes=list('lh'=cml_f5k, 'rh'=cmr_f5k), rglactions = list('no_vis'=T))

# Display the fsLR-5k surface with morphometry data.
plot_surface(cv.f5k, 'Curvature [1/mm]')

# ------------------------------------------------------------------------ #
#### fsLR-32k ####

## fsLR-32k: Pial surface
# Load fsLR-32k pial surfaces for left and right hemispheres.
f32k.pial.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii') )
f32k.pial.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii') )

# Prepare and plot the surfaces for visualization.
cml = coloredmesh.from.preloaded.data(f32k.pial.lh, morph_data = rnorm(nrow(f32k.pial.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(f32k.pial.rh, morph_data = rnorm(nrow(f32k.pial.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.pial <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.pial, 'fsLR-32k pial')

## fsLR-32k: Middle surface
# Load fsLR-32k middle surfaces for left and right hemispheres.
f32k.mid.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii') )
f32k.mid.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii') )

# Prepare and plot the surfaces for visualization.
cml = coloredmesh.from.preloaded.data(f32k.mid.lh, morph_data = rnorm(nrow(f32k.mid.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr = coloredmesh.from.preloaded.data(f32k.mid.rh, morph_data = rnorm(nrow(f32k.mid.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.mid <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.mid, 'fsLR-32k mid')

## fsLR-32k: White surface
# Load fsLR-32k white surfaces for left and right hemispheres.
f32k.wm.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii') )
f32k.wm.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii') )

# Prepare and plot the surfaces for visualization.
cml <- coloredmesh.from.preloaded.data(f32k.wm.lh, morph_data = rnorm(nrow(f32k.wm.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
cmr <- coloredmesh.from.preloaded.data(f32k.wm.rh, morph_data = rnorm(nrow(f32k.wm.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
f32k.wm <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
           rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
plot_surface(f32k.wm, 'fsLR-32k white')

# ------------------------------------------------------------------------ #
#### Native sphere ####

# Load native sphere surfaces for left and right hemispheres.
sph.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_surf-fsnative_label-sphere.surf.gii'))
sph.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_surf-fsnative_label-sphere.surf.gii'))

# Prepare colored meshes with curvature data.
cml <- coloredmesh.from.preloaded.data(sph.lh, morph_data = lfc(readgii(cv.lh)$data$unknown), hemi = 'lh', makecmap_options = list('colFn'=bw))
cmr <- coloredmesh.from.preloaded.data(sph.rh, morph_data = lfc(readgii(cv.rh)$data$unknown), hemi = 'rh', makecmap_options = list('colFn'=bw))
sph.nat <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), rglactions = list('no_vis'=T))

# Display the native sphere surface with mapped curvature data.
plot_surface(sph.nat, 'Native sphere curvature [1/mm]')

##--------------------------------------------------------------------- #
#### Superficial White Matter (SWM) in fsnative surface ####

# Define a function to visualize SWM at varying depths (1mm, 2mm, 3mm).
mesh_swm <- function(mm = '1') {
  
  # Load SWM surface data for specified depth for both hemispheres.
  f32k.swm.lh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-L_surf-fsnative_label-swm',mm,'.0mm.surf.gii') )
  f32k.swm.rh <- read.fs.surface(filepath = paste0(dir_surf, subjectID,'_hemi-R_surf-fsnative_label-swm',mm,'.0mm.surf.gii') )
  
  # Generate colored meshes using random morphometric data for visualization.
  cml = coloredmesh.from.preloaded.data(f32k.swm.lh, morph_data = rnorm(nrow(f32k.swm.lh$vertices),5,1), makecmap_options = list('colFn'=grays) )
  cmr = coloredmesh.from.preloaded.data(f32k.swm.rh, morph_data = rnorm(nrow(f32k.swm.rh$vertices),5,1), makecmap_options = list('colFn'=grays) )
  bv <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE,
             rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T))
  return(bv)
}

# Visualize SWM surfaces for depths of 1mm, 2mm, and 3mm.
plot_surface(mesh_swm('1'), 'Superficial white matter 1mm depth')
plot_surface(mesh_swm('2'), 'Superficial white matter 2mm depth')
plot_surface(mesh_swm('3'), 'Superficial white matter 3mm depth')

# ------------------------------------------------------------------------ #
####  qMRI Mapping functions  ####

# Function to aggregate qMRI surface data paths for a given type and surface.
surf_qmri <- function(qmri = '', surf = 'fsLR-32k') {
  
  # Compile file paths for qMRI data across hemispheres.
  files_lh <- paste0(dir_maps, subjectID, '_hemi-L_surf-', surf, '_label-midthickness_', qmri, '.func.gii')
  files_rh <- paste0(dir_maps, subjectID, '_hemi-R_surf-', surf, '_label-midthickness_', qmri, '.func.gii')
  
  # Combine left and right hemisphere data into a single array.
  surf_map <- c(readgii(files_lh)$data, readgii(files_rh)$data)
  
  return(surf_map)
}

# Function to visualize qMRI data on a specified brain surface.
mesh_qmri <- function(qmri = '', surf = 'fsLR-32k', label = 'pial', rq = c(0.15, 0.95), cmap = inferno) {
  
  # Retrieve qMRI data mapped to specified surface.
  map_surf <- surf_qmri(qmri, surf)
  
  # Load associated surface data for visualization.
  surf.lh <- read.fs.surface(filepath = paste0(dir_surf, '/', subjectID, '_hemi-L_space-nativepro_surf-', surf, '_label-', label, '.surf.gii' ) )
  surf.rh <- read.fs.surface(filepath = paste0(dir_surf, '/', subjectID, '_hemi-R_space-nativepro_surf-', surf, '_label-', label, '.surf.gii' ) )
  
  # Determine data range for color mapping based on specified quantiles.
  crange <- quantile(c(map_surf[[1]], map_surf[[2]]), probs = rq)
  
  # Set the color limits
  lf= limit_fun(crange[1], crange[2])
  
  # Apply color mapping based on the data range.
  cml = coloredmesh.from.preloaded.data(surf.lh, morph_data = lf(map_surf[[1]]), makecmap_options = list('colFn'=cmap) )
  cmr = coloredmesh.from.preloaded.data(surf.rh, morph_data = lf(map_surf[[2]]), makecmap_options = list('colFn'=cmap) )
  
  # Generate and return the visual representation.
  return(brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), draw_colorbar = FALSE, rglactions = list('trans_fun'=limit_fun(-1, 1), 'no_vis'=T)))
}

# Visualize T1 mapping on different standard and native surfaces.
plot_surface(mesh_qmri(qmri = 'T1map', surf = 'fsnative'))   # T1 map on fsnative surface
plot_surface(mesh_qmri(qmri = 'T1map', surf = 'fsaverage5')) # T1 map on fsaverage5 surface
plot_surface(mesh_qmri(qmri = 'T1map', surf = 'fsLR-32k'))   # T1 map on fsLR-32k surface
plot_surface(mesh_qmri(qmri = 'T1map', surf = 'fsLR-5k'))    # T1 map on fsLR-5k surface

# ------------------------------------------------------------------------ #
####  Schaefer-400 labels ####
# Generate and visualize Schaefer-400 labels on the pial surface.
schaefer.400 <- vis.subject.annot('freesurfer/', subjectID, 'schaefer-400_mics', 'both', surface='pial',
                           views=NULL, rglactions = list('no_vis'=T))
plot_surface(schaefer.400, 'Schaefer-400')

####  Economo labels ####
# Generate and visualize Economo cortical labels on the pial surface.
economo <- vis.subject.annot('freesurfer/', subjectID, 'economo_mics', 'both', surface='pial',
                           views=NULL, rglactions = list('no_vis'=T))
plot_surface(economo, 'economo', img_only=TRUE)