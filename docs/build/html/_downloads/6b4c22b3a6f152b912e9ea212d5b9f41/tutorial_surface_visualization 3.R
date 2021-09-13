# R script
#
# Tutorial 1 - Surface Visualization 
# micapipe v0.1.1
#
# Created by RRC on September 2021 (the second year of the pademic)

# Set the environment
library("RColorBrewer")
library("viridis")
library('fsbrain')

# Set the working directory to the out directory
setwd("~/tmp/micaConn/micapipe_tutorials") # <<<<<<<<<<<< CHANGE THIS PATH

# This variable will be different for each subject
subjectID <- 'sub-HC001_ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
subjectDir <- 'micapipe/sub-HC001/ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's DIRECTORY

# Here we define the atlas 
atlas <- 'schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS

# Set paths and variables
dir_FS <- paste0('freesurfer/' + subjectID)
dir_conte <- paste0(subjectDir + '/anat/surfaces/conte69/')
dir_morph <- paste0(subjectDir + '/anat/surfaces/morphology/')
dir_mpc <- paste0(subjectDir + '/anat/surfaces/micro_profiles/')
dir_func <- paste0(subjectDir + '/func/surfaces/')

# ------------------------------------------------------------------------ # 
#### Morphology ####
# Thickness
# Thickness: Native surface
# Native surface
nom <-  paste0(BIDSid,"_space-fsnative_desc-qc_morphology-thickness.png"); nom.png <- paste0(dir_QC_png,"/",nom)
th.lh <- paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-lh_thickness.mgh")
th.rh <- paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-rh_thickness.mgh")
if (file.exists(th.lh) && file.exists(th.rh) ) {
  if (!file.exists(nom.png)) {
    print(paste0("[INFO].... Creating ",nom, " of the THICKNESS on native surface"))
    rgla <- list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T);
    cm <- vis.data.on.subject(dir_fs, BIDSid, morph_data_lh=th.lh, morph_data_rh=th.rh,
                              draw_colorbar = TRUE, surface="pial", rglactions = rgla,  makecmap_options = list('colFn'=inferno))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Thickness [mm]', grid_like = FALSE, view_angles = view_angles);
    while (rgl.cur() > 0) { rgl.close() }
  } else { print(paste0("[INFO].... File exists: ", nom)) }
}



# Thickness: fsaverage5
# fsaverage5 thickness
nom <-  paste0(BIDSid,"_space-fsaverage5_desc-qc_morphology-thickness.png"); nom.png <- paste0(dir_QC_png,"/",nom)
th.lh.fs5 <- paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-lh_thickness.mgh")
th.rh.fs5 <- paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-rh_thickness.mgh")
if (file.exists(th.lh.fs5) && file.exists(th.rh.fs5)) {
  if (!file.exists(nom.png)) {
    print(paste0("[INFO].... Creating ",nom, " of the THICKNESS on fsaverage5"))
    rgla <- list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T);
    cm <- vis.data.on.subject(dir_fs, 'fsaverage5', morph_data_lh=th.lh.fs5, morph_data_rh=th.rh.fs5,
                              draw_colorbar = TRUE, surface="pial", rglactions = rgla,  makecmap_options = list('colFn'=inferno))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Thickness [mm]', grid_like = FALSE, view_angles = view_angles);
    while (rgl.cur() > 0) { rgl.close() }
  } else { print(paste0("[INFO].... File exists: ", nom)) }
}


# Thickness: conte 69



# Curvature
# Curvature: Native surface
####  Morphology - Curvature ####
nom <-  paste0(BIDSid,"_space-fsnative_desc-qc_morphology-curvature.png"); nom.png <- paste0(dir_QC_png,"/",nom)
cv.lh <- paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-lh_curvature.mgh")
cv.rh <- paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-rh_curvature.mgh")
if (file.exists(cv.lh) && file.exists(cv.rh)) {
  if (!file.exists(nom.png)) {
    print(paste0("[INFO].... Creating PNG of the CURVATURE on native surface"))
    rgla <- list('trans_fun'=limit_fun(-0.3, 0.3), 'no_vis'=T);
    cm <- vis.data.on.subject(dir_fs, BIDSid, morph_data_lh=cv.lh, morph_data_rh=cv.rh,
                              draw_colorbar = 'horizontal', surface="inflated",  rglactions=rgla, makecmap_options = list('colFn'=col.curv))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Curvature', grid_like = FALSE, view_angles = view_angles);
    while (rgl.cur() > 0) { rgl.close() }
  } else { print(paste0("[INFO].... File exists: ", nom)) }
}



# Curvature: fsaverage5
# fsaverage5 curvature
nom <-  paste0(BIDSid,"_space-fsaverage5_desc-qc_morphology-curvature.png"); nom.png <- paste0(dir_QC_png,"/",nom)
cv.lh.fs5 <- paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-lh_curvature.mgh")
cv.rh.fs5 <- paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-rh_curvature.mgh")
if (file.exists(cv.lh.fs5) && file.exists(cv.rh.fs5)) {
  if (!file.exists(nom.png)) {
    print(paste0("[INFO].... Creating Creating PNG of the CURVATURE on fsaverage5"))
    rgla <- list('trans_fun'=limit_fun(-0.3, 0.3), 'no_vis'=T);
    cm <- vis.data.on.subject(dir_fs, 'fsaverage5', morph_data_lh=cv.lh.fs5, morph_data_rh=cv.rh.fs5,
                              draw_colorbar = TRUE, surface="inflated", rglactions = rgla,  makecmap_options = list('colFn'=col.curv))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Curvature', grid_like = FALSE, view_angles = view_angles);
    while (rgl.cur() > 0) { rgl.close() }
  } else { print(paste0("[INFO].... File exists: ", nom)) }
}


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
#### Microstructural profile covariance ####
mask.lh <- ifelse(freesurferformats::read.fs.morph(paste0(dir_fs,"/",BIDSid,"/surf/lh.thickness"))<0.5,0,1)
mask.rh <- ifelse(freesurferformats::read.fs.morph(paste0(dir_fs,"/",BIDSid,"/surf/rh.thickness"))<0.5,0,1)
for (i in 1:14) {
  i0 <- sprintf("%02d", i)
  mpc.lh <- paste0(dir_mpc,"/",BIDSid,"_space-fsnative_desc-lh_MPC-",i,".mgh")
  mpc.rh <- paste0(dir_mpc,"/",BIDSid,"_space-fsnative_desc-rh_MPC-",i,".mgh")
  if (file.exists(mpc.lh) && file.exists(mpc.rh)) {
    nom <-  paste0(BIDSid,"_space-fsnative_desc-qc_MPC-",i0,".png"); nom.png <- paste0(dir_QC_png,"/",nom)
    if (!file.exists(nom.png)) {
      # Load the data
      morph <- list(lh=freesurferformats::read.fs.morph(mpc.lh)*mask.lh,
                    rh=freesurferformats::read.fs.morph(mpc.rh)*mask.rh )
      # Calculate color range (quatile 20% to 95%) >>>> rgla <- list('trans_fun'=limit_fun(1200,2000), 'no_vis'=T);
      Qt <- round(quantile(c(morph$lh[morph$lh!=0], morph$rh[morph$rh!=0]), probs = c(0.05,0.95)),1)
      print(paste0("[INFO].... Creating PNG of MPC-", i0))
      rgla <- list('trans_fun'=limit_fun(Qt[1],Qt[2]), 'no_vis'=T);
      cm <- vis.data.on.subject(dir_fs, BIDSid, morph_data_lh=morph$lh, morph_data_rh=morph$rh,
                                surface = "pial", draw_colorbar = 'horizontal', rglactions = rgla)
      img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = paste0('MPC-',i0), grid_like = FALSE, view_angles = view_angles);
      while (rgl.cur() > 0) { rgl.close() }
    } else { print(paste0("[INFO].... File exists: ", nom)) }
  }
}


# ------------------------------------------------------------------------ # 
####  Schaefer-400 labels ####
# Parcellations surfaces
atlas <- gsub("lh.", "",  list.files(path =dir_fs.label, pattern = "lh.*_mics.annot"))
atlas <- gsub(".annot", "", atlas)
for (annot in atlas) {
  nom <- paste0(BIDSid,"_atlas-",annot,"_desc-qc.png"); nom.png <- paste0(dir_QC_png,"/",nom)
  if (!file.exists(nom.png)) {
    print(paste0("[INFO].... Creating PNG of ", annot, " on native surface"))
    cm <- vis.subject.annot(dir_fs, BIDSid, annot, 'both', surface="pial", rglactions= list('no_vis'=T))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, img_only = TRUE, grid_like = FALSE, view_angles = view_angles, colorbar_legend = annot);
    while (rgl.cur() > 0) { rgl.close() }
  } else {
    print(paste0("[INFO].... File exists: ", nom))
  }
}
if (file.exists("fsbrain_cbar.png")==TRUE) {file.remove("fsbrain_cbar.png")}
