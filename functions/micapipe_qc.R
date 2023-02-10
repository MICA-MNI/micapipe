# R script
#
# Generates PNG images for micapipe Quality Check
# micapipe v0.1.5
# Extra functions
#
# Created by RRC on Feb-Apr 2021 (the second year of the pademic)
#
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      micapipe Quality Control

      Arguments:
      --BIDSid=<sub-000_ses-01>             - character
      --subject_dir=<subject derivatives>   - character
      --mica=<MICAPIPE path>                - character
      --dir_fs=<path to freesurfer>         - character
      --tracts=<number of tracts>           - character of numeric
      --help                                - print this text

      Example:
      micapipe.qc.R --BIDSid='sub-001' --subject_dir='derivatives/micapipe/sub-001' --mica='/Users/me/micapipe' --dir_fs='derivatives/freesurfer' --tracts=20M\n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}
# -----------------------------------------------------
#### Parse arguments ####
BIDSid <- argsL$BIDSid
subject_dir <- argsL$subject_dir
mica.dir <- argsL$mica
dir_fs <- argsL$dir_fs
tracts <- argsL$tracts
Filter <- "SIFT2"

# -----------------------------------------------------
#### Libraries ####
library('scales') # alpha function
library('ggplot2')
library('tidyr')
library('viridis')
library('plotly')
library('networkD3') # Sankey diagram
library('htmlwidgets') # save the widget

# -----------------------------------------------------
#### Functions ####
# --------------------------------------------------
module.col <- function() {
  # List of all the main outputs
  files <- data.frame(
    module=c(rep("proc_structural",7), rep("proc_freesurfer",2), rep("Morphology",12), rep("post_structural", 13), rep("proc_dwi", 13), rep("SC", 3), rep("proc_rsfmri", 20)),
    variables=c("t1w.nativepro", "t1w.firstout", "t1w.fast", "t1w.MNI0.8", "t1w.MNI2.0", "t1w.mask", "t1w.5tt", "fs.t1w", "fs.reconall",
                "morph.lh.thick", "morph.rh.thick", "morph.lh.thick.fsa5", "morph.rh.thick.fsa5", "morph.lh.thick.c69", "morph.rh.thick.c69",
                "morph.lh.curv", "morph.rh.curv", "morph.lh.curv.fsa5", "morph.rh.curv.fsa5", "morph.lh.curv.c69", "morph.rh.curv.c69",
                "t1w.fsnative", "t1w.cerebellum", "t1w.subcortical", "c69.lh.mid", "c69.rh.mid", "c69.lh.pial", "c69.rh.pial", "c69.lh.white",
                "c69.rh.white", "c69.lh.sphere", "c69.rh.sphere", "c69.lh.midth", "c69.rh.midth",
                "dwi_res", "dwi_corr", "T1nativepro_in_dwi", "dwi_mask", "dwi_dti", "dti_FA", "dti_ADC", "fod_wmN", "dwi_in_T1nativepro",
                "t1w_nativepro_in_dwi_NL", "dwi_5ttm", "dwi_gmwmi", "tdi_1M", "dwi_cere", "dwi_subc", "tck.tdi",
                "singleecho", "fmri_brain", "fmri_HP", "fmri_mean", "fix_output", "fmri_processed", "global_signal", "spikeRegressors","vol2surfTS.lh","out_surf_native.lh",
                "out_surf_fsa5.lh", "out_surf.lh", "vol2surfTS.rh", "out_surf_native.rh", "out_surf_fsa5.rh", "out_surf.rh", "rsfmri_subcortex", "timese_subcortex", "fmri_tSNR", "fmri.surf.c69.timeseries"
    ),
    files=c(paste0(proc_struct, "/", BIDSid, "_space-nativepro_t1w.nii.gz"),
            paste0(proc_struct, "/first/", BIDSid, "_space-nativepro_t1w_all_fast_firstseg.nii.gz"),
            paste0(proc_struct, "/", BIDSid, "_space-nativepro_t1w_brain_pve_0.nii.gz"),
            paste0(dir_warp, "/", BIDSid,  "_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz"),
            paste0(dir_warp, "/", BIDSid,  "_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_1Warp.nii.gz"),
            paste0(proc_struct, "/", BIDSid, "_space-nativepro_t1w_brain_mask.json"),
            paste0(proc_struct, "/", BIDSid, "_space-nativepro_t1w_5TT.nii.gz"),
            paste0(dir_fs, "/",BIDSid,"/mri/T1.mgz"),
            paste0(dir_fs, "/",BIDSid,"/scripts/recon-all.log"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-lh_thickness.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-rh_thickness.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-lh_thickness_10mm.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-rh_thickness_10mm.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-conte69-32k_desc-lh_thickness.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-conte69-32k_desc-rh_thickness.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-lh_curvature.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsnative_desc-rh_curvature.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-lh_curvature_10mm.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-fsaverage5_desc-rh_curvature_10mm.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-conte69-32k_desc-lh_curvature.mgh"),
            paste0(dir_morpho,"/",BIDSid,"_space-conte69-32k_desc-rh_curvature.mgh"),
            paste0(proc_struct, "/", BIDSid, "_space-fsnative_t1w.nii.gz"),
            paste0(dir_volum, "/", BIDSid, "_space-nativepro_t1w_atlas-cerebellum.nii.gz"),
            paste0(dir_volum, "/", BIDSid, "_space-nativepro_t1w_atlas-subcortical.nii.gz"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-lh_midthickness.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-rh_midthickness.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-lh_pial.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-rh_pial.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-lh_white.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-rh_white.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_lh_sphereReg.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_rh_sphereReg.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-lh_midthickness.surf.gii"),
            paste0(dir_conte69, "/", BIDSid, "_space-conte69-32k_desc-rh_midthickness.surf.gii"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-MP-PCA_residuals-dwi.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-dwi_preproc.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-t1w_nativepro.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-brain_mask.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_model-DTI.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_model-DTI_map-FA.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_model-DTI_map-ADC.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_model-CSD_map-FOD_desc-wmNorm.mif"),
            paste0(proc_struct, "/", BIDSid, "_space-nativepro_desc-dwi.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-t1w_nativepro_NL.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-5tt.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-gmwmi-mask.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-iFOD1-1M_tdi.mif"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_atlas-cerebellum.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_atlas-subcortical.nii.gz"),
            paste0(proc_dwi, "/", BIDSid, "_space-dwi_desc-iFOD2-",tracts,"_tdi.mif"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_brain.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_HP.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_mean.nii.gz"),
            paste0(rsfmri_ICA,"/filtered_func_data_clean.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_clean.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_global.txt"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_spikeRegressors_REFRMS.1D"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsnative_lh.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsnative_lh_10mm.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsaverage5_lh_10mm.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-conte69-32k_lh_10mm.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsnative_rh.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsnative_rh_10mm.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsaverage5_rh_10mm.mgh"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-conte69-32k_rh_10mm.mgh"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_subcortical.nii.gz"),
            paste0(rsfmri_volum, "/", BIDSid, "_space-rsfmri_desc-singleecho_timeseries_subcortical.txt"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_desc-tSNR.txt"),
            paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-conte69-32k_desc-timeseries_clean.txt")
    )
  )
  # Convert from factor to characters
  files[] <- lapply(files, as.character)

  # post_structural
  for (seg in parc) { files <- rbind(files, c("post_structural", paste0("t1w.annot",".",seg), paste0(dir_volum, "/", BIDSid, "_space-nativepro_t1w_atlas-",seg,".nii.gz")))}
  # GD
  for (seg in parc) { files <- rbind(files, c("GD", paste0("GD",".",seg), paste0(dir_geo,"/",BIDSid,"_space-fsnative_atlas-",seg,"_GD.txt")))}
  # MPC
  for (seg in parc) { files <- rbind(files, c("MPC", paste0("MPC",".",seg), paste0(dir_surf,"/micro_profiles/",BIDSid,"_space-fsnative_atlas-",seg,"_desc-MPC.txt")))}
  # SC
  for (seg in parc) { for (Size in c("cor", "sub", "full")) {
    files <- rbind(files, c("SC", paste0("SC",".",seg,".",Size), paste0(dwi_cnntm, "/", BIDSid, "_space-dwi_atlas-",seg,"_desc-iFOD2-",tracts,"-",Filter,"_",Size,"-connectome.txt")) )
  }}
  # proc_rsfmri
  for (seg in parc) { for (space in c("fsnative", "conte69-32k")) {
    files <- rbind(files, c("proc_rsfmri", paste0("FC",".",seg,".",space), paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-",space,"_atlas-",seg,"_desc-FC.txt")) )
  }}
  for (seg in parc) { files <- rbind(files, c("proc_rsfmri", paste0("FC",".",seg,".timeseries"), paste0(rsfmri_surf, "/", BIDSid, "_rsfmri_space-fsnative_atlas-",seg,"_desc-timeseries.txt")) )}
  # test existence of files
  files$done <- as.numeric(file.exists(files$files))

  # a2009 conte69 not supported
  files <- files[files$variables!="FC.aparc-a2009s.conte69-32k",]

  # files %>% filter(done == FALSE)
  Completed <- aggregate(done~module, data = files, mean )
  Completed <- Completed[match(nodes$name, Completed$module),]
  Completed.col <- c(Completed$done[3:9], Completed$done)
  Completed.col <- ifelse(Completed.col>0 & Completed.col<0.95,"#EE7942",Completed.col)
  Completed.col <- ifelse(Completed.col>=0.95,"#00688B",Completed.col)
  Completed.col <- ifelse(Completed.col==0,"#EEE9E9",Completed.col)
  Completed.col[1:7] <- alpha(Completed.col[1:7],0.35)
  return(Completed.col)
}

# -----------------------------------------------------
#### Varibles ####
# Paths
proc_struct <- paste0(subject_dir, "/anat")          # Structural processing directory
  dir_volum <- paste0(proc_struct, "/volumetric" )      # Cortical segmentations
  dir_surf <- paste0(proc_struct, "/surfaces")         # Structural surfaces
    dir_conte69 <- paste0(dir_surf, "/conte69")
    dir_mpc <- paste0(dir_surf, "/micro_profiles")       # MPC outputs
    dir_morpho <- paste0(dir_surf, "/morphology")        # morphometry
    dir_geo <- paste0(dir_surf, "/geo_dist")          # Geodesic distance outputs
proc_dwi <- paste0(subject_dir, "/dwi")              # dwi processing directory
  dir_eddy <- paste0(proc_dwi, "/eddy/dwi_post_eddy.eddy_") # Eddy directory
  dwi_cnntm <- paste0(proc_dwi, "/connectomes")
proc_rsfmri <- paste0(subject_dir, "/func")
  rsfmri_volum <- paste0(proc_rsfmri, "/volumetric" )
  rsfmri_surf <- paste0(proc_rsfmri, "/surfaces") # rsfMRI surfaces
  rsfmri_ICA <- paste0(proc_rsfmri, "/ICA_MELODIC")
dir_warp <- paste0(subject_dir, "/xfm")
  dir_fs.label <- paste0(dir_fs,"/",BIDSid,"/label")
dir_QC <- paste0(subject_dir, "/QC")                 # QC directory
  dir_QC_png <- paste0(subject_dir, "/QC/png")

# List of parcellations volumes (remove cerebellum,subcortical and keep the parcellation string name)
parc <- grep("nii.gz", list.files(dir_volum), value = TRUE)
parc <- grep("cerebellum", parc, invert=TRUE, value = TRUE)
parc <- gsub(".nii.gz", "", grep("subcortical", parc, invert=TRUE, value = TRUE))
parc <- unlist(lapply(parc, function(x) strsplit(x, "atlas-")[[1]][2]))

# -----------------------------------------------------
# Loads the MICs colormaps
load(file=paste0(mica.dir,"/functions/cmap_MICs.Rdata"))

# -----------------------------------------------------
#### Micapipe Workflow - sankey diagram ####
# Nodes
nodes <- data.frame(name=c("proc_structural","proc_freesurfer","proc_dwi","post_structural","proc_rsfmri","SC","MPC","GD","Morphology"),
                    color=c("#104E8B",  "#00688B",  "#668B8B",  "#7A8B8B",  "#B3362C",  "#6F4FA3",  "#3B8549",  "#3B70A2",  "#cd950c"), # colored nodes
                    group=c("proc_structural", "proc_freesurfer","proc_dwi","post_structural","proc_rsfmri","SC","MPC","GD","Morphology"))
nodes$name <- as.character(nodes$name)
nodes$color <- as.character(nodes$color)

# Links
links <- data.frame(source=c(3, 3, 3, 3, 0, 2, 1, 1, 1, 1, 1, 0, 0, 0),
                    target=c(4, 5, 6, 7, 8, 5, 3, 4, 6, 7, 8, 2, 3, 4),
                    value=c(66.66, 66.66, 66.66, 66.66, 66.66, 66.66, 50.00, 33.33, 33.33, 33.33, 33.33, 50.00, 50.00, 33.33),
                    group=as.factor(c(4, 5, 6, 7, 8, 5, 3, 4, 6, 7, 8, 2, 3, 4)) ) # colored Right to Left
                    # group=as.factor(c(3, 3, 3, 3, 0, 2, 1, 1, 1, 1, 1, 0, 0, 0)) ) # colored Left to Right && modules <- c(nodes$name[1:4], nodes$name)

# Dinamic Color change Nodes$color or Grays
modules <- c(nodes$name[3:9], nodes$name)
ColMap <- module.col()

# Colouring in java
Colors <- paste0('d3.scaleOrdinal() .domain([\"',paste(c(levels(links$group), nodes$name),collapse='\",\"'),'\"])','.range([\"',paste(ColMap,collapse='\",\"'),'\"])')
# Plot the diagram
p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", NodeGroup = "group", LinkGroup = "group",
              units = "", fontSize = 14, nodeWidth = 30, colourScale=Colors, fontFamily = "Courier")
# save the diagram
saveWidget(p, file=paste0(dir_QC,"/",BIDSid ,"_desc-qc_micapipe_workflow.html"))

# -----------------------------------------------------
#### Plot the connectomes ####
# Load the connectomes fuction
load.conn <- function(PATH.txt, sym=TRUE) {
  M <- as.matrix(read.csv(PATH.txt, sep = " ",header = FALSE,))
  if (sym==TRUE) {
    # Mirror matrix completes inferior triangle
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
  }
  return(M)
}

# Create a PNG per connectome
notFound <- matrix(data=0, nrow = 256, ncol = 256)

# Creating QC png connectomes for all modalities
for (seg in parc) { subj_id <- paste0(BIDSid,"_atlas-",seg)
    print(paste("[INFO]....  Creating PNG connectomes from parcellation:", seg))

    # Functional connectome
    conn.rsf <- paste0(rsfmri_surf,"/",BIDSid,"_rsfmri_space-fsnative_atlas-",seg,"_desc-FC.txt")
      if (file.exists(conn.rsf)==TRUE) { conn.rsf <- load.conn(conn.rsf); File <- ""; ColMap <- cmap.FC(256)
          png(paste0(dir_QC_png,"/", subj_id, "_desc-qc_FC.png"))
          image(conn.rsf, axes=FALSE, main=paste0(subj_id,"_FC"), col=ColMap ); dev.off()
          } else { conn.rsf <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30"; print(paste(File, seg, "FC")) }

    # Functional connectome conte69
    conn69.rsf <- paste0(rsfmri_surf,"/",BIDSid,"_rsfmri_space-conte69-32k_atlas-",seg,"_desc-FC.txt")
      if (file.exists(conn69.rsf)==TRUE) { conn69.rsf <- load.conn(conn69.rsf); File <- ""; ColMap <- cmap.FC(256)
          png(paste0(dir_QC_png,"/", subj_id,"_desc-qc_FC-conte69.png"))
          image(conn69.rsf, axes=FALSE, main=paste0(File, subj_id,"_FC-conte69.png"), col=ColMap ); dev.off()
      } else { conn69.rsf <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30"; print(paste(File, seg, "FC conte 69")) }

    # Geodesic distance
    conn.geo <- paste0(dir_geo,"/",BIDSid,"_space-fsnative_atlas-",seg,"_GD.txt")
      if (file.exists(conn.geo)==TRUE) { conn.geo <- load.conn(conn.geo, sym=FALSE); File <- ""; ColMap <- cmap.GD(256)
          png(paste0(dir_QC_png,"/", subj_id, "_desc-qc_GD.png"))
          image(conn.geo, axes=FALSE, main=paste0(subj_id,"_GD"), col=ColMap ); dev.off()
       } else { conn.geo <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30"; print(paste(File, seg, "GD")) }

    # Micro structural profiles
    conn.mpc <- paste0(dir_surf,"/micro_profiles/",BIDSid,"_space-fsnative_atlas-",seg,"_desc-MPC.txt")
      if (file.exists(conn.mpc)==TRUE) { conn.mpc <- load.conn(conn.mpc); File <- ""; ColMap <- cmap.MPC(256)
          png(paste0(dir_QC_png,"/", subj_id, "_desc-qc_MPC.png"))
          image(conn.mpc, axes=FALSE, main=paste0(BIDSid," ",seg,"MPC"), col=ColMap ); dev.off()
      } else { conn.mpc <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30"; print(paste(File, seg, "MPC")) }

    # Intensity profiles
    conn.int <- paste0(dir_surf,"/micro_profiles/",BIDSid,"_space-fsnative_atlas-",seg,"_desc-intensity_profiles.txt")
      if (file.exists(conn.int)==TRUE) { conn.int <- load.conn(conn.int, sym=FALSE); File <- ""; ColMap <- cmap.MPC(256)
          png(paste0(dir_QC_png,"/", subj_id, "_desc-qc_MPC-intensity.png"), height = 250)
          image(t(conn.int), axes=FALSE, main=paste0(subj_id,"_Intensity"), col=ColMap ); dev.off()
      } else { conn.int <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30"; print(paste(File, seg, "intensity")) }
}
