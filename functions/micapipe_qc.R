# R script
#
# Generates PNG images for micapipe Quality Check
# mipapipe v0.0.1
# Extra functions 
# 
# Created by RRC on Jan-Feb 2021 (the second year of the pademic)
#

# -----------------------------------------------------
#### Libraries ####
library('scales') # alpha function
library('ggplot2')
library('tidyr')
library('viridis')
library('plotly')
library('networkD3') # Sankey diagram
library('fsbrain')
library('rgl')

# -----------------------------------------------------
#### Functions ####

# Overlay Fixed image to Moving

# MRI example


# -----------------------------------------------------
#### Variables ####
id <- "HC10"
subject_dir <- "/Users/rcruces/tmp/derivatives-single/sub-HC10"

# Paths
proc_struct <- paste0(subject_dir, "/proc_struct") # Structural processing directory
dir_volm <- paste0(proc_struct, "/volumetric" )  # Cortical segmentations
dir_surf <- paste0(proc_struct, "/surfaces")     # Structural surfaces
dir_QC <- paste0(subject_dir, "/QC")              # QC directory
dir_QC_png <- paste0(subject_dir, "/QC/png")
rsfmri_surf <- paste0(subject_dir, "/proc_rsfmri/surfaces") # rsfMRI surfaces
proc_dwi <- paste0(subject_dir, "/proc_dwi")       # dwi processing directory
dir_mpc <- paste0(subject_dir, "/QC")              # MPC
dir_morpho <- paste0(subject_dir, "/QC")              # morphometry
dir_eddy <- paste0(proc_dwi, "/eddy/dwi_post_eddy.eddy_") # Eddy directory

# Directories of connectomes
dir_mpc <- paste0(rsfmri_surf, "/micro_profiles")    # MPC outputs
dir_geo <- paste0(rsfmri_surf, "/geo_dist")          # Geodesic distance outputs

# List of parcellations (remove cerebellum,subcortical and keep the parcellation string name)
parc <- grep("nii.gz", list.files(dir_volm), value = TRUE)
parc <- grep("cerebellum", parc, invert=TRUE, value = TRUE)
parc <- gsub(".nii.gz", "", grep("subcortical", parc, invert=TRUE, value = TRUE))
parc <- unlist(lapply(parc, function(x) strsplit(x, "nativepro_")[[1]][2]))

# -----------------------------------------------------
# Loads the MICs colormaps
#mica.dir <- paste0(argsL$mica,"/functions/")
mica.dir <- "/Users/rcruces/git_here/micapipe"
load(file=paste0(mica.dir,"/functions/cmap_MICs.Rdata"))

# -----------------------------------------------------
# Micapipe Workflow
# MICAPIPE sankey diagram
links <- read.csv("~/Desktop/mica_sankey-links.csv")
Nodes <- read.csv("~/Desktop/mica_sankey-nodes.csv", as.is = TRUE)

# Node group 
Nodes$group <- as.factor(Nodes$name)
links$group <- as.factor(links$source)

# Dinamic Color change Nodes$color or Grays
modules <- c(Nodes$name[1:4], Nodes$name)
process <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
ColMap.full <- c(alpha(Nodes$color[1:4],0.3), Nodes$color)
ColMap <- ifelse(process==TRUE, ColMap.full, "#EEE9E9")
# Colouring in java
Colors <- paste0('d3.scaleOrdinal() .domain([\"',paste(c(levels(links$group), Nodes$name),collapse='\",\"'),'\"])','.range([\"',paste(ColMap,collapse='\",\"'),'\"])')
# Plot the diagram
p <- sankeyNetwork(Links = links, Nodes = Nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", NodeGroup = "group", LinkGroup = "group",
              units = "", fontSize = 14, nodeWidth = 30, colourScale=Colors, fontFamily = "Courier")

# save the widget
library(htmlwidgets)
saveWidget(p, file="/Users/rcruces/tmp/derivatives-single/sub-HC10/QC/micapipe_workflow.html")

# -----------------------------------------------------
# Surface with values
# Visualize: https://cran.r-project.org/web/packages/freesurferformats/vignettes/freesurferformats.html
# Reads: http://127.0.0.1:20380/library/fsbrain/doc/fsbrain.html
# rgla <- list('trans_fun'=clip.data); # fit color to data
library('fsbrain')
library('viridis')
library('rgl')

# Variables
subject_dir = '/Users/rcruces/tmp/derivatives-single/micapipe_v0.0.1/sub-HC10/proc_struct/surfaces'
subject = 'HC10'
dir_QC_png ='/Users/rcruces/tmp/micapipe_qc/'

# proc_freesurfer - Annotation files from FREESURFER: aparc.a2009s
for (annot in c('aparc', 'aparc.a2009s')) {
  nom <- paste0(subject,"_space-fsnative_",annot,".png"); nom.png <- paste0(dir_QC_png, nom)
  if (!file.exists(nom.png)) {
    print(paste0("INFO.... Creating PNG of ", annot, " on native surface"))
    cm <- vis.subject.annot(subject_dir, subject, annot, 'both', surface="pial", rglactions= list('no_vis'=T))
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, img_only = TRUE);
    while (rgl.cur() > 0) { rgl.close() }
  } else {
    print(paste0("INFO.... File exists: ", nom))
  }
}

####  fs - Sulcal depth #### 
nom <-  paste0(subject,"_space-fsnative_sulc.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating PNG of the SULCAL DEPTH on native surface"))
  rgla <- list('trans_fun'=limit_fun(-5, 5), 'no_vis'=T);
  cm <- vis.subject.morph.native(subject_dir, subject, 'sulc', cortex_only = FALSE, 
                                 draw_colorbar = TRUE, surface="white", rglactions = rgla, makecmap_options = list('colFn'=cividis))
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Sulcal depth [mm]', quality = 1L);
  while (rgl.cur() > 0) { rgl.close() } 
} else { print(paste0("INFO.... File exists: ", nom)) }


####  fs - Cortical thickness #### 
nom <-  paste0(subject,"_space-fsnative_th.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating PNG of the THICKNESS on native surface"))
  rgla <- list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T);
  cm <- vis.subject.morph.native(subject_dir, subject, 'thickness', cortex_only = FALSE, 
                                 draw_colorbar = TRUE, surface="pial", rglactions = rgla, makecmap_options = list('colFn'=inferno))
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Thickness [mm]');
  while (rgl.cur() > 0) { rgl.close() }
} else { print(paste0("INFO.... File exists: ", nom)) }


#### fs - Curvature #### 
col.curv <- colorRampPalette(c('darkolivegreen3', 'gray50', 'firebrick3'));
nom <-  paste0(subject,"_space-fsnative_curv.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating PNG of the CURVATURE on native surface"))
  rgla <- list('trans_fun'=limit_fun(-0.2, 0.2), 'no_vis'=T);
  cm <- vis.subject.morph.native(subject_dir, subject, 'curv', cortex_only = FALSE, 
                                 draw_colorbar = TRUE, surface="white", rglactions = rgla, makecmap_options = list('colFn'=col.curv))
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Curvature');
  while (rgl.cur() > 0) { rgl.close() }
} else { print(paste0("INFO.... File exists: ", nom)) }


# -----------------------------------------------------
#### Morphology - Cortical thickness #### 
morpho <- "/Users/rcruces/tmp/derivatives-single/micapipe_v0.0.1/sub-HC10/proc_struct/surfaces/morphology"

# Native surface inflated
nom <-  paste0(subject,"_morphology_space-fsnative_th.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating ",nom, " of the THICKNESS on native surface"))
  rgla <- list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T);
  cm <- vis.data.on.subject(subject_dir, subject, morph_data_lh=paste0(morpho,"/lh_thickness.mgh"), morph_data_rh=paste0(morpho,"/rh_thickness.mgh"), 
                            draw_colorbar = TRUE, surface="inflated", rglactions = rgla,  makecmap_options = list('colFn'=inferno))
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Thickness [mm]');
  while (rgl.cur() > 0) { rgl.close() }
} else { print(paste0("INFO.... File exists: ", nom)) }

# fsaverage5 surface
nom <-  paste0(subject,"_morphology_space-fsaverage5_th.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating ",nom, " of the THICKNESS on fsaverage5"))
  rgla <- list('trans_fun'=limit_fun(1.5, 4), 'no_vis'=T);
  cm <- vis.data.on.subject(subject_dir, 'fsaverage5', morph_data_lh=paste0(morpho,"/lh.thickness_fsa5.mgh"), morph_data_rh=paste0(morpho,"/rh.thickness_fsa5.mgh"), 
                            draw_colorbar = TRUE, surface="pial", rglactions = rgla,  makecmap_options = list('colFn'=inferno))
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Thickness [mm]');
  while (rgl.cur() > 0) { rgl.close() }
} else { print(paste0("INFO.... File exists: ", nom)) }


# -----------------------------------------------------
####  Morphology - Curvature ####
nom <-  paste0(subject,"_morphology_space-fsnative_curv.png"); nom.png <- paste0(dir_QC_png, nom)
if (!file.exists(nom.png)) {
  print(paste0("INFO.... Creating PNG of the CURVATURE on native surface"))
  rgla <- list('trans_fun'=limit_fun(-0.3, 0.3), 'no_vis'=T);
  cm <- vis.data.on.subject(subject_dir, subject, morph_data_lh=paste0(morpho,"/lh_curv.mgh"), morph_data_rh=paste0(morpho,"/rh_curv.mgh"), 
                            draw_colorbar = 'horizontal',  rglactions=rgla, makecmap_options = list('colFn'=col.curv), surface="inflated")
  img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = 'Curvature');
  while (rgl.cur() > 0) { rgl.close() }
} else { print(paste0("INFO.... File exists: ", nom)) }


# -----------------------------------------------------
#### Microstructural profile covariance #### 
mcp <- "/Users/rcruces/tmp/derivatives-single/micapipe_v0.0.1/sub-HC10/proc_struct/surfaces/micro_profiles"
for (i in 1:14) {
  i0 <- sprintf("%02d", i)
  nom <-  paste0(subject,"_space-fsnative_mpc-",i0,".png"); nom.png <- paste0(dir_QC_png, nom)
  if (!file.exists(nom.png)) {
    print(paste0("INFO.... Creating PNG of MPC-", i0))
    rgla <- list('trans_fun'=limit_fun(1200,2000), 'no_vis'=T);
    cm <- vis.data.on.subject(subject_dir, subject, morph_data_lh=paste0(mcp,"/lh.",i,".mgh"), morph_data_rh=paste0(mcp,"/rh.",i,".mgh"), surface = "inflated", draw_colorbar = 'vertical', rglactions = rgla)
    img <- vis.export.from.coloredmeshes(cm, output_img = nom.png, colorbar_legend = paste0('MPC-',i0));
    while (rgl.cur() > 0) { rgl.close() }
  } else { print(paste0("INFO.... File exists: ", nom)) }
}

# -----------------------------------------------------
# Parcellations
# vis.subject.annot(subject_dir, subject, 'economo_mics', 'both', surface="pial")
# vis.subject.annot(subject_dir, subject, 'schaefer-100_mics', 'both', surface="pial")

# -----------------------------------------------------
#### Plot the curvature & save figure #### 
# rgla <- list('trans_fun'=clip.data, 'snapshot_png'="~/Desktop/fsbrain.png");
# vis.subject.morph.native(subject_dir, subject, 'curv', cortex_only = TRUE, draw_colorbar = 'horizontal',  
#                          rglactions=rgla, makecmap_options = list('colFn'=magma), surface="pial")
# rgl.close()

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
for (seg in parc) { subj_id <- paste0(id,"_",seg)
    print(paste("[INFO]....  Creating PNG connectomes from parcellation:", subj_id))

    # Functional connectome
    conn.rsf <- paste0(rsfmri_surf,"/",id,"_rsfMRI-connectome_",seg,"_clean.txt")
      if (file.exists(conn.rsf)==TRUE) { conn.rsf <- load.conn(conn.rsf); File <- ""; ColMap <- cmap.FC(256)
          } else { conn.rsf <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30" }
      png(paste0(dir_QC_png,"/", subj_id, "_FC.png"))
      image(conn.rsf, axes=FALSE, main=paste0(subj_id,"_FC"), col=ColMap ); dev.off()

    # Functional connectome conte69
    conn69.rsf <- paste0(rsfmri_surf,"/",id,"_rsfMRI-connectome_",seg,"_conte69_clean.txt")
      if (file.exists(conn69.rsf)==TRUE) { conn69.rsf <- load.conn(conn69.rsf); File <- ""; ColMap <- cmap.FC(256)
          } else { conn69.rsf <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30" }
      png(paste0(dir_QC_png,"/", subj_id,"_FC-conte69.png"))
      image(conn69.rsf, axes=FALSE, main=paste0(File, subj_id,"_FC-conte69.png"), col=ColMap ); dev.off()

    # Geodesic distance
    conn.geo <- paste0(dir_surf,"/geo_dist/",subj_id,"_GD.txt")
      if (file.exists(conn.geo)==TRUE) { conn.geo <- load.conn(conn.geo, sym=FALSE); File <- ""; ColMap <- cmap.GD(256)
          } else { conn.geo <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30" }
      png(paste0(dir_QC_png,"/", subj_id, "_GD.png"))
      image(conn.geo, axes=FALSE, main=paste0(subj_id,"_GD"), col=ColMap ); dev.off()

    # Micro structural profiles
    conn.mpc <- paste0(dir_surf,"/micro_profiles/mpc_",seg,"_mics.txt")
      if (file.exists(conn.mpc)==TRUE) { conn.mpc <- load.conn(conn.mpc); File <- ""; ColMap <- cmap.MPC(256)
          } else { conn.mpc <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30" }
      png(paste0(dir_QC_png,"/", subj_id, "_MPC.png"))
      image(conn.mpc, axes=FALSE, main=paste0(subj_id,"_MPC"), col=ColMap ); dev.off()

    # Intensity profiles
    conn.int <- paste0(dir_surf,"/micro_profiles/intensity_profiles_",seg,"_mics.txt")
      if (file.exists(conn.int)==TRUE) { conn.int <- load.conn(conn.int, sym=FALSE); File <- ""; ColMap <- cmap.MPC(256)
          } else { conn.int <- notFound; File <- "NOT FOUND - "; ColMap <- "gray30" }
      png(paste0(dir_QC_png,"/", subj_id, "_MPC-intensity.png"), height = 250)
      image(t(conn.int), axes=FALSE, main=paste0(subj_id,"_Intensity"), col=ColMap ); dev.off()
}


# -----------------------------------------------------
#### DWI - FSL EDDY parameters #### 

# Read the gradient table information
grad <- read.table(paste0(proc_dwi,"/",id,"_dwi_corr-grad.txt"), sep = "",header = FALSE)
colnames(grad) <- c("x", "y", "z", "shell")
# Round shells
grad$shell <- round(grad$shell/100)*100
grad$volume <- ave(grad$shell, grad$shell, FUN=seq_along)

# RMS relative to the first volume
# second column the RMS relative the previous volume
RMS <- read.table(paste0(dir_eddy, "movement_rms"), sep = "",header = FALSE)
colnames(RMS) <- c("Absolute", "Relative")
RMS$shell <- grad$shell
RMS$Volume <- grad$volume
RMS$seq <- 1:dim(RMS)[1]

RMS <- gather(RMS, "Displacement", "mm", 1:2)
# For each shell (non 0)
ggplot(RMS, aes(x=seq, y=mm, group=Displacement)) +
  geom_line(aes(color=Displacement))+
  ggtitle("Estimated mean displacement")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Volume") + ylab("Displacement [mm]")


# One row per scan, one column per slice. b0s set to zero
# The numbers denote how many standard deviations off the square root of the mean squared difference between observation and prediction is.
M <- as.matrix(read.table(paste0(dir_eddy, "outlier_n_sqr_stdev_map"), sep = "", header = FALSE, skip = 1))
rownames(M) <- 1:dim(M)[1]
colnames(M) <- 1:dim(M)[2]
plot(M[,1], type='line', main="outlier_n_sqr_stdev_map", bty='n', ylab = "")
gather(M)
apply(M[-1,], 2, function(x) lines(x, col=alpha("gray65", 0.5)))
image(M, x = 1:dim(M)[1], y = 1:dim(M)[2], ylab = "Slice", xlab = "Volume", main="No. std.devs from mean slice-difference", col = viridis(50), breaks = seq(-1,4,0.1))

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6264528/
# This is a text file with one row for each volume in --imain and one column for each parameter.
# The first six columns correspond to subject movement starting with three translations followed by three rotations.
# The remaining columns pertain to the EC-induced fields and the number and interpretation of them will depend of which EC model was specified.
eddy <- read.table(paste0(dir_eddy, "parameters"), sep = "", header = FALSE, skip = 1)
colnames(eddy)[1:9] <- c("x.tra", "y.tra", "z.tra", "x.rot", "y.rot", "z.rot", "x.ec", "y.ec", "z.ec")

# Eddy currents
plot(eddy$x.ec, type='line', col="red", bty='n', main="Eddy currents linear term (x-axis)", xlab="Volume", ylab="Hz/mm")
plot(eddy$y.ec, type='line', col="red", bty='n', main="Eddy currents linear term (y-axis)", xlab="Volume", ylab="Hz/mm")
plot(eddy$z.ec, type='line', col="red", bty='n', main="Eddy currents linear term (z-axis)", xlab="Volume", ylab="Hz/mm")

# Color by sitance from origin
D <- apply(eddy[,7:9], 1, function(x) dist(rbind(c(0,0,0), x)))
plot_ly(x = cumsum(eddy$x.ec), y = cumsum(eddy$y.ec), z = cumsum(eddy$z.ec), type = 'scatter3d', mode = 'lines',
               opacity = 1, line = list(width = 6, color = D, reverscale = FALSE),
              marker = list(size = 3.5, color = D))


D <- apply(eddy[,4:6], 1, function(x) dist(rbind(c(0,0,0), x)))
plot_ly(x = cumsum(eddy$x.rot), y = cumsum(eddy$y.rot), z = cumsum(eddy$z.rot), type = 'scatter3d', mode = 'lines',
        opacity = 1, line = list(width = 6, color = D, reverscale = FALSE,colorscale = 'Viridis'),
        marker = list(size = 3.5, color = D, colorscale = 'Viridis'))
