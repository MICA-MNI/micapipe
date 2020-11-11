
# -----------------------------------------------------
# Variables
id <- "HC10"
ses <- "ses-pre"
out <- "/Users/rcruces/tmp/derivatives"
mica <- 

# Resting state path and name
subject_dir <- paste0(out, "/sub-", id, "/", ses)  
  proc_struct <- paste0(subject_dir, "/proc_struct") # Structural processing directory
    dir_volm <- paste0(proc_struct, "/volumetric" )  # Cortical segmentations
    dir_surf <- paste0(proc_struct, "/surfaces")     # Structural surfaces
  dir_QC <- paste0(subject_dir, "/QC/")
    rsfmri_surf <- paste0(subject_dir, "/proc_rsfmri/surfaces")

# Directories of connectomes
dir_mpc <- paste0(rsfmri_surf, "/micro_profiles")
dir_geo <- paste0(rsfmri_surf, "/geo_dist")

# List of parcellations (remove cerebellum,subcortical and keep the parcellation string name)
parc <- grep("cerebellum", list.files(dir_volm), invert=TRUE, value = TRUE)
parc <- gsub(".nii.gz", "", grep("subcortical", parc, invert=TRUE, value = TRUE))
parc <- unlist(lapply(parc, function(x) strsplit(x, "nativepro_")[[1]][2]))

# -----------------------------------------------------
# Loads the MICs colormaps
# mica.dir <- paste0(argsL$mica,"/functions/")
mica.dir <- "/Users/rcruces/git_here/micapipe/functions/"
load(file=paste0(mica.dir,"cmap_MICs.Rdata"))

# Load the connectomes 
load.conn <- function(PATH.txt, sym=TRUE) {
  M <- as.matrix(read.table(PATH.txt, sep = " ",header = FALSE))
  if (sym==TRUE) {
    # Mirror matrix completes inferior triangle
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
  }
  return(M)
}

# Create a PNG per connectome
for (seg in parc) {
  # Functional connectome
  conn.rsf <- load.conn(paste0(rsfmri_surf,"/",id,"_rsfMRI-connectome_",seg,"_clean.txt"))
  png(paste0(dir_QC, id, "_", seg, "_FC", ".png")) 
  image(conn.rsf, axes=FALSE, main=paste0(id,"_",seg,"_FC"), col=cmap.FC(256)); dev.off()
  # Geometric distance
  conn.geo <- load.conn(paste0(dir_surf,"/geo_dist/",id,"_",seg,"_GD.txt"), sym=FALSE)
  png(paste0(dir_QC, id, "_", seg, "_GD", ".png"))
  image(conn.geo, axes=FALSE, main=paste0(id,"_",seg,"_GD"), col=cmap.GD(256)); dev.off()
  # Micro structural profiles
  conn.mpc <- load.conn(paste0(dir_surf,"/micro_profiles/mpc_",seg,"_mics.txt"))
  png(paste0(dir_QC, id, "_", seg, "_MPC", ".png"))
  image(conn.mpc, axes=FALSE, main=paste0(id,"_",seg,"_MPC"), col=cmap.MPC(256)); dev.off()
  # Intensity profiles
  conn.int <- load.conn(paste0(dir_surf,"/micro_profiles/intensity_profiles_",seg,"_mics.txt"))
  png(paste0(dir_QC, id, "_", seg, "_Intensity", ".png"), height = 250)
  image(t(conn.int), axes=FALSE, main=paste0(id,"_",seg,"_Intensity"), col=cmap.MPC(256)); dev.off()
}


# -----------------------------------------------------
# Sankey Diagram
# https://www.r-graph-gallery.com/sankey-diagram.html

  
# -----------------------------------------------------
# Eddy parameters
"Total movement RMS"
# RMS relative to the first volume 
# second column the RMS relative the previous volume
dwi_post_eddy.eddy_movement_rms

"Outliers"
dwi_post_eddy.eddy_outlier_report

""
# One row per scan, one column per slice. b0s set to zero
dwi_post_eddy.eddy_outlier_n_sqr_stdev_map

dwi_post_eddy.eddy_outlier_map
# All numbers are either 0, meaning that scan-slice is not an outliers, or 1 meaning that is is.

dwi_post_eddy.eddy_outlier_n_stdev_map
# The numbers denote how many standard deviations off the mean difference between observation and prediction is.

dwi_post_eddy.eddy_outlier_n_sqr_stdev_map
# The numbers denote how many standard deviations off the square root of the mean squared difference between observation and prediction is.

eddy_parameters
This is a text file with one row for each volume in --imain and one column for each parameter. The first six columns correspond to subject movement starting with three translations followed by three rotations. The remaining columns pertain to the EC-induced fields and the number and interpretation of them will depend of which EC model was specified.