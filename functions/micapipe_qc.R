# Libraries
library('scales')
library('ggplot2')
library('tidyr')
library('viridis')
library('plotly')

# -----------------------------------------------------
# Variables
id <- "HC10"
ses <- "ses-pre"
out <- "/Users/rcruces/tmp/derivatives"


# Resting state path and name
subject_dir <- paste0(out, "/sub-", id, "/", ses)
  proc_struct <- paste0(subject_dir, "/proc_struct") # Structural processing directory
    dir_volm <- paste0(proc_struct, "/volumetric" )  # Cortical segmentations
    dir_surf <- paste0(proc_struct, "/surfaces")     # Structural surfaces
  dir_QC <- paste0(subject_dir, "/QC")              # QC directory
  dir_QC_png <- paste0(subject_dir, "/QC/png")
    rsfmri_surf <- paste0(subject_dir, "/proc_rsfmri/surfaces") # rsfMRI surfaces
  proc_dwi <- paste0(subject_dir, "/proc_dwi")       # dwi processing directory
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
# mica.dir <- paste0(argsL$mica,"/functions/")
mica.dir <- "/Users/rcruces/git_here/micapipe/functions"
load(file=paste0(mica.dir,"/cmap_MICs.Rdata"))

# Load the connectomes
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
# Sankey Diagram
# https://www.r-graph-gallery.com/sankey-diagram.html


# -----------------------------------------------------
# FSL EDDY parameters

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
