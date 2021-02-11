# -----------------------------------------------------
#
# Generates PNG images for micapipe Quality Check
# Derived from a NIFTI
# mipapipe v0.0.1
# Extra functions 
# 
# Created by RRC on Feb 2021 (the second year of the pademic)
#
# -----------------------------------------------------
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      micapipe PNG creator

      Arguments:
      --img=<path to image file>           - character
      --out=<path to output image PNG>     - character
      --overlay=<NIFTI>                    - character
      --overlayThr <3rd Qt by default>     - numeric
      --overlayCol <purple by default>     - factor: 'viridis', 'cividis', 'plasma', 'magma', 'inferno', 'set100'
      --color <grayscale by default>       - factor: 'viridis', 'cividis', 'plasma', 'magma', 'inferno', 'set100'
      --help                               - print this text

      Example of a registration:
      micapipe_qc-png.R --in='sub-HC62_ses-01_space-t1nativepro_t1w.nii.gz' --out='sub-HC62_ses-01_space-t1nativepro_FA-on-t1w.png' --overlay='sub-HC62_ses-01_space-t1nativepro_FA.nii.gz' 
      \n\n")
        
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
arg.miss <- function(arg) {
  if(is.null(arg)) {
    cat("---------------------------------------------------------
    ERROR ... A mandatory argument is missing:
      ---------------------------------------------------------
      ")
    print(paste0("img       : ",argsL$img))
    print(paste0("out      : ",argsL$out))
    exit()
  }
}

check.file <- function(arg) {
  if(!file.exists(arg)) {
    cat("---------------------------------------------------------
        ")
    print(paste0("ERROR ... File does not exist: ",arg))
    cat("---------------------------------------------------------
        ")
    exit()
  }
}

# -----------------------------------------------------
# Libraries
require("oro.nifti")
require("neurobase")
require("scales")
require("viridis")

# Colormaps
colmaps <- list('viridis'=viridis(100), 'cividis'=cividis(100), 'plasma'=plasma(100), 'magma'=magma(100), 'inferno'=inferno(100), 
                'set100'=colorRampPalette(c("black", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"))(100))
'%!in%' <- function(x,y)!('%in%'(x,y)) #not in function

# Main image color
if (is.null(argsL$color)) { Col <- gray(0:64/64) 
} else {
  if(as.character(argsL$color) %!in% names(colmaps) ) {
    cat("--overlayCol is not supported")
    Col <- gray(0:64/64)
  } else {
    Col <- unlist(colmaps[argsL$color])
  }
}

# Overlay image color
if (is.null(argsL$overlayCol)) { Col.y <- alpha("darkorchid1", 0.45)
} else {
  if(as.character(argsL$overlayCol) %!in% names(colmaps) ) {
    cat("--overlayCol is not supported")
    Col.y <- alpha("darkorchid1", 0.45)
  } else {
    Col.y <- alpha(unlist(colmaps[argsL$overlayCol]), 0.45)
  }
}

## Case: if argument is missing
arg.miss(argsL$img)
arg.miss(argsL$out)

## Case if file does not exist
check.file(argsL$img)

# Print arguments
print(paste0("--img         : ",argsL$img))
print(paste0("--out         : ",argsL$out))
if (!is.null(argsL$overlay)) { print(paste0("--overlay     : ",argsL$overlay)) }
if (!is.null(argsL$overlay.thr)) { print(paste0("--overlay.thr : ", argsL$overlay.thr)) }
if (!is.null(argsL$color)) { print(paste0("--color : ", argsL$color)) }
if (!is.null(argsL$overlayCol)) { print(paste0("--overlayCol : ", argsL$overlayCol)) }

# out name
nom <- gsub(pattern = ".png", replacement = "", argsL$out)

# -----------------------------------------------------
# First case only prints image
if ( is.null(argsL$overlay) & file.exists(argsL$img) ) {
  print(paste0("INFO.... Creating PNG with orthographic view of a NIFTI"))
  # Load Image
  nii <- readNIfTI(argsL$img, reorient = FALSE)

  # Plot NIFTI and Save png for QC purposes
  print(paste0("INFO.... Saving nifti image as PNG"))
  png(paste0(nom, ".png"))
  orthographic(nii, col.crosshairs = alpha("aliceblue", 0.2), axes=TRUE, col=Col)
  dev.off()
}

# -----------------------------------------------------
# Second case prints NIFTI and overlay
if ( !is.null(argsL$overlay) & file.exists(argsL$img) ) {
    check.file(argsL$overlay)
    # Load Images
    nii <- readNIfTI(argsL$img, reorient = FALSE)
    over <- readNIfTI(argsL$overlay, reorient = FALSE)
    # Set the threshold for the overlay image
    if (is.null(argsL$overlay.thr)){ thr <- summary(over@.Data[over@.Data!=0])[5] }
    else { thr <- argsL$overlay.thr}
    print(paste0("INFO.... Creating PNG with orthographic view of a NIFTI with overlay, thr=",thr))    
    # Plot NIFTI and overlay  and Save png for QC purposes
    print(paste0("INFO.... Saving nifti and overlay image as PNG"))
    png(paste0(nom, ".png"))
    ortho2(nii, y = over > thr, col.y = Col.y, col.crosshairs = alpha("aliceblue", 0.2), col=Col)
    dev.off()
}
