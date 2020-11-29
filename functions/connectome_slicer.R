# -----------------------------------------------------
# Slices a connectome matrix based on lookup tables
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
      micapipe connectome slicer

      Arguments:
      --conn=<path to connectivty matrix>   - character
      --lut1=<path to LUT>                  - character
      --lut2=<path to LUT>                  - character
      --mica=<MICAPIPE path>                - character
      --help                                - print this text

      Example:
      connectome_slicer.R --conn='HC10_10M_glasser-360_full-connectome.txt' --lut1='lut_subcortical-cerebellum_mics.csv' --lut2='lut_glasser-360_mics.csv' \n\n")

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
arg.miss <- function() {
  cat("---------------------------------------------------------
    ERROR ... Connectome slicer: A mandatory argument is missing
      ---------------------------------------------------------
      ")
  print(paste0("Connectome     : ",argsL$conn))
  print(paste0("LUT-1          : ",argsL$lut1))
  print(paste0("LUT-2          : ",argsL$lut2))
  print(paste0("MICAPIPE-path  : ",argsL$mica))
}

## Arg1 default
if(is.null(argsL$conn)) {
  arg.miss()
  exit()
}

## Arg2 default
if(is.null(argsL$lut1)) {
  arg.miss()
  exit()
}

## Arg3 default
if(is.null(argsL$lut2)) {
  arg.miss()
  exit()
}

## Arg3 default
if(is.null(argsL$mica)) {
  arg.miss()
  exit()
}

# -----------------------------------------------------

# Subcortical LUT
lut1 <- read.csv(argsL$lut1)

# Cortical LUT
lut2 <- read.csv(argsL$lut2)
indx <- sort(c(lut1$mics, lut2$mics))

# Connectivity matrix
# Read the matrix
conn <- argsL$conn
M <- as.matrix(read.table(conn, sep = " ",header = FALSE))
print(paste0("INFO.... Connectome dimensions: ",dim(M)[1]," x ",dim(M)[2]))
if (dim(M)[1]==length(indx)) { print("Connectome is already sliced!"); # exit()
} else {
  print(paste0("INFO.... Connectome new dimensions: ",length(indx)," x ",length(indx)))
  M <- M[indx, indx]
  # Overwrites the connectome sliced
  print(paste0("    Path=",conn))
  write.table(M, conn, sep = " ", row.names = FALSE, col.names = FALSE)
   }

# Loads the MICs colormaps
mica.dir <- paste0(argsL$mica,"/functions/")
load(file=paste0(mica.dir,"cmap_MICs.Rdata"))

# Save png for QC purposes
nom <-  gsub(".txt","", strsplit(conn, split = "connectomes/")[[1]][2])
qc <- paste0(strsplit(conn, split = "proc_dwi/")[[1]][1],"QC/png/")
print(paste0("INFO.... Saving Connectome as png for QC"))
print(paste0("    Path=",qc, nom, ".png"))
png(paste0(qc, nom, ".png"))
  # Mirror matrix completes inferior triangle
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  image(log(M), axes=FALSE, main=nom, col=cmap.SC(256))
dev.off()
