# starsolo_soupx.R
## R script to perform ambient RNA removal on STARsolo count matrices with SoupX

# Get command line arguments
args=(commandArgs(TRUE))


# Load libraries  & helper functions ----
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
# library(future) #TODO add multicore for Seurat w/ `future` for `getClusterIDs()`

suppressPackageStartupMessages(library(SoupX))

# write_sparse requirements
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(R.utils))

# source("~/DWM_utils/sc_utils/seurat_helpers/seutils.R")

SOLO_DIR = args[1] # Directory to either `Gene` or `GeneFull` matrix outputs
NCORES = args[2]

# Helper function(s) ----
## Quick function to find clusters
getClusterIDs <- function(
  toc, 
  verbose=F
){
  seu <- CreateSeuratObject(toc)
  seu <- seu %>%
    NormalizeData(verbose=verbose)%>%
    ScaleData(verbose=verbose) %>%
    FindVariableFeatures(verbose=verbose)%>%
    RunPCA(verbose=verbose) %>%
    FindNeighbors(verbose=verbose) %>%
    FindClusters(verbose=verbose)

  return(
    setNames(Idents(seu), Cells(seu))
  )
}

# Borrowed/adapted from the Marioni Lab, DropletUtils package (https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R)
#   (Had R version issues getting it to work as a dependency)
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
write_sparse <- function(
    path, # name of new directory
    x, # matrix to write as sparse
    barcodes=NULL, # cell IDs, colnames
    features=NULL, # gene IDs, rownames
    overwrite=F,
    verbose=F
){
  require(utils, quietly=T)
  require(Matrix, quietly=T)
  require(R.utils, quietly=T)
  
  if(!dir.exists(path)){
    dir.create(
      path,
      showWarnings=verbose,
      recursive = T
    )
  }
  
  if(is.null(barcodes)){
    barcodes=colnames(x)
  }
  if(is.null(features)){
    features=rownames(x)
  }

  # Remove old files if overwriting...
  if(overwrite){
    if(verbose){message("Overwriting old files if they exist...")}
    
    if(file.exists(paste0(path, "/matrix.mtx.gz"))){ # check matrix
      file.remove(paste0(path, "/matrix.mtx.gz")) 
      if(verbose){message("matrix file removed...")}
    }
    if(file.exists(paste0(path, "/matrix.mtx.gz.tmp"))){ # check matrix tmp file
      file.remove(paste0(path, "/matrix.mtx.gz.tmp")) 
      if(verbose){message("matrix tmp file removed...")}
    }
    if(file.exists(paste0(path, "/barcodes.tsv.gz"))){ # check barcodes
      file.remove(paste0(path, "/barcodes.tsv.gz"))
      if(verbose){message("barcodes file removed...")}
    }
    if(file.exists(paste0(path, "/features.tsv.gz"))){ # check features
      file.remove(paste0(path, "/features.tsv.gz"))
      if(verbose){message("features file removed...")}
    }
  }

  # gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)
  # gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))

  # Open file connections
  mhandle <- file.path(
    path, 
    "matrix.mtx"
  )
  bhandle <- gzfile(
    file.path(path, "barcodes.tsv.gz"), 
    open="w"
  )
  fhandle <- gzfile(
    file.path(path, "features.tsv.gz"),
    open="w"
  )
  on.exit({
    close(bhandle)
    close(fhandle)
  })
    
  # Write matrix...
  tryCatch(
    {
      if(verbose){message("Writing new matrix file to `",mhandle,"`")}
      writeMM(
        x, 
        file=mhandle
      )
    },
    error = function(e) {
      cat("Error: ", e$message, "\n")
    },
    warning = function(w) {
      cat("Warning: ", w$message, "\n")
    }
  )

  # Write barcodes...
  tryCatch(
    {
      if(verbose){message("Writing new barcodes file to `",bhandle,"`!")}
      write(
        barcodes,
        file = bhandle
      )
    },
    error = function(e) {
      cat("Error: ", e$message, "\n")
    },
    warning = function(w) {
      cat("Warning: ", w$message, "\n")
    }
  )

  # Write features...
  tryCatch(
    {
      if(verbose){message("Writing new features file to `",fhandle,"`!")}
      write.table(
        features, 
        file=fhandle, 
        row.names=FALSE, 
        col.names=FALSE, 
        quote=FALSE, 
        sep="\t"
      )
    },
    error = function(e) {
      cat("Error: ", e$message, "\n")
    },
    warning = function(w) {
      cat("Warning: ", w$message, "\n")
    }
  )
  
  # Annoyingly, writeMM doesn't take connection objects.
  gzip(mhandle)
  if(verbose){message("New matrix file gzipped!")}
  
  # return(NULL)
}

# Read in raw and filtered matrices ----
message(paste0("Reading in raw/filtered matrices...\n"))
tod = Seurat::Read10X(paste0(SOLO_DIR,'/raw')) #droplets
toc = Seurat::Read10X(paste0(SOLO_DIR,'/filtered')) # cells

#     SoupX ----
# https://github.com/constantAmateur/SoupX

#initialize soup objects from STARsolo outputs
soup <- SoupChannel(
  tod=tod, # table_of_droplets
  toc=toc  # table_of_counts
)

# set cluster IDs for cells before soup estimations
## quick preprocessing/clustering to get cluster IDs for cells
message("Preprocessing to get cluster IDs...\n")
tmp.clusters <- getClusterIDs(toc=soup$toc)

soup <- setClusters(
  soup,
  tmp.clusters
)

# Memory clean-up
rm(
  toc, 
  tod, 
  tmp.clusters
)
gc(verbose=F)

# Run soupx
message(paste0("Running soupx...\n"))
soup <- autoEstCont(soup)
adj.mat <- adjustCounts(soup)
message(paste0("Adjusted matrix has ", ncol(adj.mat)," cells and ", nrow(adj.mat), " features...\n"))

# Save adjusted matrices to disk ----
message(paste0("Saving adjusted count matrix to `",SOLO_DIR,"/soupx","`...\n"))
if(!is.null(adj.mat)){
  write_sparse(
    path=paste0(SOLO_DIR,"/soupx"), # name of new directory
    x=adj.mat,                      # matrix to write as sparse
    barcodes=NULL,                  # cell IDs, colnames
    features=NULL,                  # gene IDs, rownames
    overwrite=TRUE,
    verbose=TRUE
  )
}else{
  message(paste0("Adjusted matrix is NULL for `", SOLO_DIR,"`...\n"))
}

# save Rho values
# rhos <- list()
# for(i in 1:length(soup)){
#   rhos[[i]] <- mean(soup[[i]]$metaData$rho)
# }
# rhos <- do.call(rbind,rhos)
