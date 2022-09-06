# Run analysis on single nucleus ATAC data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Works through Signac

# Load libraries ----------------------------------------------------------

library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

# Global parameters -------------------------------------------------------

## frequently modified
projectName <- ""
workingdir <- "./"
regression.param <- 0
cum.var.thresh <- 80
resolution <- 0.5

## infrequently modified
path_to_data <- "PancDB_data/cellranger_scRNA/HPAP-108"
run.jackstraw <- TRUE
min.cells <- 3
min.features <- 200

##load local functions
sourceable.functions <- list.files(path = "RFunctions", pattern = "*.R", full.names = TRUE)
invisible(sapply(sourceable.functions, source))
