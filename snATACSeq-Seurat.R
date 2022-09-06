# Run analysis on single cell RNA data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Available sequencing data: single cell [RNA_rep1, ATAC_alpha/beta(?)]; bulk [RNA, ATAC]

# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

# Global parameters -------------------------------------------------------

## frequently modified
projectName <- "hpap108"
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
