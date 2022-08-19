# Run analysis on single cell RNA data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Patient HPAP-051 
# Demographics: female, T2D, BMI >31
# Available sequencing data: single cell [RNA_rep1, ATAC_alpha/beta(?)]; bulk [RNA, ATAC]


# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)

# Global parameters -------------------------------------------------------

## frequently modified
projectName <- "hpap051_scRNA"
workingdir <- "./"

## infrequently modified
min.cells <- 3
min.features <- 200

# Load data ---------------------------------------------------------------

try(setwd(workingdir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(projectName, "_sessionInfo.txt"))

seurat.object <- Read10X_h5("PancDB_data/HPAP-051_scRNA/outs/filtered_feature_bc_matrix.h5")
seurat.object<- CreateSeuratObject(seurat.object, project = projectName, min.cells = min.cells, min.features = min.features)
seurat.object


# QC ----------------------------------------------------------------------

seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "MT-")

#filter by feature count and perent.mito
seurat.object <- subset(seurat.object, 
												subset = nFeature_RNA >= 200 &
												nFeature_RNA <= 2500 &
												percent.mt <= 5)

#plots
VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2








