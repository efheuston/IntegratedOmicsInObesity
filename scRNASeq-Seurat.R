# Run analysis on single cell RNA data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Patient HPAP-051 
# Demographics: female, T2D, BMI >31
# Available sequencing data: single cell [RNA_rep1, ATAC_alpha/beta(?)]; bulk [RNA, ATAC]


##NB: 
# currently only Scale data (not sctransform) is included. Add "regression.param <- not0" and sctransform equation to regress variables


# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

# Global parameters -------------------------------------------------------

## frequently modified
projectName <- "ObesityL_scRNA"
workingdir <- "./"
path_to_data <- "PancDB_data/cellranger_scRNA"
regression.vars <- NULL
cum.var.thresh <- 80
resolution <- 0.5

## infrequently modified
do.sctransform <- TRUE
run.jackstraw <- TRUE
min.cells <- 3
min.features <- 200

##load local functions
sourceable.functions <- list.files(path = "RFunctions", pattern = "*.R", full.names = TRUE)
invisible(sapply(sourceable.functions, source))




# Load data ---------------------------------------------------------------

try(setwd(workingdir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(projectName, "_sessionInfo.txt"))




##load seurat object

sc.data <- sapply(list.dirs(path = path_to_data, recursive = FALSE, full.names = TRUE), 
									basename, 
									USE.NAMES = TRUE)


object.list <- c()
for(i in 1:length(sc.data)){
	object.list[[i]] <- Read10X_h5(paste0(names(sc.data)[i], "/outs/filtered_feature_bc_matrix.h5"))
	object.list[[i]] <- CreateSeuratObject(object.list[[i]], 
																			 project = projectName, 
																			 min.cells = min.cells, 
																			 min.features = min.features)
	print(paste("finished", sc.data[[i]]))
}

seurat.object <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list), project = projectName)

##add metadata
metadata.df <- read.csv(file = "PancDB_data/20220801MetaData.csv", header = TRUE, row.names = 1)
seurat.object <- AssignMetadata(metadata.df = metadata.df, seurat.object = seurat.object) # created "AssignMetadata" function

# QC ----------------------------------------------------------------------

seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "MT-")

##filter by feature count and perent.mito
seurat.object <- subset(seurat.object, 
												subset = nFeature_RNA >= 200 &
												nFeature_RNA <= 2500 &
												percent.mt <= 5)
filtered.cells <- length(rownames(seurat.object@meta.data))-length(rownames(seurat.object@meta.data))
print(paste("filtered out", filtered.cells, "cells"))
seurat.object

##plot qc stats
VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalize and scale data ----------------------------------------------------------


if(do.sctransform == FALSE){ # standard method
	print("Performing standard normalization and scaling")
	if(length(regression.vars) >1){
		print("HEY YOU! You're performing standard scaling on more than 1 regression variable. You should probably be doing SCTransform. Set `do.sctransform` to TRUE")
	}

	##normalize
	seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
	
	##find HVG
	seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
	top10hvg <- head(VariableFeatures(seurat.object), 10)
	plot1 <- VariableFeaturePlot(seurat.object)
	plot2	<- LabelPoints(plot = plot1, points = top10hvg, repel = TRUE)
	plot1 + plot2
	top10hvg
	
	##scale (a linear transformation)
	all.genes <- rownames(seurat.object)
	
	seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress = regression.vars)
	
} else if(do.sctransform == TRUE){
	print("Performing SCTransform")
	seurat.object <- SCTransform(seurat.object, method = "glsGamPoi", vars.to.regress = regression.vars, verbose = TRUE)
} else {
	print("Must set do.sctransform to logical")
}

# Linear dimensional reduction --------------------------------------------

seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object), )

print(seurat.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat.object, dims = 1:2, reduction = "pca")
DimPlot(seurat.object, reduction = "pca")
DimHeatmap(seurat.object, dims = 1:2, cells = 500, balanced = TRUE)


# Determine dimensionality ------------------------------------------------

if(run.jackstraw == TRUE){
	seurat.object <- JackStraw(seurat.object, num.replicate = 100)
	seurat.object <- ScoreJackStraw(seurat.object, dims = 1:40)
	JackStrawPlot(seurat.object, dims = 1:40)
}

ElbowPlot(seurat.object)

# account for variance

tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
paste0("Num pcs for 80% variance: ", length(which(cumsum(tot.var) <= 80)))
paste0("Num pcs for 85% variance: ", length(which(cumsum(tot.var) <= 85)))
paste0("Num pcs for 90% variance: ", length(which(cumsum(tot.var) <= 90)))
paste0("Num pcs for 95% variance: ", length(which(cumsum(tot.var) <= 95)))

cluster.dims <- 0
if(cum.var.thresh > 0){
	cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
}


# Louvain cluster ---------------------------------------------------------

##cluster cells
seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims)
seurat.object <- FindClusters(seurat.object, resolution = resolution)

##umap
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims)
DimPlot(seurat.object, reduction = "umap", cols = color.palette, label = T, label.size = 7, repel = T)

##saveRDS
saveRDS(seurat.object, file = paste0(workingdir, projectName, "-seuratObject.RDS"))


# Find cluster biomarkers -------------------------------------------------

##find positively expressed markers for all clusters compared to all remaining clusters

markers.seurat.pos <- FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
testing <- markers.seurat.pos %>% 
	group_by(cluster) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)

markers.seurat.all <- FindAllMarkers(seurat.object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.all %>%
	group_by(cluster) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
	
##create workbook
markers.table <- openxlsx::createWorkbook()

##write positive markers to table
openxlsx::addWorksheet(markers.table, sheetName = "PosMarkers")
openxlsx::writeData(markers.table, sheet = "PosMarkers", x = markers.seurat.pos,startCol = 1, startRow = 1, colNames = TRUE)

##write all markers to table
openxlsx::addWorksheet(markers.table, sheetName = "AllMarkers")
openxlsx::writeData(markers.table, sheet = "AllMarkers", x = markers.seurat.all, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(projectName, "_seuratMarkers.xlsx"), overwrite = FALSE, returnValue = TRUE)

