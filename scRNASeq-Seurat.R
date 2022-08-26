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

# Load data ---------------------------------------------------------------

##load seurat object
try(setwd(workingdir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(projectName, "_sessionInfo.txt"))

raw.seurat.object <- Read10X_h5(paste0(path_to_data,"/outs/filtered_feature_bc_matrix.h5"))
raw.seurat.object<- CreateSeuratObject(raw.seurat.object, project = projectName, min.cells = min.cells, min.features = min.features)
raw.seurat.object
##add metadata


# QC ----------------------------------------------------------------------

raw.seurat.object[["percent.mt"]] <- PercentageFeatureSet(raw.seurat.object, pattern = "MT-")

##plots
VlnPlot(raw.seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
plot1<- FeatureScatter(raw.seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw.seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##filter by feature count and perent.mito
seurat.object <- subset(raw.seurat.object, 
												subset = nFeature_RNA >= 200 &
												nFeature_RNA <= 2500 &
												percent.mt <= 5)
filtered.cells <- length(rownames(raw.seurat.object@meta.data))-length(rownames(seurat.object@meta.data))
print(paste("filtered out", filtered.cells, "cells"))
seurat.object

VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalize and scale data ----------------------------------------------------------

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
seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress = "nCount_RNA")


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
DimPlot(seurat.object, reduction = "umap")

##saveRDS
saveRDS(seurat.object, file = paste0(workingdir, projectName, "-seuratObject.RDS"))




# Add metadata ------------------------------------------------------------

metadata.df <- read.csv(file = "PancDB_data/20220801MetaData.csv", header = TRUE, row.names = 1)



mydf <- seurat.object@meta.data
head(mydf)
head(pt.data)
dim(pt.data)


for(col.pos in 1:dim(metadata.df)[2]){
	col.id <- colnames(metadata.df)[col.pos]
	for(row.id in 1:dim(metadata.df)[1]){
		metadata.info <- metadata.df[row.id, col.id]
		simplified.pt.id <- stringr::str_replace_all(rownames(metadata.df)[row.id], "[^[:alnum:]]", "")
		simplified.meta.id <- stringr::str_replace_all(seurat.object@meta.data$orig.ident, "[^[:alnum:]]", "")
		
		seurat.object@meta.data[[col.id]][grepl(simplified.pt.id, simplified.meta.id, ignore.case = TRUE) == TRUE] <- metadata.info
	}
}









seurat.object <- readRDS("hpap108-seuratObject.RDS")

seurat.object <- AssignMetadata(pt.data, seurat.object = seurat.object)
head(seurat.object)

for(datacol in 1:dim(pt.data)[2]){
	newcol <- colnames(pt.data)[datacol]
	for(pt.id in 1:dim(pt.data)[1]){
		pt.info <- pt.data[pt.id, newcol]
		simplified.pt.id <- stringr::str_replace_all(rownames(pt.data)[pt.id], "[^[:alnum:]]", "")
		simplified.meta.id <- stringr::str_replace_all(mydf$orig.ident, "[^[:alnum:]]", "")
			
		mydf[[newcol]][grepl(simplified.pt.id, simplified.meta.id, ignore.case = TRUE) == TRUE] <- pt.info
	}
}








