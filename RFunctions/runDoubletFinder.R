# Notes -------------------------------------------------------------------

# script to run doublet finder
# based on tutorial from https://www.youtube.com/watch?v=NqvAS4HgmrE & https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_doublets.R <Biomagician>


# Best practices ----------------------------------------------------------

# Run on filtered data
# Run on each sample separately
# Samples need to be normalized and plotted

# DoubletFinder -----------------------------------------------------------

library(DoubletFinder)
library(ggplot2)

runDoubletFinder <- function(seurat.object = NULL, sctransformed = NULL, tot.var = NULL, predicted.doubletRate = 0.05){
	
	print("running DoubletFinder")
	#Run PCA	
	seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

	# determine dimensionality
	tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
	print(paste0("Num pcs for 80% variance: ", length(which(cumsum(tot.var) <= 80))))
	print(paste0("Num pcs for 85% variance: ", length(which(cumsum(tot.var) <= 85))))
	print(paste0("Num pcs for 90% variance: ", length(which(cumsum(tot.var) <= 90))))
	print(paste0("Num pcs for 95% variance: ", length(which(cumsum(tot.var) <= 95))))
	
	cluster.dims <- 0
	if(cum.var.thresh > 0){
		cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
	}
	
	seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims)
	seurat.object <- FindClusters(seurat.object, resolution = resolution)
	seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims)

	## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
	sweep.res.list <- paramSweep_v3(seurat.object, PCs = 1:cluster.dims, sct = sctransformed)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats)
	
	ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
		geom_point() +
		geom_line()
	
	pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
		filter(BCmetric == max(BCmetric)) %>%
		select(pK) 
	pK <- as.numeric(as.character(pK[[1]]))
	print(paste0("Found pK = ", pK))
	
	## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
	annotations <- seurat.object@meta.data$seurat_clusters
	homotypic.prop <- modelHomotypic(annotations) 
	nExp_poi <- round(nrow(seurat.object@meta.data) * predicted.doubletRate) 
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seurat.object <- doubletFinder_v3(seurat.object,
																	PCs = 1:cluster.dims, 
																	pK = pK, 
																	nExp = nExp_poi, 
																	reuse.pANN = FALSE,
																	sct = sctransformed)
names(seurat.object@meta.data)[grep("DF.cl", names(seurat.object@meta.data))] <- "DF.classifications"
names(seurat.object@meta.data)[grep("pANN", names(seurat.object@meta.data))] <- "pANN"

return(seurat.object)
}
