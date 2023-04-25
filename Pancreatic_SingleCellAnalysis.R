# Set up ------------------------------------------------------------------

# Define variables for scATAC analysis in ArchR
scProject <- "Panc_SingleCellAnalysis"
nThreads <- parallelly::availableCores()
res <- 0.5

# Define variables for scRNA analysis in Seurat
regression.vars <- c("sequencerID", "SampleSex", "SampleAge")
cum.var.thresh <- 90
resolution <- 0.5
do.sctransform <- "each" # one of FALSE, each, pooled

## infrequently modified
do.doubletFinder <- TRUE
min.cells <- 3
min.features <- 200
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05



# Load libraries, functions, and metadata ---------------------------------------------------------------

library(ArchR)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)


## load user-supplied functions
invisible(sapply(list.files(path = "RFunctions/", pattern = "*.R", full.names = TRUE), source))

## load metadata
metadata <- read.table(file = "HPAPMetaData.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$DonorID # set donorID as rownaes


# Set run parameters
addArchRGenome("hg38")
addArchRThreads(threads = nThreads)

# Generate arrowFiles -----------------------------------------------------

# Create sample list
for(i in path_to_data){
	ifelse(file.exists(paste0(i, "/outs")),"", path_to_data <-path_to_data[!path_to_data %in% i])
}
names(path_to_data) <- sapply(path_to_data, basename)

# Create arrow files
arrowfiles <- createArrowFiles(
	inputFiles = paste0(path_to_data, "/outs/fragments.tsv.gz"),
	sampleNames = names(path_to_data),
	minTSS = 4, 
	minFrags = 1000, 
	addTileMat = TRUE,
	addGeneScoreMat = TRUE, force = TRUE
)
saveRDS(arrowfiles, paste0(scProject, "-ArrowFiles.RDS"))


# Identify doublets -------------------------------------------------------

dbltScores <- addDoubletScores(input = (arrowfiles), k = 10, knnMethod = "UMAP", LSIMethod = 1)

# Create arch.proj -----------------------------------------------------

arch.proj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = FALSE)

# scATAC QC Plots ----------------------------------------------------------------

# plot unique nuclear fragments vs TSS enrichment score
cellcoldata <- getCellColData(arch.proj, select = c("log10(nFrags)", "TSSEnrichment"))
ggPoint(
	x = cellcoldata[,1],
	y = cellcoldata[,2],
	colorDensity = TRUE,
	continuousSet = "sambaNight", 
	xlabel = "Log10 Unique Fragments",
	ylabel = "TSS Enrichment",
	xlim = c(log10(500), quantile(cellcoldata[,1], probs = 0.99)),
	ylim = c(0, quantile(cellcoldata[,2], probs = 0.99))) +
	geom_hline(yintercept = 4, lty = "dashed") +
	geom_vline(xintercept = 3, lty = "dashed")

# plot TSSEnrichment ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
) +
	geom_vline(xintercept = 8, lty = "dashed") 

# plot nFrags ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges") +
	geom_vline(xintercept = log10(1000), lty = "dashed") +
	geom_vline(xintercept = log10(30000), lty = "dashed")

arch.proj@cellColData[,names(metadata)] <- lapply(names(metadata), function(x){
	arch.proj@cellColData[[x]] <- metadata[match(vapply(strsplit(as.character(arch.proj$Sample), "_"), `[`, 1, FUN.VALUE = character(1)), metadata$DonorID), x]
	}
)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)

# Filter  ---------------------------------------------------------

arch.proj <- filterDoublets(ArchRProj = arch.proj, cutEnrich = 1, filterRatio = 1.5) 
arch.proj <- arch.proj[which(arch.proj$TSSEnrichment > 8 & 
														 	arch.proj$nFrags > 1000 & 
														 	arch.proj$nFrags<40000 &
														 	arch.proj$BlacklistRatio < 0.03)]

arch.proj <- addIterativeLSI(
	ArchRProj = arch.proj,
	useMatrix = "TileMatrix", 
	name = "IterativeLSI", 
	iterations = 10, 
	clusterParams = list( #See Seurat::FindClusters
		resolution = c(0.3), 
		sampleCells = 10000, 
		n.start = 10
	), 
	varFeatures = 25000, 
	dimsToUse = 1:30,
	force = TRUE)

# Run Harmony -------------------------------------------------------------

arch.proj@cellColData <- arch.proj@cellColData[,!names(arch.proj@cellColData) %in% c("TissueSource", "scRNA", "scATAC", "scMultiome", "BulkRNA", "BulkATAC")] # Remove extra columns from metadata

arch.proj$SampleAge <- as.factor(arch.proj$SampleAge) # factorize regression columns

arch.proj <- addHarmony(ArchRProj = arch.proj, 
												reducedDims = "IterativeLSI", 
												name = "Harmony", 
												groupBy = c("Sample", "SampleSex", "SampleAge"), 
												max.iter.harmony = 20,
												force = TRUE)

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)

arch.proj <- addUMAP(ArchRProj = arch.proj,
										 reducedDims = "Harmony",
										 name = "UMAP_harmony",
										 nNeighbors = 30,
										 minDist = 0.5,
										 metric = "cosine",
										 force = TRUE)

arch.proj <- addClusters(input = arch.proj, reducedDims = "Harmony", method = "Seurat", name = paste0("Harmony_res", as.character(res)), resolution = res, force = TRUE)

plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "Obesity", embedding = "UMAP_harmony")
plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "BMI", embedding = "UMAP_harmony", plotAs = "points")
p1 <- plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "Obesity", embedding = "UMAP_harmony", randomize = TRUE)
p2 <- plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = paste0("Harmony_res", as.character(res)), embedding = "UMAP_harmony")
ggAlignPlots(p1, p2, type = "h")

table(getCellColData(ArchRProj = arch.proj, select = paste0("Harmony_res", as.character(res))))
cM <- confusionMatrix(paste0(arch.proj$Harmony_res0.5), paste0(arch.proj$obesity)) # Could not automate this line
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)

cM <- confusionMatrix(paste0(arch.proj@cellColData[,paste0("Harmony_res", as.character(res))]), paste0(arch.proj$obesity))
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)



# Identify marker genes ---------------------------------------------------

markergenes <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markergenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markergenes, file = paste0(scProject, "-markergenes.RDS"))

# Calling peaks -----------------------------------------------------------

pathToMacs2 <- findMacs2()
BSgenome.Hsapiens.UCSC.hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

arch.proj <- addGroupCoverages(arch.proj, groupBy = paste0("Harmony_res", as.character(res)))
arch.proj <- addReproduciblePeakSet(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), pathToMacs2 = pathToMacs2)
arch.proj <- addPeakMatrix(arch.proj)


markerPeaks <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(markerPeaks, file = paste0(scProject, "-MarkerPeaks.RDS"))
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

heatmapPeaks <- plotMarkerHeatmap(seMarker = markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", name = "encode", force = TRUE)

enrichEncode <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)

# Motif Deviations --------------------------------------------------------

if("Motif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", name = "Motif")
}
arch.proj <- addBgdPeaks(arch.proj)
arch.proj <- addDeviationsMatrix(arch.proj, peakAnnotation = "Motif", force = TRUE, )
plotVarDev <- getVarDeviations(arch.proj, name = "MotifMatrix", plot = TRUE)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)



# Footprinting ------------------------------------------------------------

motifPositions <- getPositions(arch.proj)


# Deviant Motifs ----------------------------------------------------------


seGroupMotif <- getGroupSE(arch.proj, useMatrix = "MotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
saveRDS(seGroupMotif, file = paste0(scProject, "_seGroupMotif.RDS"))
corGSM_MM <- correlateMatrices(arch.proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "Harmony")
saveRDS(corGSM_MM, file = paste0(scProject, "_corGSM_MM.RDS"))

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


# Integrate scRNA object (Seurat) -----------------------------------------

seurat.object <- readRDS("/Users/heustonef/Desktop/Obesity/scRNA/Obesity_scRNA-Anchored-NW-OB-90pctvar.RDS")
#check import
colnames(seurat.object@meta.data)
seurat.object$integrated_snn_res.0.5 <- paste0("RNA", seurat.object$integrated_snn_res.0.5)

arch.proj <- addGeneIntegrationMatrix(  # step takes ~95min
	ArchRProj = arch.proj,
	useMatrix = "GeneScoreMatrix",
	matrixName = "GeneIntegrationMatrix", 
	reducedDims = "Harmony",
	seRNA = seurat.object,
	addToArrow = FALSE,
	groupRNA = "integrated_snn_res.0.5",
	nameCell = "predictedCell_Un",
	nameGroup = "predictedGroup_Un",
	nameScore = "predictedScore_Un"
	
)

saveArchRProject(arch.proj, outputDirectory = working.dir, load = TRUE)
pal <- paletteDiscrete(values = seurat.object$integrated_snn_res.0.5)
plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un", pal = pal)

cM <- as.matrix(confusionMatrix(arch.proj$Harmony_res0.5, arch.proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
saveArchRProject(arch.proj, outputDirectory = working.dir, load = TRUE)


# Trajectory --------------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)

pal <- paletteDiscrete(values = seurat.object$integrated_snn_res.0.5)

plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un") +
	theme_ArchR(legendTextSize = 10)
plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "Harmony_res0.5",) +
	theme_ArchR(legendTextSize = 10)



# Heatmaps ----------------------------------------------------------------

arch.markers <- readRDS(paste0(working.dir, "Obesity_scHPAP-markergenes.RDS"))
heatmap.islets <- plotMarkerHeatmap(seMarker = arch.markers, 
																		cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																		labelMarkers = unlist(panc.markers),
																		transpose = TRUE)


heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = paste0(scProject, "-UMAP_harmony-res", as.character(res), "-AllMarkersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
plot(heatmap.plot)
dev.off()

for(i in 1:length(panc.markers)){
	chart.name <- names(panc.markers[i])
	gene.set <- panc.markers[i]
	
	heatmap.islets <- plotMarkerHeatmap(seMarker = arch.markers, 
																			cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																			labelMarkers = unlist(gene.set),
																			transpose = TRUE)
	
	
	heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")
	
	png(filename = paste0(scProject, "-UMAP_harmony-res", as.character(res), "-", as.character(chart.name), "Markersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
	plot(heatmap.plot)
	dev.off()
	
	
}


motifPositions <- getPositions(arch.subset)

motifs <- c("FOS", "STAT3", "FOSL2", "NFE2L2", "MAFK", "FOXA1", "FOXA2")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs


seFoot <- getFootprints(
	ArchRProj = arch.subset, 
	positions = motifPositions[markerMotifs], 
	groupBy = "Harmony_res0.5"
)
plotFootprints(
	seFoot = seFoot,
	ArchRProj = arch.subset, 
	normMethod = "Subtract",
	plotName = "Footprints-Subtract-Bias",
	addDOC = FALSE,
	smoothWindow = 5
)







# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), paste0(rnaProject, "_sessionInfo.txt"))

