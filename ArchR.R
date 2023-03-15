
# Set up ------------------------------------------------------------------

projectName <- "Obesity_scHPAP"
path_to_data <- c(list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams/", full.names = TRUE, recursive = FALSE))
working.dir <- "/Users/heustonef/Desktop/Obesity/snATAC"
records.dir <- "~/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/"
nThreads <- parallelly::availableCores() - 4
res <- 0.5
testable.factors <- c("BMI", "obesity") # factors to query during Harmony regression


# Libraries ---------------------------------------------------------------

library(ArchR)
library(harmony)
addArchRGenome("hg38")
addArchRThreads(threads = nThreads)


# Load data ---------------------------------------------------------------


setwd(working.dir)
sink(paste0(records.dir, projectName, "_sessionInfo.txt"))
sessionInfo()
sink()



# Generate arrowFiles -----------------------------------------------------

for(i in path_to_data){
	ifelse(file.exists(paste0(i, "/outs")),"", path_to_data <-path_to_data[!path_to_data %in% i])
}
names(path_to_data) <- sapply(path_to_data, basename)

arrowfiles <- createArrowFiles(
	inputFiles = paste0(path_to_data, "/outs/fragments.tsv.gz"),
	sampleNames = names(path_to_data),
	filterTSS = 4, 
	filterFrags = 1000, 
	addTileMat = TRUE,
	addGeneScoreMat = TRUE
)
saveRDS(arrowfiles, paste0(projectName, "-ArrowFiles.RDS"))



# Subset Arrow Files ------------------------------------------------------

arrowfiles <- readRDS("ArchRFiles-ArrowFiles.RDS")
arrowfiles <- sort(arrowfiles)
metadata <- read.table(file = paste0(records.dir, "HPAPMetaData.txt"), header = TRUE, sep = "\t")
rownames(metadata) <- metadata$DonorID

# Define subset
# as of 2023.03.10 we are taking all "NoDM" donors regardless of BMI
# > only a little annoyed because it took me far too long to figure out the logic for the original cutoffs

for(i in arrowfiles){
	x <- strsplit(i, "_")[[1]][1]
	if(!(metadata[x, "SimpDisease"] == "NoDM" & metadata[x, "scATAC"] >0)){
		arrowfiles <- arrowfiles[arrowfiles!=i]}
}
arrowfiles <- paste0("ArrowFiles/", arrowfiles)

# Identify doublets -------------------------------------------------------


dbltScores <- addDoubletScores(input = (arrowfiles), k = 10, knnMethod = "UMAP", LSIMethod = 1) #25 samples = 40min


# Create arch.proj -----------------------------------------------------

arch.proj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = TRUE)
arch.proj
paste0("Memory Size = ", rounds(object.size(arch.proj)/10^6, 3), " MB")
getAvailableMatrices(arch.proj)

cellcoldata <- getCellColData(arch.proj, select = c("log10(nFrags)", "TSSEnrichment"))

# plot unique nuclear fragments vs TSS enrichment score
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

# -- TSS plots show some samples have 2 TSS enrichment at 2 different scores. Will set cutoff @ 8 to maintain a single enriched peak.

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
# arch.proj$obesity <- NA
# arch.proj$obesity[arch.proj$BMI >=30] <- 'obese'
# arch.proj$obesity[arch.proj$BMI <= 25] <- 'nonobese'
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
# loadArchRProject(path = working.dir)
# saveRDS(arch.proj, paste0(projectName, "_noFilters.RDS"))

# Filter  ---------------------------------------------------------

arch.proj <- filterDoublets(ArchRProj = arch.proj, cutEnrich = 1, filterRatio = 1.5) # see notes
arch.proj <- arch.proj[which(arch.proj$TSSEnrichment > 8 & arch.proj$nFrags > 1000 & arch.proj$nFrags<40000)]
	# >1000 : https://www.nature.com/articles/s41586-021-03604-1#Sec9

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


harmony.groups <- colnames(arch.proj@cellColData)
harmony.groups <- harmony.groups[!(harmony.groups %in% testable.factors)]
arch.proj <- addHarmony(ArchRProj = arch.proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = , force = TRUE) # addHarmony "grouby" defines variables to correct for

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
saveRDS(arch.proj, paste0(projectName, "_Harmony.RDS"))


arch.proj <- addUMAP(ArchRProj = arch.proj,
										 reducedDims = "Harmony",
										 name = "UMAP_harmony",
										 nNeighbors = 30,
										 minDist = 0.5,
										 metric = "cosine",
										 force = TRUE)
arch.proj <- addClusters(input = arch.proj, reducedDims = "Harmony", method = "Seurat", name = paste0("Harmony_res", as.character(res)), resolution = res, force = TRUE)


plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "obesity", embedding = "UMAP_harmony")
p1 <- plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "obesity", embedding = "UMAP_harmony", randomize = TRUE)
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

#Without MAGIC
markergenes <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markergenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markergenes, file = paste0(projectName, "-markergenes.RDS"))
markerList$C1

#With MAGIC
arch.proj <- addImputeWeights(ArchRProj = arch.proj, reducedDims = "Harmony")
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


# Calling peaks -----------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)
pathToMacs2 <- findMacs2()

arch.proj <- addGroupCoverages(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), force = TRUE)
arch.proj <- addReproduciblePeakSet(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), pathToMacs2 = pathToMacs2, )

arch.proj <- addPeakMatrix(arch.proj)


markerPeaks <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(markerPeaks, file = paste0(projectName, "-MarkerPeaks.RDS"))
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
heatmapPeaks <- markerHeatmap(seMarker = markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)


arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", name = "Motif")
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


arch.proj <- addArchRAnnotations(ArchRProj = arch.proj, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

enrichATAC <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "ATAC", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapATAC <- plotEnrichmentHeatmap(enrichATAC, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")




# CUSTOM ENRICHMENT!!!! ---------------------------------------------------

# https://www.archrproject.com/bookdown/custom-enrichment.html


# Motif Deviations --------------------------------------------------------

if("Motif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", name = "Motif")
}
arch.proj <- addBgdPeaks(arch.proj)
arch.proj <- addDeviationsMatrix(arch.proj, peakAnnotation = "Motif", force = TRUE)
plotVarDev <- getVarDeviations(arch.proj, name = "MotifMatrix", plot = TRUE)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)



# Footprinting ------------------------------------------------------------

motifPositions <- getPositions(arch.proj)


# Deviant Motifs ----------------------------------------------------------


seGroupMotif <- getGroupSE(arch.proj, useMatrix = "MotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
saveRDS(seGroupMotif, file = paste0(projectName, "_seGroupMotif.RDS"))
corGSM_MM <- correlateMatrices(arch.proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "Harmony")
saveRDS(corGSM_MM, file = paste0(projectName, "_corGSM_MM.RDS"))

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


corGSM_MM <- readRDS(paste0(projectName, "_corGSM_MM.RDS"))




# Integrate scRNA object (Seurat) -----------------------------------------

seurat.object <- readRDS("/Users/heustonef/Desktop/Obesity/")
#check import
colnames(seurat.object@meta.data)
seurat.object$SCT_snn_res.0.5 <- paste0("SCT", seurat.object$SCT_snn_res.0.5)

arch.proj <- addGeneIntegrationMatrix(  # step takes ~95min
	ArchRProj = arch.proj,
	useMatrix = "GeneScoreMatrix",
	matrixName = "GeneIntegrationMatrix", 
	reducedDims = "Harmony",
	seRNA = seurat.object,
	addToArrow = FALSE,
	groupRNA = "SCT_snn_res.0.5",
	nameCell = "predictedCell_Un",
	nameGroup = "predictedGroup_Un",
	nameScore = "predictedScore_Un"
	
)

saveArchRProject(arch.proj, outputDirectory = working.dir, load = TRUE)
pal <- paletteDiscrete(values = seurat.object$SCT_snn_res.0.5)
plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un", pal = pal)

cM <- as.matrix(confusionMatrix(arch.proj$Harmony_res0.5, arch.proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

# Trajectory --------------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)

pal <- paletteDiscrete(values = seurat.object$SCT_snn_res.0.5)

plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un", pal = pal) +
theme_ArchR(legendTextSize = 10)





panc.markers <- list()
panc.markers$alpha <- c(
	"GCG", 
	"TTR", 
	"IRX2", 
	"HIGD1A", 
	"GLS", 
	"TM4SF4",
	"FAP", 
	"GPX3", 
	"SLC7A2", 
	"GC")
panc.markers$prog.beta <- c(
	'MAFB'
)
panc.markers$beta <- c(
	"INS",
	"IAPP", 
	"G6PC2", 
	"ADCYAP1", 
	"ERO1B", 
	"DLK1", 
	"NPTX2", 
	"GSN", 
	"INS-IGF2", 
	"HADH",
	"MAFA")
panc.markers$gamma <- c(
	"PPY", 
	"ID2", 
	"GCNT3", 
	"FXYD2", 
	"STMN2", 
	"THSD7A", 
	"SLITRK6", 
	"SERTM1", 
	"TM4SF4", 
	"ETV1")
panc.markers$delta <- c(
	"SST", 
	"RBP4", 
	"LEPR", 
	"RGS2", 
	"SEC11C", 
	"PRG4", 
	"BCHE", 
	"ADGRL2", 
	"HHEX", 
	"SLC38A1")
panc.markers$prog.acinar <- c(
	'PTF1A',
	'NR5A2'
)
panc.markers$acinar <- c(
	'PTF1A',
	'AMY1A',
	'CPA1'
)
panc.markers$definitive.endoderm <- c(
	'FOX2A',
	'SOX17',
	'GATA4',
	'CXCR4',
	'KIT'
)
panc.markers$panc.endoderm <- c(
	'PDX1',
	'PTF1A',
	'NXK2-2',
	'NGN3',
	'IA1',
	'ISL1',
	'PAX6',
	"MNX1",
	"ONECUT1"
)
panc.markers$prog.duct.panc <- c(
	'SOX9',
	'NKX6-1',
	'CHGA'
)
panc.markers$duct <- c(
	'CK19',
	'CFTR',
	'HNF1B'
)
panc.markers$prog.endocrine <- c(
	'NGN3',
	'NEUROG3',
	'NEUROD1'
)

# Heatmaps ----------------------------------------------------------------

arch.markers <- readRDS("ArchRFiles-markergenes.RDS")
heatmap.islets <- plotMarkerHeatmap(seMarker = arch.markers, 
																		cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																		labelMarkers = unlist(panc.markers),
																		transpose = TRUE)


heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = paste0(projectName, "-UMAP_harmony-res", as.character(res), "-AllMarkersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
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
	
	png(filename = paste0(projectName, "-UMAP_harmony-res", as.character(res), "-", as.character(chart.name), "Markersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
	plot(heatmap.plot)
	dev.off()

	
}
