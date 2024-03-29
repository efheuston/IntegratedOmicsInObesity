# Notes
# Need to exclude HPAP-043. Has 32 cells

# Set up ------------------------------------------------------------------

atacProject <- "Obesity_snATAC-ArchR_NW-OB"
nThreads <- parallelly::availableCores()
res <- 0.5
testable.factors <- c("BMI", "obesity") # factors to query during Harmony regression
comp.type <- "macbookPro" # one of macbookPro, biowulf, or workPC

# Directories -------------------------------------------------------------

if(comp.type == "macbookPro"){
	working.dir <- "/Users/heustonef/Desktop/Obesity/snATAC/"
	path_to_data <- c(list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams/", full.names = TRUE, recursive = FALSE))
	records.dir <- "~/OneDrive/SingleCellMetaAnalysis/GitRepositories/PancObesity/"
	metadata.location <- "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/"
	functions.path <- "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/GitRepositories/RFunctions/"
} else if(comp.type == "biowulf"){
	working.dir <- "/data/CRGGH/heustonef/hpapdata/cellranger_snATAC"
	records.dir <- working.dir
	path_to_data <- c(list.dirs("/data/CRGGH/heustonef/hpapdata/cellranger_snATAC/cellrangerOuts/", full.names = TRUE, recursive = FALSE))
	metadata.location <- "/data/CRGGH/heustonef/hpapdata/"
	funtions.path <- "/data/CRGGH/heustonef/hpapdata/RFunctions/"
	# library(vctrs, lib.loc = "/data/heustonef/Rlib_local/")
	# library(purrr, lib.loc = "/data/heustonef/Rlib_local/")
}


# Libraries ---------------------------------------------------------------

library(ArchR)
library(harmony)
addArchRGenome("hg38")
addArchRThreads(threads = nThreads)


# Load data ---------------------------------------------------------------


setwd(working.dir)
sink(paste0(records.dir, atacProject, "_sessionInfo.txt"))
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
	minTSS = 4, 
	minFrags = 1000, 
	addTileMat = TRUE,
	addGeneScoreMat = TRUE, force = TRUE
)
saveRDS(arrowfiles, paste0(atacProject, "-ArrowFiles.RDS"))



# Subset Arrow Files ------------------------------------------------------

arrowfiles <- readRDS("ArchRFiles-ArrowFiles.RDS")
arrowfiles <- sort(arrowfiles)
# metadata <- read.table(file = paste0(records.dir, "HPAPMetaData.txt"), header = TRUE, sep = "\t")
metadata <- read.table(file = "/data/CRGGH/heustonef/hpapdata/cellranger_snATAC/HPAPMetaData.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$DonorID

# Define subset
# As of 4.3.23 were are going to limit to obese and lean only
# At this point exclude HPAP-092

for(i in arrowfiles){
	x <- strsplit(i, "_")[[1]][1]
	if(!(metadata[x, "SimpDisease"] == "NoDM" & metadata[x, "scATAC"] >0)){
		arrowfiles <- arrowfiles[arrowfiles!=i]}
	arrowfiles <- arrowfiles[!(arrowfiles %in% c("HPAP-092_FGC2061.arrow", "HPAP-092_FGC2381.arrow"))]
}

# Identify doublets -------------------------------------------------------
dbltScores <- addDoubletScores(input = (arrowfiles), k = 10, knnMethod = "UMAP", LSIMethod = 1) #25 samples = 40min


# Create arch.proj -----------------------------------------------------

arch.proj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = FALSE)
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
# saveRDS(arch.proj, paste0(atacProject, "_noFilters.RDS"))

# Filter  ---------------------------------------------------------

arch.proj <- filterDoublets(ArchRProj = arch.proj, cutEnrich = 1, filterRatio = 1.5) # see notes
arch.proj <- arch.proj[which(arch.proj$TSSEnrichment > 8 & 
														 	arch.proj$nFrags > 1000 & 
														 	arch.proj$nFrags<40000 &
														 	arch.proj$BlacklistRatio < 0.03)]
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


# harmony.groups <- colnames(arch.proj@cellColData)
# harmony.groups <- harmony.groups[!(harmony.groups %in% testable.factors)]
arch.proj@cellColData <- arch.proj@cellColData[,!names(arch.proj@cellColData) %in% c("TissueSource", "scRNA", "scATAC", "scMultiome", "BulkRNA", "BulkATAC")]

# factorize regression columns
arch.proj$SampleAge <- as.factor(arch.proj$SampleAge)
arch.proj <- addHarmony(ArchRProj = arch.proj, 
												reducedDims = "IterativeLSI", 
												name = "Harmony", 
												groupBy = c("Sample", "SampleSex", "SampleAge"), 
												max.iter.harmony = 20, #did not converge after 10
												force = TRUE) # addHarmony "grouby" defines variables to correct for

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
saveRDS(arch.proj, paste0(atacProject, "_Harmony.RDS"))


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

#Without MAGIC
markergenes <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markergenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markergenes, file = paste0(atacProject, "-markergenes.RDS"))
markerList$C1

#With MAGIC
# arch.proj.magic <- addImputeWeights(ArchRProj = arch.proj, reducedDims = "Harmony")
# saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


plotEmbedding(
	ArchRProj = arch.proj, 
	colorBy = "GeneScoreMatrix", 
	name = "PECAM1", 
	plotAs = "points",
	embedding = "UMAP_harmony",
	quantCut = c(0.01, 0.95),
	imputeWeights = NULL)

plotEmbedding(
	ArchRProj = arch.proj, 
	colorBy = "MotifMatrix", 
	name = "z:STAT3_777", 
	plotAs = "points",
	embedding = "UMAP_harmony",
	quantCut = c(0.01, 0.95),
	imputeWeights = getImputeWeights(arch.proj))

for(i in names(markerList)){
	df <- data.frame(markerList[i])
	df <- df %>% select(5, 7, 8, 9)
	write.table(x = df, file = paste0("snATAC_", i, "-markerGeneEnrichment.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
}


# Calling peaks -----------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)
library(BSgenome.Hsapiens.UCSC.hg38)
pathToMacs2 <- findMacs2()
BSgenome.Hsapiens.UCSC.hg18 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

arch.proj <- addGroupCoverages(arch.proj, groupBy = paste0("Harmony_res", as.character(res)))
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
arch.proj <- addReproduciblePeakSet(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), pathToMacs2 = pathToMacs2)

arch.proj <- addPeakMatrix(arch.proj)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


markerPeaks <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(markerPeaks, file = paste0(atacProject, "-MarkerPeaks.RDS"))

markerPeaks <- readRDS(paste0(atacProject, "-MarkerPeaks.RDS"))

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
heatmapPeaks <- plotMarkerHeatmap(seMarker = markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)

arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", name = "cisbp", force = TRUE)
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", name = "encode", force = TRUE) # 306 min
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "homer", name = "homer", force = TRUE)

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)



arch.proj <- addArchRAnnotations(ArchRProj = arch.proj, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#enrichATAC <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "ATAC", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
#heatmapATAC <- plotEnrichmentHeatmap(enrichATAC, n = 7, transpose = TRUE)
#ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")




# CUSTOM ENRICHMENT!!!! ---------------------------------------------------

# https://www.archrproject.com/bookdown/custom-enrichment.html


# Motif Deviations --------------------------------------------------------

if("Motif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", name = "Motif")
}
arch.proj <- addBgdPeaks(arch.proj)
arch.proj <- addDeviationsMatrix(arch.proj, peakAnnotation = "Motif", force = TRUE, )
plotVarDev <- getVarDeviations(arch.proj, name = "MotifMatrix", plot = TRUE)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)

## subset archR project
for(i in c("C2", "C3", "C4")){
	arch.cell.subset <- BiocGenerics::which(arch.proj$Harmony_res0.5 %in% i)
	cellsSample <- arch.proj$cellNames[arch.cell.subset]
	arch.subset <- arch.proj[cellsSample, ]
	
	arch.subset <- addBgdPeaks(arch.subset)
	arch.subset <- addDeviationsMatrix(arch.subset, peakAnnotation = "Motif", force = TRUE, )
	plotVarDev <- getVarDeviations(arch.subset, name = "MotifMatrix", plot = TRUE)
	saveArchRProject(ArchRProj = arch.subset, outputDirectory = "/Users/heustonef/Desktop/Obesity/snATAC/varmotifsubset/", load = TRUE)
	
	motifPositions <- getPositions(arch.subset)
	
	# Deviant Motifs ----------------------------------------------------------

	seGroupMotif <- getGroupSE(arch.subset, useMatrix = "MotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
	saveRDS(seGroupMotif, file = paste0(atacProject, "_seGroupMotif.RDS"))
	corGSM_MM <- correlateMatrices(arch.subset, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "Harmony")
	saveRDS(corGSM_MM, file = paste0(atacProject, "_corGSM_MM.RDS"))
	
	saveArchRProject(ArchRProj = arch.subset, outputDirectory = "/Users/heustonef/Desktop/Obesity/snATAC/varmotifsubset/", load = TRUE)
	# corGSM_MM <- readRDS(paste0(atacProject, "_C", as.character(i), "-corGSM_MM.RDS"))
}



# Footprinting ------------------------------------------------------------

motifPositions <- getPositions(arch.proj)


# Deviant Motifs ----------------------------------------------------------


seGroupMotif <- getGroupSE(arch.proj, useMatrix = "MotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
saveRDS(seGroupMotif, file = paste0(atacProject, "_seGroupMotif.RDS"))
corGSM_MM <- correlateMatrices(arch.proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "Harmony")
saveRDS(corGSM_MM, file = paste0(atacProject, "_corGSM_MM.RDS"))

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


corGSM_MM <- readRDS(paste0(atacProject, "_corGSM_MM.RDS"))




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

png(filename = paste0(atacProject, "-UMAP_harmony-res", as.character(res), "-AllMarkersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
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
	
	png(filename = paste0(atacProject, "-UMAP_harmony-res", as.character(res), "-", as.character(chart.name), "Markersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
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





# Export for UCSC browser -------------------------------------------------

#export big wig files
getGroupBW(
	ArchRProj = arch.subset,
	groupBy = "Harmony_res0.5",
	normMethod = "ReadsInTSS", 
	tileSize = 100, 
	maxCells = 1000,
	ceiling = 4,
	verbose = TRUE, 
	threads = getArchRThreads(),
	logFile = createLogFile("getGroupBW")
)

getGroupBW(
	ArchRProj = arch.subset,
	groupBy = "Obesity",
	normMethod = "ReadsInTSS", 
	tileSize = 100, 
	maxCells = 1000,
	ceiling = 4,
	verbose = TRUE, 
	threads = getArchRThreads(),
	logFile = createLogFile("getGroupBW")
)

#export peak files
peak.lists <- list.files(path = "~/Desktop/Obesity/snATAC/PeakCalls", full.names = TRUE, pattern = ".rds", ignore.case = TRUE)
for(i in peak.lists){
	cluster <- str_split(basename(i), "-")[[1]][1]
	gr <- readRDS(i)
	df <- df <- data.frame(seqnames=seqnames(gr),
												 starts=start(gr)-1,
												 ends=end(gr),
												 names=c(rep(".", length(gr))),
												 scores=c(rep(".", length(gr))),
												 strands=strand(gr)
	)
	write.table(df, file=paste0("~/Desktop/Obesity/snATAC/", cluster, "_peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

#create track files

track.file <- paste0("UCSCBrowserTrackFile-", track.type, ".txt")
track.extension <- ".bw"
track.type <- "bigWig"
track.name <- c("C", "ATAC_") # pos[1] is the target, pos[2] is the replacement
track.descriptor <- "ATAC_bw"
track.url <- "https://hpc.nih.gov/~heustonef/"
file.list <- list.files(path = "~/Desktop/Obesity/snATAC/GroupBigWigs/Harmony_res0.5", full.names = FALSE, pattern = ".bw", ignore.case = TRUE)
track.list <-c()
for(i in file.list){
	track.groupID <- str_split(i, "-T")[[1]][1]
	track.command <- paste0("track type=", track.type, 
													" name=", gsub(track.name[1], track.name[2], track.groupID), 
													" description=", track.descriptor, 
													" bigDataUrl=", track.url, i)
	track.list <- c(track.list, track.command)
}
lapply(track.list, write, track.file, append=TRUE)
track.list


file.list <- list.files(path = "~/Desktop/Obesity/snATAC/PeakCalls/", full.names = FALSE, pattern = ".bed", ignore.case = FALSE, recursive = FALSE)
rack.extension <- ".bed"
track.type <- "BED"
track.name <- c("C", "peaks_") # pos[1] is the target, pos[2] is the replacement
track.descriptor <- "ATAC_bed"
track.url <- "https://hpc.nih.gov/~heustonef/"
track.file <- paste0("UCSCBrowserTrackFile-", track.type, ".txt")
track.list <-c()
for(i in file.list){
	track.groupID <- str_split(i, "-T")[[1]][1]
	track.command <- paste0("track type=", track.type, #Track type not required for BED (only for BED details) 
													" name=", gsub(track.name[1], track.name[2], track.groupID), 
													" description=", track.descriptor, 
													" bigDataUrl=", track.url, i)
	track.list <- c(track.list, track.command)
}
lapply(track.list, write, track.file, append=TRUE)
track.list


p <- plotEmbedding(
	ArchRProj = arch.proj, 
	colorBy = "MotifMatrix", 
	name = "z:STAT3_777", 
	embedding = "UMAP_harmony",
	imputeWeights = getImputeWeights(arch.proj)
)
plotVarDev <- getVarDeviations(arch.proj, name = "MotifMatrix", plot = TRUE, threads = 1)

