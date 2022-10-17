
# Set up ------------------------------------------------------------------

# Starting new package tutorial on ArchR. This streamlines (?) scRNA/snatac integration
# Full tutorial: https://www.archrproject.com/bookdown/getting-set-up.html
# Brief tutorial: https://www.archrproject.com/articles/Articles/tutorial.html
# Note that loading `library(ArchR)` automatically sets default number of threads to availableWorkers()*0.5
# This can be changed manually with `addArchRThreads(threads = 8)`



# Samples -----------------------------------------------------------------

# Samples were selected from HPAP database. 
# Requirements: BMI > 30, Caucasian or Black, T2D or NoDisease, matched scRNA sample



# Set global variables ----------------------------------------------------

projectName <- "T2D_ArchR"
path_to_data <- list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams", full.names = TRUE, recursive = FALSE)
working.dir <- "/Users/heustonef/Desktop/PancDB_Data/ArchR/"
records.dir <- "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/"
res = 0.5
nthreads <- parallel::detectCores()

# Libraries ---------------------------------------------------------------

library(ArchR)
library(pheatmap)
addArchRThreads(threads = nthreads)
addArchRGenome("hg38")


# Load data ---------------------------------------------------------------


setwd(working.dir)
sink(paste0(records.dir, projectName, "_sessionInfo.txt"))
sessionInfo()
sink()

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

# In the future ADD METADATA HERE

saveRDS(arrowfiles, paste0(projectName, "-ArrowFiles.RDS"))


# Identify doublets -------------------------------------------------------

dbltScores <- addDoubletScores(input = arrowfiles, k = 10, knnMethod = "UMAP", LSIMethod = 1) #25 samples = 40min


# Create diab.project -----------------------------------------------------

diab.proj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = TRUE)
diab.proj
paste0("Memory Size = ", round(object.size(diab.proj)/10^6, 3), " MB")
getAvailableMatrices(diab.proj)

cellcoldata <- getCellColData(diab.proj, select = c("log10(nFrags)", "TSSEnrichment"))

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
	ArchRProj = diab.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
) +
	geom_vline(xintercept = 8, lty = "dashed") 

# -- TSS plots show some samples have 2 TSS enrichment at 2 different scores. Will set cutoff @ 8 to maintain a single enriched peak.

# plot nFrags ridge plot per sample
plotGroups(
	ArchRProj = diab.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges") +
	geom_vline(xintercept = log10(1000), lty = "dashed") +
	geom_vline(xintercept = log10(40000), lty = "dashed")

metadata <- read.table(file = "../HPAPMetaData-20220912.txt", header = TRUE, sep = "\t")

diab.proj@cellColData[,names(metadata)] <- lapply(names(metadata), function(x){
	diab.proj@cellColData[[x]] <- metadata[match(as.character(diab.proj$Sample), metadata$SampleID), x]
	}
)

saveArchRProject(ArchRProj = diab.proj, outputDirectory = working.dir, load = TRUE)
# loadArchRProject(path = working.dir)
saveRDS(diab.proj, paste0(projectName, "_noFilters.RDS"))
# diab.proj <- readRDS("T2D_ArchR_noFilters.RDS")

# Filter  ---------------------------------------------------------

diab.proj <- filterDoublets(ArchRProj = diab.proj, cutEnrich = 1, filterRatio = 1.5) # see notes
diab.proj <- diab.proj[which(diab.proj$TSSEnrichment > 8 & diab.proj$nFrags > 1000 & diab.proj$nFrags<40000)]
	# >1000 : https://www.nature.com/articles/s41586-021-03604-1#Sec9

diab.proj <- addIterativeLSI(
	ArchRProj = diab.proj,
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



# addHarmony "grouby" defines variables to correct for

saveArchRProject(ArchRProj = diab.proj, outputDirectory = working.dir, load = TRUE)

diab.proj <- addClusters(input = diab.proj, reducedDims = "IterativeLSI", method = "Seurat", name = paste0("LSI_res", as.character(res)), resolution = res, force = TRUE)

table(getCellColData(ArchRProj = diab.proj, select = paste0("LSI_res", as.character(res))))
cM <- confusionMatrix(paste0(diab.proj$LSI_res0.5), paste0(diab.proj$Sample)) # Could not automate this line
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM), 
	color = paletteContinuous("whiteBlue"), 
	border_color = "black"
)

cM <- confusionMatrix(paste0(diab.proj$LSI_res0.5), paste0(diab.proj$Simp.Disease))
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM), 
	color = paletteContinuous("whiteBlue"), 
	border_color = "black"
)
diab.proj <- addUMAP(ArchRProj = diab.proj,
										 reducedDims = "IterativeLSI",
										 name = "UMAP_LSI",
										 nNeighbors = 30,
										 minDist = 0.5, 
										 metric = "cosine")

plotEmbedding(ArchRProj = diab.proj, colorBy = "cellColData", name = "Simp.Disease", embedding = "UMAP_LSI")
saveArchRProject(ArchRProj = diab.proj, outputDirectory = working.dir, load = TRUE)

temp.proj <- diab.proj
diab.proj <- addHarmony(ArchRProj = diab.proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = c("SequencerID", "SampleEthnicity", "SampleAge", "BMI", "SampleSex"), force = TRUE)

temp.proj <- addHarmony(ArchRProj = temp.proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = c("SequencerID", "SampleEthnicity", "SampleAge", "BMI", "SampleSex"), force = TRUE)


diab.proj <- addUMAP(ArchRProj = diab.proj,
										 reducedDims = "Harmony",
										 name = "UMAP_harmony",
										 nNeighbors = 30,
										 minDist = 0.5,
										 metric = "cosine",
										 force = TRUE)
diab.proj <- addClusters(input = diab.proj, reducedDims = "Harmony", method = "Seurat", name = paste0("Harmony_res", as.character(res)), resolution = res, force = TRUE)


plotEmbedding(ArchRProj = diab.proj, colorBy = "cellColData", name = "Simp.Disease", embedding = "UMAP_harmony")
p1 <- plotEmbedding(ArchRProj = diab.proj, colorBy = "cellColData", name = "Simp.Disease", embedding = "UMAP_harmony", randomize = TRUE)
p2 <- plotEmbedding(ArchRProj = diab.proj, colorBy = "cellColData", name = paste0("Harmony_res", as.character(res)), embedding = "UMAP_harmony")
ggAlignPlots(p1, p2, type = "h")
png(filename = paste0(projectName, "-UMAP_harmony-res", as.character(res), ".png"), height = 800, width = 1200, bg = "transparent")
ggAlignPlots(p1, p2, type = "h")
dev.off()

plotEmbedding(ArchRProj = diab.proj, colorBy = "cellColData", name = "DiseaseStatus", embedding = "UMAP_harmony")




saveArchRProject(ArchRProj = diab.proj, outputDirectory = working.dir, load = TRUE)

# make sure devtools::install_github("immunogenomics/presto") is installed!
diab.markers <- getMarkerFeatures(ArchRProj = diab.proj,
																	useMatrix = "GeneScoreMatrix",
																	groupBy = "Harmony_res0.5",
																	bias = c("TSSEnrichment", "log10(nFrags)"),
																	testMethod = "wilcoxon",
																	threads = 1)
saveArchRProject(ArchRProj = diab.proj, outputDirectory = working.dir, load = TRUE)
saveRDS(diab.markers, paste0(projectName, "_HarmonyMarkerFeatures.RDS"))

diab.markerList <- getMarkers(diab.markers, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
diab.markerList


# Cell type markers -------------------------------------------------------

# load cell type markers. These are accessed from __Generation of human islet cell type-specific identity genesets, Nature Communications, 2022__.  
islet.markers <- list()
islet.markers$islets.alpha <- read.table(file = "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/Generation of human islet cell type-specific identity genesets GeneSetsDownload - 2022/herrera_alpha_hs.txt")[,1]
islet.markers$islets.beta <- read.table(file = "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/Generation of human islet cell type-specific identity genesets GeneSetsDownload - 2022/herrera_beta_hs.txt")[,1]
islet.markers$islets.gamma <- read.table(file = "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/Generation of human islet cell type-specific identity genesets GeneSetsDownload - 2022/herrera_gamma_hs.txt")[,1]
islet.markers$islets.delta <- read.table(file = "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/Generation of human islet cell type-specific identity genesets GeneSetsDownload - 2022/herrera_delta_hs.txt")[,1]

#Load top markers from same publication
source("/Users/heustonef/Desktop/PancDB_Data/ArchR/PancreaticCellMarkers.R")


# Heatmaps ----------------------------------------------------------------

heatmap.islets <- plotMarkerHeatmap(seMarker = diab.markers, 
																cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																labelMarkers = unlist(top.markers),
																transpose = TRUE)

heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = paste0(projectName, "-UMAP_harmony-res", as.character(res), "-AllMarkersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
plot(heatmap.plot)
dev.off()


length(intersect(unlist(diab.markerList)$name, unlist(top.markers))) # this tells us the number of "top.markers" found in the diab.markerList
intersect(unlist(diab.markerList)$name, unlist(top.markers)) # this tells us which "top markers" are in the diab.markerList













