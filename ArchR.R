
# Set up ------------------------------------------------------------------

# Set global variables ----------------------------------------------------

projectName <- "Obesity_scHPAP"
path_to_data <- list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams", full.names = TRUE, recursive = FALSE)
working.dir <- "/Users/heustonef/Desktop/PancDB_Data/ArchR/"
records.dir <- "/Users/heustonef/OneDrive - National Institutes of Health/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/"

# Libraries ---------------------------------------------------------------

library(ArchR)
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


# Create archproject -----------------------------------------------------

archproj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = TRUE)
archproj
paste0("Memory Size = ", round(object.size(archproj)/10^6, 3), " MB")
getAvailableMatrices(archproj)

cellcoldata <- getCellColData(archproj, select = c("log10(nFrags)", "TSSEnrichment"))

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
	ArchRProj = archproj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
) +
	geom_vline(xintercept = 8, lty = "dashed") 

# -- TSS plots show some samples have 2 TSS enrichment at 2 different scores. Will set cutoff @ 8 to maintain a single enriched peak.

# plot nFrags ridge plot per sample
plotGroups(
	ArchRProj = archproj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges") +
	geom_vline(xintercept = log10(1000), lty = "dashed") +
	geom_vline(xintercept = log10(40000), lty = "dashed")

metadata <- read.table(file = "../HPAPMetaData-20220912.txt", header = TRUE, sep = "\t")

archproj@cellColData[,names(metadata)] <- lapply(names(metadata), function(x){
	coldata[[x]] <- metadata[match(as.character(coldata$Sample), metadata$SampleID), x]
	}
)


saveArchRProject(ArchRProj = archproj, outputDirectory = working.dir, load = TRUE)
# loadArchRProject(path = working.dir)

# Filter  ---------------------------------------------------------

archproj <- filterDoublets(ArchRProj = archproj, cutEnrich = 1, filterRatio = 1.5) # see notes
archproj <- archproj[which(archproj$TSSEnrichment > 8 & archproj$nFrags > 1000 & archproj$nFrags<40000)]
	# >1000 : https://www.nature.com/articles/s41586-021-03604-1#Sec9

archproj <- addIterativeLSI(
	ArchRProj = archproj,
	useMatrix = "TileMatrix", 
	name = "IterativeLSI", 
	iterations = 6, 
	clusterParams = list( #See Seurat::FindClusters
		resolution = c(0.5), 
		sampleCells = 10000, 
		n.start = 10
	), 
	varFeatures = 25000, 
	dimsToUse = 1:30,
	force = TRUE)


archproj <- addHarmony(ArchRProj = archproj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = c("Sample", "nFrags", "SampleID"), force = TRUE)

# addHarmony "grouby" defines variables to correct for

saveArchRProject(ArchRProj = archproj, outputDirectory = working.dir, load = TRUE)

archproj <- addClusters(input = archproj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 0.5, force = TRUE)

table(archproj$Clusters)
cM <- confusionMatrix(paste0(archproj$Clusters), paste0(archproj$Sample))
cM
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
	mat = as.matrix(cM), 
	color = paletteContinuous("whiteBlue"), 
	border_color = "black"
)
p
