
# Set up ------------------------------------------------------------------

# Starting new package tutorial on ArchR. This streamlines (?) scRNA/snatac integration
# Full tutorial: https://www.archrproject.com/bookdown/getting-set-up.html
# Brief tutorial: https://www.archrproject.com/articles/Articles/tutorial.html
# Note that loading `library(ArchR)` automatically sets default number of threads to availableWorkers()*0.5
# This can be changed manually with `addArchRThreads(threads = 8)`



# Set global variables ----------------------------------------------------

projectName <- "T2D_ArchR"
path_to_data <- list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams", full.names = TRUE, recursive = FALSE)
working.dir <- "/Users/heustonef/Desktop/PancDB_Data/ArchR/"


# Libraries ---------------------------------------------------------------

library(ArchR)
addArchRGenome("hg38")



# Load data ---------------------------------------------------------------


setwd(working.dir)
sink(paste0(projectName, "_sessionInfo.txt"))
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
saveRDS(arrowfiles, paste0(projectName, "-ArrowFiles.RDS"))


# Identify doublets -------------------------------------------------------

dbltScores <- addDoubletScores(input = arrowfiles, k = 10, knnMethod = "UMAP", LSIMethod = 1) #25 samples = 40min


# Create ArchRProject -----------------------------------------------------

archrproj <- ArchRProject(ArrowFiles = arrowfiles, outputDirectory = working.dir, copyArrows = TRUE)
archrproj
paste0("Memory Size = ", round(object.size(archrproj)/10^6, 3), " MB")
getAvailableMatrices(archrproj)

cellcoldata <- getCellColData(archrproj, select = c("log10(nFrags)", "TSSEnrichment"))

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
	ArchRProj = archrproj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
)

# plot TSSEnrichment ridge plot per sample
plotGroups(
	ArchRProj = archrproj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges"
)

saveArchRProject(ArchRProj = archrproj, outputDirectory = working.dir, load = TRUE)


# Filter Doublets ---------------------------------------------------------

# ?filterDoublets

archrproj <- filterDoublets(ArchRProj = archrproj, cutEnrich = 1, filterRatio = 1.5) # see notes


archrproj <- addIterativeLSI(
	ArchRProj = archrproj,
	useMatrix = "TileMatrix", 
	name = "IterativeLSI", 
	iterations = 4, 
	clusterParams = list( #See Seurat::FindClusters
		resolution = c(0.3), 
		sampleCells = 10000, 
		n.start = 10
	), 
	varFeatures = 25000, 
	dimsToUse = 1:30
)







