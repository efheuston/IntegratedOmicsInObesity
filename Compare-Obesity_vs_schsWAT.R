# Define global variables
working.dir <- "~/Desktop/Obesity/scRNA/"
fc.thresh <- 1.0

# Load libraries
library(dplyr)
library(ggVennDiagram)


# Set working dir
setwd(working.dir)
# Load data

#load Obesity data: positive FindAllMarkers
panc <- readRDS("Obesity_scRNA-Anchored-NW-OB-posmarkers-90pctvar.rds")

#require p_val_adj < 0.05 & avg_log2FC > fc.thresh
panc.thresh <- panc %>%
  filter(p_val_adj <= 0.05) %>%
  filter(avg_log2FC >= fc.thresh)
View(panc.thresh)

#load schsWAT data: positive FindAllMarkers (by celltype)
celltype.wat <- readRDS("schsWAT-posmarkers-95pctvar.rds")

#require p_val_adj < 0.05 & avg_log2FC > fc.thresh
celltype.thresh <- celltype.wat %>%
  filter(p_val_adj <= 0.05) %>%
  filter(avg_log2FC >= fc.thresh)
#create a new column with more general categories
celltype.thresh$tissue <- celltype.thresh$cluster
celltype.thresh$tissue <- gsub("adipocyte|ASPC|pericyte", "adipose", celltype.thresh$tissue)
celltype.thresh$tissue <- gsub("macrophage|monocyte|dendritic_cell|mast_cell|neutrophil|b_cell|nk_cell|t_cell", "immune", celltype.thresh$tissue)
celltype.thresh$tissue <- as.factor(celltype.thresh$tissue)

cluster.v.cluster <- function(x.levels, x.comparison, y.levels, y.comparison){
  venn.list <- c()
  x.name <- "x.comparison"
  y.name <- "y.comparison"
  for(x.level in levels(x.levels)){
    x.compare <- toupper(unique(x.comparison[x.levels == x.level]))
    x.level.name <- paste(x.name, as.character(x.level), sep = ".level_")
    for(y.level in levels(y.levels)){
      y.compare <- toupper(unique(y.comparison[y.levels == y.level]))
      y.level.name <- paste(y.name, as.character(y.level), sep = ".level_")
      venn.compare <- list(x.level.name = x.compare, 
                                y.level.name= y.compare)
      if(length(intersect(x.compare, y.compare))>0){
        names(venn.compare) <- c(x.level.name, y.level.name)
        venn.list <- append(venn.list, list(venn.compare))
      }
    }
  }
  return(venn.list)
}

venn.cluster<- cluster.v.cluster(x.levels = panc.thresh$cluster, 
                                    x.comparison = panc.thresh$gene, 
                                    y.levels = celltype.thresh$cluster, 
                                    y.comparison = celltype.thresh$gene)
lapply(venn.cluster, ggVennDiagram)


venn.tissue <- cluster.v.cluster(x.levels = panc.thresh$cluster, 
                                 x.comparison = panc.thresh$gene, 
                                 y.levels = celltype.thresh$tissue, 
                                 y.comparison = celltype.thresh$gene)
lapply(venn.tissue, ggVennDiagram)
