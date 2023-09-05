# IPA tables comparing scWAT and Obesity_scRNA-Anchored-NW-OB-allmarkers-90pctvar (grouped by cell-type defined via scRNA+snATAC data)
library(dplyr)
library(Hmisc)
library(corrplot)
library(stringr)
				
setwd("~/Desktop/Obesity/scRNA/")

color.palatte <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

#Rename fxn

rename.cols <- function(x){
	sample.names <- colnames(x)
	sample.names <- unlist(lapply(sample.names, function(x) str_split(x, "\\.\\.\\.", simplify = TRUE)[1]))
	sample.names <- sapply(sample.names, function(x) gsub(pattern = "\\.allmarkers\\.95pctvar", replacement = "", x = x), simplify = TRUE)
	sample.names <- sapply(sample.names, function(x) gsub(pattern = "_scRNA\\.Anchored\\.NW\\.OB\\.allmarkers\\.90pctvar", replacement = "", x = x), simplify = TRUE)
	colnames(x) <-sample.names
	return(x)
}

ipa.table <- read.table("WAT_immune-Vs-All-ToxFxn.txt", skip = 2, header = TRUE, sep = "\t", na.strings = "N/A", row.names = 1)
ipa.table <- ipa.table[,!names(ipa.table) %in% c("X26DEM_IPA.FDR0.05.adj_age.sex...2023.07.21.01.49.PM")]
ipa.table <- rename.cols(ipa.table)

ipa.table <- ipa.table %>%
	filter(if_any(everything(), ~!is.na(.)))


ipa.rcorr <- rcorr(as.matrix(ipa.table), type = "pearson")
corrplot(corr=ipa.rcorr$r,p.mat = ipa.rcorr$P, type = "full", insig = "blank", sig.level =.1, pch.cex = .9, col = color.palatte(200))


ipa.table <- read.table("WAT_immune-Vs-All-ToxFxn.txt", skip = 2, header = TRUE, sep = "\t", na.strings = "N/A", row.names = 1)
ipa.table <- ipa.table[,!names(ipa.table) %in% c("X26DEM_IPA.FDR0.05.adj_age.sex...2023.07.21.01.49.PM")]

ipa.table <- ipa.table %>%
	filter(if_any(everything(), ~!is.na(.)))


ipa.rcorr <- rcorr(as.matrix(ipa.table), type = "pearson")
corrplot(corr=ipa.rcorr$r,p.mat = ipa.rcorr$P, type = "full", insig = "blank", sig.level =.1, pch.cex = .9, col = color.palatte(200))
