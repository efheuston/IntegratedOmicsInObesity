AssignMetadata <- function(metadata.df, seurat.object){
	for(col.pos in 1:dim(metadata.df)[2]){
		col.id <- colnames(metadata.df)[col.pos]
		for(row.id in 1:dim(metadata.df)[1]){
			metadata.info <- metadata.df[row.id, col.id]
			simplified.pt.id <- stringr::str_replace_all(rownames(metadata.df)[row.id], "[^[:alnum:]]", "")
			simplified.meta.id <- stringr::str_replace_all(seurat.object@meta.data$orig.ident, "[^[:alnum:]]", "")
			
			seurat.object@meta.data[[col.id]][grepl(simplified.pt.id, simplified.meta.id, ignore.case = TRUE) == TRUE] <- metadata.info
		}
	}
	return(seurat.object)
}
