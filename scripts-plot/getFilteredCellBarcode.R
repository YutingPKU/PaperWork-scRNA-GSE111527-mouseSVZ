
pbmc <- readRDS('DataSetA.seurat.pca40.umap.tsne.annotation.rds')
id <- split(rownames(pbmc@meta.data), f = pbmc@meta.data$orig.ident)
getBarcode <- function(ls){
  unlist(lapply(ls, function(vec){
    unlist(strsplit(vec, "_"))[2]
  }))
}
barcode.ls <- lapply(id, getBarcode)
