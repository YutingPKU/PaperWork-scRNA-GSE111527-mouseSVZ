##################################################
## Project: mouse SVZ scRNA paper
## Script purpose: filtering cells and merge replicates
## Date: 2020-05-12
## Author: Yuting Liu
##################################################


## Section: defining envs
##################################################
rm(list = ls())
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-humanSVZ/")
library(Seurat)
library(dplyr)

## Section: defining filtering function
##################################################
getCount <- function(countfile){
  ########################################
  ############ Setup the Seurat Object
  #####################################
  # Load the 10X dataset
  pbmc.data <- read.table(countfile, row.names = 1, stringsAsFactors = F, header = T )
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "dropseq_SVZ")
  
  ##################################################
  ###### QC and selecting cells for further analysis
  #####################################
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^Mt")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 10)
  mat <- as.matrix(GetAssayData(pbmc, slot = "counts"))
  return(mat)
  
}

## Section: getting filtered counts per replicates
##################################################
file.ls <- list.files(path = 'data/DEG/', pattern = "txt", full.names = T)
count.ls <- lapply(file.ls, getCount)

gene.ls <- lapply(count.ls, function(mat){
  return(rownames(mat))
})
comm.gene <- Reduce(intersect, gene.ls)

mat.ls <- lapply(count.ls, function(mat){
  return(mat[match(comm.gene, rownames(mat)),])
})

mat.a <- cbind(mat.ls[[1]], mat.ls[[2]], mat.ls[[3]], mat.ls[[4]], mat.ls[[5]])
mat.b <- cbind(mat.ls[[6]], mat.ls[[7]], mat.ls[[8]], mat.ls[[9]])

saveRDS(mat.a, file = 'data/DataSetA.DGE.filtered.merged.rds')
saveRDS(mat.b, file = 'data/DataSetB.DGE.filtered.merged.rds')
