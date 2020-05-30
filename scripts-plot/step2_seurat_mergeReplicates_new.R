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
  name <- unlist(strsplit(countfile, "/|_"))[5]
  ########################################
  ############ Setup the Seurat Object
  #####################################
  # Load the 10X dataset
  pbmc.data <- read.table(countfile, row.names = 1, stringsAsFactors = F, header = T )
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = name)
  
  ##################################################
  ###### QC and selecting cells for further analysis
  #####################################
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^Mt")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 10)
  return(pbmc)
  
}

## Section: getting filtered counts per replicates
##################################################
file.ls <- list.files(path = 'data/DEG/', pattern = "txt", full.names = T)
count.ls <- lapply(file.ls, getCount)
sample.ls <- unlist(lapply(file.ls, function(vec){
  ls <- unlist(strsplit(vec, "/|_"))
  return(ls[5])
}))

datasetA <- merge(x = count.ls[[1]], y = c(count.ls[[2]], count.ls[[3]], count.ls[[4]], count.ls[[5]]), project = 'DataSetA',
                  merge.data = FALSE, add.cell.ids = sample.ls[1:5])
datasetB <- merge(x = count.ls[[6]], y = c(count.ls[[7]], count.ls[[8]], count.ls[[9]]), project = 'DataSetB',
                  merge.data = FALSE, add.cell.ids = sample.ls[6:9])
saveRDS(datasetA, file = 'data/DataSetA.filtered.seurat.merged.rds')
saveRDS(datasetB, file = 'data/DataSetB.filtered.seurat.merged.rds')
