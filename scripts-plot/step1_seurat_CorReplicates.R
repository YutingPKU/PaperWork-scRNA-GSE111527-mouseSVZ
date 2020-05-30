##################################################
## Project: mouse SVZ paper work
## Script purpose: calculate the correlation of replicates based on normalized UMI
## Date: 2020-05-27
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################

setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-humanSVZ/")
library(Seurat)
library(dplyr)
#library(Matrix)


## Section: get norm-umi from DEG counts
##################################################
getNormUMI <- function(countfile){
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
  VlnPlot(pbmc, features = c('nFeautre_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
  FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 10)
  
  ########################################
  ####### Normalizing the data
  ##################################
  # After removing unwanted cells from the dataset, the next step is to normalize the data. 
  # By default, we employ a global-scaling normalization method "LogNormalize" 
  # that normalizes the gene expression measurements for each cell by the total expression, 
  # multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  nUMI <- pbmc[['RNA']]@data
  nUMI <- data.frame(id = rownames(nUMI), value =  rowMeans(nUMI))
  return(nUMI)
}

file.ls <- list.files(path = 'data/DEG/', pattern = "txt", full.names = T)
norm.umi <- lapply(file.ls, getNormUMI)

gene.ls <- lapply(norm.umi, function(mat){
  return(rownames(mat))
})
comm.gene <- Reduce(intersect, gene.ls)

mat.ls <- lapply(norm.umi, function(mat){
  return(mat[match(comm.gene, rownames(mat)),2])
})

mat <- matrix(unlist(mat.ls), ncol = 9, byrow = F)
sample.ls <- lapply(file.ls, function(vec){
  ls <- unlist(strsplit(vec, "/|_"))
  return(ls[5])
})
colnames(mat) <- unlist(sample.ls)
cor(mat, method = 'pearson')
