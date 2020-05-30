##################################################
## Project: mouse SVZ paper work
## Script purpose: plot the pseudotime on tsne, pseudotime is calculated by scanpy
## Date: 2020-05-26
## Author: Yuting Liu
##################################################

## Section: set env
##################################################
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-mouseSVZ")
library(Seurat)

## Section: load the data
##################################################
pbmc <- readRDS('data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.rds')
ann <- read.csv('pycharm-scanpy/test.csv')

## Section: plot pseudotime on tsne
##################################################
type <- data.frame(ann$dpt_pseudotime)
rownames(type) <- colnames(pbmc)
pbmc[['pseudotime']] <- type
#Idents(pbmc) <- pbmc[['pseudotime']]

FeaturePlot(pbmc, features = 'pseudotime', reduction = 'tsne')

## Section: plot gene expression, cells ordered by pseudotime, expression are smoothed with loess regression
##################################################
df <- data.frame(cbind(value = matrix(unlist(pbmc[['RNA']][grep('Slc1a3',rownames(pbmc)),]))[,1], time = unlist(ann$dpt_pseudotime)))
df$value <- as.numeric(as.character(df$value))
df$time <- as.numeric(as.character(df$time))

ggplot(df, aes(time, value)) +
  geom_point() +
  geom_smooth()
