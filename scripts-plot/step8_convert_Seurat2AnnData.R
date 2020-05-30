##################################################
## Project: mouse SVZ paper work
## Script purpose: conversion from seurat object to annData
## Date: 2020-05-24
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-mouseSVZ")
library(Seurat)


## Section: load the data
##################################################
pbmc <- readRDS('data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.rds')
sub.pbmc <- subset(pbmc, subset= anno1 != 'OLCs')

## Section: save to loom
##################################################
lfile <- as.loom(x = sub.pbmc, filename = 'data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.noOLCs.loom')
lfile$close_all()
