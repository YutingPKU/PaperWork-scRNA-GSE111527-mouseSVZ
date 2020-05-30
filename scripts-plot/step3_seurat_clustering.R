##################################################
## Project: mouse SVZ scRNA paper
## Script purpose: clustering 
## Date: 2020-05-12
## Author: Yuting Liu
##################################################


## Section: defining envs
##################################################
rm(list = ls())
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-humanSVZ/")
library(Seurat)
library(dplyr)

## Section: generating seurat object
##################################################
pbmc.data <- readRDS('data/DataSetA.DGE.filtered.merged.rds')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "mSVZ", min.cells = 3, min.features = 200)

## Section: QC and filtering cells
##################################################
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^Mt")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 500)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Section: find highly variable genes
##################################################
pbmc <- FindVariableFeatures(pbmc, selection.method = "mvp", mean.cutoff = c(0.01, 3), 
                             dispersion.cutoff  = c(1,Inf))
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Section: scaling the data
##################################################
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Section: dimensional reduction
##################################################
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:25, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

## Section: clustering the cells
##################################################
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution =1 )
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "tsne")
saveRDS(pbmc, file = "data/DataSetA.seurat.pca20.umap.tsne.rds")

