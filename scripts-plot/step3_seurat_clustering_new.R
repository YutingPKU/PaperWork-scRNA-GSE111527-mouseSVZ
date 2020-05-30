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
pbmc <- readRDS('data/DataSetA.filtered.seurat.merged.rds')

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
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:50)
ElbowPlot(pbmc, ndims = 50)

## Section: clustering the cells
##################################################
pbmc <- FindNeighbors(pbmc, dims = 1:40)
pbmc <- FindClusters(pbmc, resolution =1 )
pbmc <- RunUMAP(pbmc, dims = 1:40)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, dims = 1:40)
DimPlot(pbmc, reduction = "tsne")
saveRDS(pbmc, file = "data/DataSetA.seurat.pca40.umap.tsne.rds")


## Section: filter four clusters which do not contain all five replicates
##################################################
pbmc <- readRDS('data/DataSetA.seurat.pca40.umap.tsne.rds')
Idents(pbmc) <- "seurat_clusters"
pbmc <- subset(pbmc, idents = c(0:20,22:24))

## Section: find markers
##################################################
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

## Section: add cell type annotation
##################################################
new.cluster.ids <- c("Astrocytes", "NBs",'Microglia','MFOLs2','TAPs','Endothelial','MFOLs2','NSCs','D1 MSNs','Astrocytes',
                     'D2 MSNs','MFOLs1','OPCs','MOLs','MOLs','NBs','Pericytes','SMCs','Astrocytes','Ependymal','COPs',
                     'Astrocytes','Pericytes','PVMs')
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pdf('results/Fig1_clustering/tSNE_CellType.pdf', width = 6, height = 4.5)
DimPlot(pbmc, reduction = "tsne", label = TRUE , pt.size = 0.8,
        cols= c('Astrocytes' = rgb(61,104,233, maxColorValue = 255),'NSCs' = rgb(0,253,254, maxColorValue = 255),
                'OPCs' = rgb(169,127,254,maxColorValue = 255), 'TAPs' = rgb(254,214,2,maxColorValue = 255),
                'NBs' = rgb(253,162,4, maxColorValue = 255), 'MFOLs1' = rgb(216,110,210, maxColorValue = 255),
                'MFOLs2' = rgb(156,47,239,maxColorValue = 255), 'MOLs' = rgb(82,25,137,maxColorValue = 255),
                'COPs' = rgb(92,69,137, maxColorValue = 255),'D2 MSNs' = rgb(138,0,0,maxColorValue = 255),
                'D1 MSNs' = rgb(253,69,0,maxColorValue = 255), 'Ependymal' = rgb(253,103,177, maxColorValue = 255),
                'PVMs' = rgb(0,0,0,maxColorValue = 255),'Microglia' = rgb(187,187,187,maxColorValue = 255),
                'SMCs' = rgb(27,138,31, maxColorValue = 255), 'Pericytes' = rgb(121,248,5, maxColorValue = 255),
                'Endothelial' = rgb(120,202,122, maxColorValue = 255))) 

dev.off()

pdf('results/Fig1_clustering/UMAP_CellType.pdf')
DimPlot(pbmc, reduction = "umap", label = TRUE , pt.size = 0.5,
        cols= c('Astrocytes' = rgb(61,104,233, maxColorValue = 255),'NSCs' = rgb(0,253,254, maxColorValue = 255),
                'OPCs' = rgb(169,127,254,maxColorValue = 255), 'TAPs' = rgb(254,214,2,maxColorValue = 255),
                'NBs' = rgb(253,162,4, maxColorValue = 255), 'MFOLs1' = rgb(216,110,210, maxColorValue = 255),
                'MFOLs2' = rgb(156,47,239,maxColorValue = 255), 'MOLs' = rgb(82,25,137,maxColorValue = 255),
                'COPs' = rgb(92,69,137, maxColorValue = 255),'D2 MSNs' = rgb(138,0,0,maxColorValue = 255),
                'D1 MSNs' = rgb(253,69,0,maxColorValue = 255), 'Ependymal' = rgb(253,103,177, maxColorValue = 255),
                'PVMs' = rgb(0,0,0,maxColorValue = 255),'Microglia' = rgb(187,187,187,maxColorValue = 255),
                'SMCs' = rgb(27,138,31, maxColorValue = 255), 'Pericytes' = rgb(121,248,5, maxColorValue = 255),
                'Endothelial' = rgb(120,202,122, maxColorValue = 255))) + NoLegend()

dev.off()

