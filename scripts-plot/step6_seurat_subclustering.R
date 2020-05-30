##################################################
## Project: mouse SVZ paper work
## Script purpose: subclustering of neurogenic lineage
## Date: 2020-05-22
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-mouseSVZ")
library(Seurat)

## Section: load the data
##################################################
pbmc <- readRDS('data/DataSetA.seurat.pca40.umap.tsne.annotation.rds')
pbmc <- subset(pbmc, idents = c('NSCs','TAPs','NBs'))


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
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8 )
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
#saveRDS(pbmc, file = "data/DataSetA.neuronal.lineage.subcluster.seurat.pca40.umap.tsne.rds")

## Section: find markers
##################################################
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

## Section: add cell type annotation
##################################################
new.cluster.ids <- c('early NBs','qNSCs','mTAPs','late NBs','TAPs','late NBs','aNSCs','OLCs')
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

pdf('results/Fig2_subClustering/neurogenicLineage_tSNE_CellType.pdf', width = 6, height = 4.5)
DimPlot(pbmc, reduction = "tsne", label = TRUE , pt.size = 0.8,
        cols= c('qNSCs' = rgb(0,255,255, maxColorValue = 255), 'aNSCs' = rgb(146,252,12, maxColorValue = 255),
                'TAPs' = rgb(255,250,0, maxColorValue = 255), 'mTAPs' = rgb(255,218,3, maxColorValue = 255),
                'early NBs' = rgb(255,173,0, maxColorValue = 255), 'late NBs' = rgb(241,146,0, maxColorValue = 255),
                'OLCs' = rgb(176,68,234, maxColorValue = 255))) 

dev.off()

pdf('results/Fig2_subClustering/neurogenicLineage_UMAP_CellType.pdf')
DimPlot(pbmc, reduction = "umap", label = TRUE , pt.size = 0.5,
        cols= c('qNSCs' = rgb(0,255,255, maxColorValue = 255), 'aNSCs' = rgb(146,252,12, maxColorValue = 255),
                'TAPs' = rgb(255,250,0, maxColorValue = 255), 'mTAPs' = rgb(255,218,3, maxColorValue = 255),
                'early NBs' = rgb(255,173,0, maxColorValue = 255), 'late NBs' = rgb(241,146,0, maxColorValue = 255),
                'OLCs' = rgb(176,68,234, maxColorValue = 255))) + NoLegend()

dev.off()

type <- data.frame(pbmc@active.ident)
rownames(type) <- colnames(pbmc)
pbmc[['anno1']] <- type
saveRDS(pbmc, file = 'data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.rds')
