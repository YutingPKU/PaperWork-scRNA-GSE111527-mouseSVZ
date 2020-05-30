##################################################
## Project: mouse SVZ scRNA paper
## Script purpose: makers expression by vlnplot, feature plot and heatmap
## Date: 2020-05-15
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
rm(list = ls())
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-humanSVZ/")
library(Seurat)
library(dplyr)
library("viridis") 


## Section: defining stacked vlnplot functions
##################################################
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_blank(), 
          plot.margin = plot.margin ) +
    scale_fill_manual(values = c('Astrocytes' = rgb(61,104,233, maxColorValue = 255),'NSCs' = rgb(0,253,254, maxColorValue = 255),
                                    'OPCs' = rgb(169,127,254,maxColorValue = 255), 'TAPs' = rgb(254,214,2,maxColorValue = 255),
                                    'NBs' = rgb(253,162,4, maxColorValue = 255), 'MFOLs1' = rgb(216,110,210, maxColorValue = 255),
                                    'MFOLs2' = rgb(156,47,239,maxColorValue = 255), 'MOLs' = rgb(82,25,137,maxColorValue = 255),
                                    'COPs' = rgb(92,69,137, maxColorValue = 255),'D2 MSNs' = rgb(138,0,0,maxColorValue = 255),
                                    'D1 MSNs' = rgb(253,69,0,maxColorValue = 255), 'Ependymal' = rgb(253,103,177, maxColorValue = 255),
                                    'PVMs' = rgb(0,0,0,maxColorValue = 255),'Microglia' = rgb(187,187,187,maxColorValue = 255),
                                    'SMCs' = rgb(27,138,31, maxColorValue = 255), 'Pericytes' = rgb(121,248,5, maxColorValue = 255),
                                    'Endothelial' = rgb(120,202,122, maxColorValue = 255))) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

## Section: vlnplot
##################################################
pbmc <- readRDS('data/DataSetA.seurat.pca40.umap.tsne.annotation.rds')
levels(x = pbmc) <- c('Astrocytes','NSCs','TAPs','NBs','OPCs','COPs','MFOLs1','MFOLs2','MOLs','Ependymal','Endothelial',
                      'Pericytes','SMCs','Microglia','PVMs','D1 MSNs','D2 MSNs')
gn.ls <- c('Slc1a3','Thbs4','Ascl1','Dcx', 'Pdgfra','Fyn','Mog','Klk6','Foxj1','Flt1','Vtn','Acta2','P2ry12','Mrc1','Tac1','Penk')

pdf('results/Fig1_clustering/MakerGene_stacked_Vlnplot.pdf', width = 6, height = 6.5)
StackedVlnPlot(pbmc, features = rev(gn.ls))
dev.off()

## Section: heatmap
##################################################
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

g <- DoHeatmap(pbmc, features = top10$gene, angle = 90, group.bar = T, raster = F) + scale_fill_viridis(option = "D") 
pdf('results/Fig1_clustering/MakerGene_heatmap.pdf')
g
dev.off()

## Section: featureplot
##################################################
pdf('results/Fig1_clustering/MakerGene_FeaturePlot.pdf', width =12, height = 9 )
FeaturePlot(pbmc, features = c('Slc1a3','Aqp4','S100b','Thbs4','Ascl1','Mki67','Dcx'), reduction = 'tsne',
            cols = c('yellow','red'))
dev.off()
