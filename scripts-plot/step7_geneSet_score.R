##################################################
## Project: mouse SVZ paper work
## Script purpose: calcualte gene set scores and violnplot
## Date: 2020-05-22
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-mouseSVZ")
library(Seurat)
library(openxlsx)
library(scales)

## Section: load the data
##################################################
pbmc <- readRDS('data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.rds')
gnmat <- read.xlsx('data/public/mmc3.xlsx', sheet = 1)
gn.ls <- split(gnmat$gene_symbol, f = gnmat$cluster)
gn.ls <- gn.ls[-1]


## Section: defing calculation of gene set scores
##################################################
expmat <- pbmc[["RNA"]]@data
rownames(expmat) <- toupper(rownames(expmat))
getGS <- function(gn){
  loci <- match(gn, rownames(expmat))
  loci <- loci[which(!is.na(loci))]
  df <- expmat[loci, ]
  df <- apply(df, MARGIN = 1, FUN = function(x){ rescale(x, to = c(0,1))})
  scores <- rowMeans(df)
  return(scores)
}
GS.ls <- lapply(gn.ls, getGS)

## Section: defining stacked vlnplot functions
##################################################
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(score.ls, pt.size = 0,  plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) {
  df <- data.frame(cbind(value= unlist(score.ls), type = as.character(pbmc@meta.data$anno1)))
  df <- df[ which(df$type != 'OLCs'),]
  df$value <- as.numeric(as.character(df$value))
  df$type <- factor(df$type, levels = c('qNSCs','aNSCs','TAPs','mTAPs','early NBs','late NBs'))
  
  p<- ggplot(df, aes(x=type, y=value, fill = type))+ geom_violin(trim = T, scale = 'width') + 
    xlab("") + ylab("") + ggtitle("") + theme_classic()+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)*1.4), 
          plot.margin = plot.margin ) +
    scale_fill_manual(values = c('qNSCs' = rgb(0,255,255, maxColorValue = 255), 'aNSCs' = rgb(146,252,12, maxColorValue = 255),
                                 'TAPs' = rgb(255,250,0, maxColorValue = 255), 'mTAPs' = rgb(255,218,3, maxColorValue = 255),
                                 'early NBs' = rgb(255,173,0, maxColorValue = 255), 'late NBs' = rgb(241,146,0, maxColorValue = 255),
                                 'OLCs' = rgb(176,68,234, maxColorValue = 255))) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
 # ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  ymax <- 0.4
  return(0.4)
}


## main function
StackedVlnPlot<- function(score.ls,pt.size = 0, 
                          plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) {
  
  plot_list<- purrr::map(score.ls, function(x) modify_vlnplot(x))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, size = rel(1)*1.5), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



pdf('results/Fig2_subClustering/GeneSetsScores_from_Llorens-Bobadilla_violinplot.pdf', width = 5, height = 7)
StackedVlnPlot(GS.ls)
dev.off()