##################################################
## Project: mouse SVZ paper work
## Script purpose: using velocyto to estimate RNA velocity
## Date: 2020-05-19
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("/lustre/user/liclab/liuyt/SP/public-data/GSE111527-mouseSVZ")
library(velocyto.R)
library(pagoda2)
library(stringr)
library(BiocParallel)
library(mgsub)

## Section: load loom file
##################################################
#ldat <- read.loom.matrices("data/merge.loom")
ldat <- readRDS('data/merge.loom.rds')

## Section:pagoda2 processing
##################################################
pbmc <- readRDS('data/DataSetA.seurat.pca40.umap.tsne.annotation.rds')
# change the cell ids due to different cell ids used when running velocyto.py
id <- colnames(pbmc[["RNA"]])
id <- mgsub(id, pattern = c('an002_','an003F_','an003L_','an008_','an009_'),
            replacement = c('SRR6814853_NMRAK:',"SRR6814854_Y0E0H:","SRR6814855_EXYY0:",
                            "SRR6814856_99EJH:",  "SRR6814857_K3LEO:"))
id <- paste0(id,"x")
pbmc <- RenameCells(pbmc, new.names = id) 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000,mean.cutoff = c(0.01, 3), 
                             dispersion.cutoff  = c(1,Inf))
gn <- pbmc@assays$RNA@var.features

## Section: filtering data
##################################################
emat <- ldat$spliced
nmat <- ldat$unspliced
egn <- unlist(lapply(rownames(emat), function(vec){
  unlist(strsplit(vec,".", fixed = T))[1]
}))
ngn <- unlist(lapply(rownames(nmat), function(vec){
  unlist(strsplit(vec,".", fixed = T))[1]
}))

comm <- unique(intersect(egn, ngn))
emat <- emat[ match(comm, egn ),]
nmat <- nmat[ match(comm, ngn),]
rownames(emat) <- unlist(lapply(rownames(emat), function(vec){
  unlist(strsplit(vec,".", fixed = T))[1]
}))
rownames(nmat) <- unlist(lapply(rownames(nmat), function(vec){
  unlist(strsplit(vec,".", fixed = T))[1]
}))

loci <- match(gn, rownames(emat))
loci <- loci[-which(is.na(loci))]
emat <- emat[loci,]
nmat <- nmat[loci,]

## Section: velocity estimate
##################################################
pca.mat <- Embeddings(pbmc, reduction = "pca")
cell <- colnames(emat)
pca.mat <- pca.mat[match(cell, rownames(pca.mat)),]
cell.dist <- as.dist(1-armaCor(t(pca.mat)))

registered()
register(MulticoreParam(workers = 10, timeout = 144000000, log = F  ))
registered()
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=200, 
                                            n.cores = 10,fit.quantile=fit.quantile)


## Section: visualize on exiting embedding
##################################################
emb <- Embeddings(pbmc[["tsne"]])
emb <- emb[match(cell, rownames(emb)),]
cell.colors <- as.character(pbmc@active.ident)[match(cell, rownames(emb))]
col.ls <- mgsub::mgsub(cell.colors, pattern = c('Astrocytes' ,'NSCs' ,'OPCs' , 'TAPs','NBs', 'MFOLs1','MFOLs2' , 'MOLs' ,'COPs' ,'D2 MSNs',
                                                'D1 MSNs' , 'Ependymal','PVMs' ,'Microglia' ,'SMCs' , 'Pericytes' ,'Endothelial' ),
                       replacement = c(rgb(61,104,233, maxColorValue = 255), rgb(0,253,254, maxColorValue = 255),
                                       rgb(169,127,254,maxColorValue = 255), rgb(254,214,2,maxColorValue = 255),
                                       rgb(253,162,4, maxColorValue = 255),  rgb(216,110,210, maxColorValue = 255),
                                       rgb(156,47,239,maxColorValue = 255),  rgb(82,25,137,maxColorValue = 255),
                                       rgb(92,69,137, maxColorValue = 255), rgb(138,0,0,maxColorValue = 255),
                                       rgb(253,69,0,maxColorValue = 255), rgb(253,103,177, maxColorValue = 255),
                                       rgb(0,0,0,maxColorValue = 255), rgb(187,187,187,maxColorValue = 255),
                                       rgb(27,138,31, maxColorValue = 255), rgb(121,248,5, maxColorValue = 255),
                                       rgb(120,202,122, maxColorValue = 255))
                       )
cell.colors <- setNames(col.ls, names(pbmc@active.ident)[match(cell, rownames(emb))])


arrow.scale=3; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
pdf('results/velocyto.tsne.allcells.DataSetA.quantile002.pdf')
#show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',
#                               cell.colors=cell.colors,
#                               show.grid.flow = TRUE, min.grid.cell.mass=0.5,grid.n=100, 
#                               cex=.5,arrow.scale=3,arrow.lwd=1)
show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',
                              cell.colors=cell.colors,
                              show.grid.flow = TRUE, min.grid.cell.mass=0.5,grid.n=100, 
                              cex=.5,arrow.scale=3,arrow.lwd=1)
dev.off()

pdf('results/velocyto.pca.allcells.DataSetA.v2.pdf')
pca.velocity.plot(rvel.cd,
                  nPcs=2,
                  plot.cols=1,
                  cell.colors=cell.colors,
                  pc.multipliers=c(1,-1), ## adjust as needed to orient pcs
                  show.grid.flow = TRUE, 
                  grid.n=60 ## adjust as needed
)
dev.off()