#!/bin/bash/evn python
# -*- coding:utf-8 -*-

#------------------------------------
# @author: Yuting Liu
# @email: lyt17@pku.edu.cn
# @date: 2020/5/24
#------------------------------------

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.pyplot import plot,savefig
import scanpy as sc



# fig setting
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './testnoOLCs.h5ad'
sc.settings.set_figure_params(dpi=300)

# loading the dataset
adata = sc.read_loom('data/DataSetA.neurogenicSubset.seurat.pca10.umap.tsne.annotation.noOLCs.loom')

# pre-process, diffmap and trajectory inference
adata.uns['iroot'] = np.flatnonzero(adata.obs['ClusterName']  == 'qNSCs')[10]
sc.pp.neighbors(adata, n_neighbors=20, use_rep='X', method='gauss')
sc.tl.diffmap(adata,n_comps=10)
sc.tl.dpt(adata, n_branchings=0, n_dcs=10)
#adata = read_h5ad('./test.h5ad')

adata.obsm['tsne'] = adata.obsm['tsne_cell_embeddings']
fig, ax = pl.subplots(figsize=(12, 12))
axs = sc.pl.tsne(adata, color=['dpt_pseudotime', 'ClusterName'], show=False)
pl.show()
pl.savefig('./testnoOLCs.pdf')
adata.write(results_file)
adata.obs.to_csv( './testnoOLCs.csv')

