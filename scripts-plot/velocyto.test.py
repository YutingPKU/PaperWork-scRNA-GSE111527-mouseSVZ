#!/bin/bash/evn python
# -*- coding:utf-8 -*-

#------------------------------------
# @author: Yuting Liu
# @email: lyt17@pku.edu.cn
# @date: 2020/5/28
#------------------------------------

import velocyto as vcy
import numpy as np
import matplotlib.pyplot as plt

# load data
vlm = vcy.VelocytoLoom('../data/merge.loom')

# filter gene with low expression
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=5)
vlm.filter_genes(by_detection_levels=True)

# select highly variable genes
vlm.score_cv_vs_mean(1000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)

# normalized by total molecular count
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

# preparation for gamma fit
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=200, balanced=True, b_sight=1600, b_maxl=800, n_jobs=10)

# gamma fit
vlm.fit_gammas()

# velocity
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


# loading tsne embedding
tsn = np.genfromtxt('../data/merge.loom.tSNE.embedding.csv',delimiter=',',skip_header=1)
vlm.ca['tsne1'] = tsn[:,0]
vlm.ca['tsne2'] = tsn[:,1]
vlm.ts = np.column_stack([vlm.ca['tsne1'], vlm.ca['tsne2']])

# setting colors
l1 = ['Astrocytes' ,'NSCs' ,'OPCs' , 'TAPs','NBs', 'MFOLs1','MFOLs2' , 'MOLs' ,'COPs' ,'D2 MSNs','D1 MSNs' , 'Ependymal','PVMs' ,'Microglia' ,'SMCs' , 'Pericytes' ,'Endothelial']
l2 =[np.array([0.23921568627451,0.407843137254902,0.913725490196078]),np.array([0,0.992156862745098,0.996078431372549]),np.array([0.662745098039216,0.498039215686275,0.996078431372549]),np.array([0.996078431372549,0.83921568627451,0.00784313725490196]),np.array([0.992156862745098,0.635294117647059,0.0156862745098039]),np.array([0.847058823529412,0.431372549019608,0.823529411764706]),np.array([0.611764705882353,0.184313725490196,0.937254901960784]),np.array([0.32156862745098,0.0980392156862745,0.537254901960784]),np.array([0.36078431372549,0.270588235294118,0.537254901960784]),np.array([0.541176470588235,0,0]),np.array([0.992156862745098,0.270588235294118,0]),np.array([0.992156862745098,0.403921568627451,0.694117647058824]),np.array([0,0,0]),np.array([0.733333333333333,0.733333333333333,0.733333333333333]),np.array([0.105882352941176,0.541176470588235,0.12156862745098]),np.array([0.474509803921569,0.972549019607843,0.0196078431372549]),np.array([0.470588235294118,0.792156862745098,0.47843137254902])]
colls = dict(zip(l1,l2))

# loading cluster id
cid = np.genfromtxt('../data/merge.loom.metadata.csv',  dtype="str", delimiter=',', skip_header=1)
vlm.set_clusters(cid, cluster_colors_dict=colls)

# filter cells and setting colors
#vlm.filter_cells(bool_array=np.in1d(vlm.cluster_labels, ['Astrocytes','NSCs','TAPs','NBs']))
#vlm.ts = np.column_stack([vlm.ca['tsne1'], vlm.ca['tsne2']])
#vlm.set_clusters(vlm.cluster_labels, cluster_colors_dict=colls)


# projection to tsne
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=200, knn_random=True)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=200)


#vlm.to_hdf5('velcyto_analysis.hdf5')
plt.figure(None,(10,10))
vlm.plot_grid_arrows()
plt.savefig("vectorfield2.pdf")
