#python version: 3.6.10

import os
import sys
import glob
import palantir #ver 0.2.6
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
import random
import pickle

%matplotlib inline

data_dir = "/home/dgcha/project/chicken/palantir/"
os.chdir(data_dir)
palantir_dir = os.path.expanduser(os.getcwd())

norm_df_allgene = pd.read_csv(os.path.join(data_dir, "logcd_female_allgene.csv"), sep=",", index_col=0)

pca_max = pd.read_csv(os.path.join(data_dir, "Female_200PCAembeddings_750hvgs.csv"), sep = ",", index_col = 0)

pcadim = 20
dm = 10
neighbor = 30

pca_projections = pca_max.iloc[:, 0:pcadim]
pca_projections.columns = list(range(pcadim))

dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm, knn = neighbor)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

tsne = palantir.utils.run_tsne(ms_data)
fig, ax = palantir.plot.plot_tsne(tsne)

imp_df = palantir.utils.run_magic_imputation(norm_df_allgene, dm_res)

start_cell = np.argmax(imp_df['CXCR4'])

palantir.plot.highlight_cells_on_tsne(tsne, start_cell)

pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=750)
palantir.plot.plot_palantir_results(pr_res, tsne)

