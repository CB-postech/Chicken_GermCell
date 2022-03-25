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

data_dir = "/home/dgcha/project/hs/multispecies/Human_RData"
os.chdir(data_dir)
palantir_dir = os.path.expanduser(os.getcwd())

hs_norm = pd.read_csv(os.path.join(data_dir, "hs_male_logcd_hvg.csv"), sep=",", index_col=0)
hs_norm_allgene = pd.read_csv(os.path.join(data_dir, "hs_male_logcd.csv"), sep=",", index_col=0)

pcadim = 100
dm = 50
neighbor = 30

pca_projections, _ = palantir.utils.run_pca(hs_norm, n_components=pcadim)
hs_dm = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm, knn = neighbor)
hs_ms = palantir.utils.determine_multiscale_space(hs_dm)

hs_tsne = palantir.utils.run_tsne(hs_ms, perplexity = 200)
fig, ax = palantir.plot.plot_tsne(hs_tsne)

hs_imp = palantir.utils.run_magic_imputation(hs_norm_allgene, hs_dm)

start_cell = np.argmax(hs_imp['POU5F1'])

palantir.plot.highlight_cells_on_tsne(hs_tsne, start_cell)

pr_res = palantir.core.run_palantir(hs_ms, start_cell, num_waypoints=200)
palantir.plot.plot_palantir_results(pr_res, hs_tsne)

