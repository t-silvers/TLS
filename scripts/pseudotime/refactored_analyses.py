# TODO: It's unclear from original which analyses are actually performed
#      in the paper. This script is a placeholder for now demonstrating how
#      to load the data and perform the analyses using the refactored code.
#
#       I try to *roughly* separate the analyses into sections that are
#       demarcated by `# --`.
#
#       Once I can figure out what's happening in the original code, I'll
#       revise and clean up this script. I'll also update the configs.
#       - TRS 2023-06-06
#
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import scanpy as sc

from TLS.configs.config_manager import config
from TLS.data.assets.assets_manager import data_assets


# TODO:
CLUSTER_REF = dict(
    somite_1='Seurat_0',
    somite_2='Seurat_1',
    # Seurat_2
    somite_3='Seurat_3',
    meso='Seurat_4',
    # Seurat_5
    neuro2='Seurat_6',
    neuro_progenitors='Seurat_7',
    somite_4='Seurat_8',
    nmps='Seurat_9',
    # Seurat_10
    # Seurat_11
    # Seurat_12
    # Seurat_13
)
HOX_CLUSTERS = dict(
    nmps='Seurat_9',
    meso='Seurat_4',
    neuro_progenitors='Seurat_7',
)
SOMITE_CLUSTERS = dict(
    somite_1='Seurat_0',
    somite_2='Seurat_1',
    somite_3='Seurat_3',
    somite_4='Seurat_8'
)
ALL_TPS = ['Organoid_96h', 'Organoid_108h', 'Organoid_120h']
MYSTERY_CLUSTERS = ['Seurat_1', 'Seurat_3', 'Seurat_8', 'Seurat_0']
MYSTERY_IROOT = 'Seurat_4'
MYSTERY_COMPONENTS = '1,2'
MYSTERY_PATHS = [('p1', ["Seurat_4", "Seurat_7", "Seurat_2", "Seurat_5","Seurat_1-Seurat_3-Seurat_8-Seurat_0,0"])]
MYSTERY_NAME = 'p1'
REMOVE_CLUSTERS_96 = ['Seurat_0-Seurat_1-Seurat_8-Seurat_3,2','Seurat_0-Seurat_1-Seurat_8-Seurat_3,3']


# TODO: Factor out these common palettes to `src` and `config`
reds_cmap = mpl.colormaps['Reds'](np.linspace(0, 1, 128))
greys_r_cmap = mpl.colormaps['Greys_r'](np.linspace(0.7, 0.8, 20))
red_grey_r = np.vstack([greys_r_cmap, reds_cmap])
red_grey_r_cmap = LinearSegmentedColormap.from_list('red_grey_r_cmap', red_grey_r)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str,
                        default=config.data['TLS AnnData']['path'],
                        help='Path to the data file')
    parser.add_argument('--from-cache',
                        default=True,
                        action=argparse.BooleanOptionalAction)

    args = parser.parse_args()
    
    return vars(args)


if __name__ == '__main__':
    
    # --
    args = parse_args()
    figures_dir = Path(config.plotting['output'])
    tls_adata = sc.read(args['data'], sparse=True, cache=args['from_cache'])    
    
    
    # -- TLS 120h data processing
    tls_adata_120h = (
        tls_adata
        .tls
        .copy()
        .query_timepoints(['Organoid_120h'])
        .paga(groups='louvain')
        .louvain(restrict_to=['louvain', MYSTERY_CLUSTERS], resolution=0, key_added='louvain2')
        .paga(groups='louvain2')
        .assign_uns(key='iroot', value_func=lambda ad: np.flatnonzero(ad.obs['louvain2'] == MYSTERY_IROOT)[0])
        .dpt() # Adds observation 'dpt_pseudotime'
        .assign_obs(key='iroot', value_func=lambda ad: ad.obs['dpt_pseudotime'])
        .scale()
    )


    # -- TLS 120h 'diffmap'
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        components=MYSTERY_COMPONENTS,
        color='dpt_pseudotime'
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_120h.diffmap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'some_diffmap.pdf')
    
    
    # -- TLS 120h UMAP
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        legend_loc='on data',
        color=['dpt_pseudotime']
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_120h.umap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'some_umap.pdf')
    
    
    # -- Somitic trajectory
    # TODO:
    _FROM_CONFIG = dict(
        color_map=red_grey_r_cmap
    )
    (
        tls_adata_120h
        .pagapath_hmap(MYSTERY_PATHS[0][1],
                       data_assets.somitic_trajectory_genes,
                       plot_kwargs=_FROM_CONFIG,
                       height=6,
                       save_df=True,
                       name=MYSTERY_NAME, # TODO: shouldn't this be 'somitic'?
                       )
    )
    # plt.savefig(figures_dir / 'Figure_5A_PAGApaths.pdf')


    # -- Neural trajectory
    # TODO:
    _FROM_CONFIG = dict(
        color_map=red_grey_r_cmap
    )
    (
        tls_adata_120h
        .pagapath_hmap(MYSTERY_PATHS[0][1],
                       data_assets.neural_trajectory_genes,
                       plot_kwargs=_FROM_CONFIG,
                       height=6,
                       save_df=True,
                       name=MYSTERY_NAME, # TODO: shouldn't this be 'neural'?
                       )
    )
    # plt.savefig(figures_dir / 'Figure_5A_PAGApaths.pdf')


    # -- Plot Hox Genes pseudotime for NMPs, meso and neuro progenitors
    tls_adata_all_hox = (
        tls_adata
        .tls
        .copy()
        .query_clusters(method='louvain', clusters=HOX_CLUSTERS.values())
        .paga(groups='donor')
        
        # TODO: Unclear what's happening here:
        .assign_uns(key='iroot', value_func=lambda ad: np.flatnonzero(ad.obs['donor'] == 'Organoid_96h')[0])
        
        .dpt()
    )
    

    # -- TLS all Hox 'diffmap'
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        components=MYSTERY_COMPONENTS,
        color='dpt_pseudotime'
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_all_hox.diffmap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'hox_diffmap.pdf')
    

    # -- TLS all Hox umap
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        legend_loc='on data',
        color=['dpt_pseudotime']
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_all_hox.umap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'hox_umap.pdf')


    # -- TLS all Hox PAGA path hmap
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        legend_loc='on data',
        color=['dpt_pseudotime']
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    (
        tls_adata_all_hox
        .umap(**_FROM_CONFIG)
    )
    # plt.savefig(figures_dir / 'hox_umap.pdf')


    # -- example HoxC
    # TODO:
    _FROM_CONFIG = dict(
        color_map='magma',
        color_maps_annotations={'distance': 'viridis'},
        title='Time path'
    )
    (
        tls_adata_all_hox
        .scale()
        .pagapath_hmap(ALL_TPS,
                       data_assets.hoxc_genes,
                       plot_kwargs=_FROM_CONFIG,
                       height=6,
                       save_df=True,
                       name='Time path', # TODO: shouldn't this be 'HoxC'?
                       )
    )
    # plt.savefig(figures_dir / 'hoxc_paga_paths.pdf')
    
    
    # -- Subclustering Neural Tube Cluster
    tls_adata_neuro = (
        tls_adata
        .tls
        .copy()
        .query_clusters(method='louvain', clusters=[CLUSTER_REF.neuro2])
    )
    
    # Clusters UMAP
    tls_adata_neuro.umap(color=['louvain'])
    
    # More data processing
    tls_adata_120h_neuro = (
        tls_adata_neuro
        .query_timepoints(['Organoid_120h'])
        .louvain(restrict_to=('louvain', [CLUSTER_REF.neuro2]),
                 resolution=0.65, key_added='louvain_sub')
    )
    
    # Subclusters UMAP
    tls_adata_120h_neuro.umap(color=['louvain_sub'])
    
    # Cell cycle matrix plot
    tls_adata_120h_neuro.matrixplot(var_names=data_assets.cell_cycle_genes,
                                    groupby='louvain_sub',
                                    standard_scale='var',
                                    cmap='Reds')

    # -- Subclustering Somite Cluster 120h
    tls_adata_somite = (
        tls_adata
        .tls
        .copy()
        .query_clusters(method='louvain', clusters=SOMITE_CLUSTERS.values())
        
        # TODO: Unclear if calling this before rendering initial UMAP is correct:
        .louvain(restrict_to=('louvain', SOMITE_CLUSTERS.values()),
                 resolution=0.3,
                 key_added='louvain_sub')   
    )
    
    # Clusters UMAP
    tls_adata_somite.umap(color=['louvain'])
    
    # Subclusters UMAP
    tls_adata_somite.umap(color=['louvain_sub'])
    
    # ----- Restrict to 120h
    tls_adata_120h_somite = (
        tls_adata_somite
        .query_timepoints(['Organoid_120h'])
    )
    
    # Rank genes plots
    tls_adata_120h_somite.rank_genes_groups(groupby='louvain_sub',
                                            key_added='test',
                                            rankby_abs=False,
                                            use_raw=True)

    tls_adata_120h_somite.rank_genes_groups_matrixplot(key='test',
                                                       dendrogram=False,
                                                       standard_scale='var',
                                                       cmap='Reds')

    # ----- Restrict to 96h
    tls_adata_96h_somite = (
        tls_adata_somite
        .query_timepoints(['Organoid_96h'])
        .louvain(restrict_to=('louvain', SOMITE_CLUSTERS.values()),
                 resolution=0.3, key_added='louvain_sub')
    )
    
    # Subclusters UMAP
    tls_adata_96h_somite.umap(color='louvain_sub')
    
    # Subclusters UMAP w some clusters removed
    (
        tls_adata_96h_somite
        .exclude_clusters(method='louvain_sub', clusters=REMOVE_CLUSTERS_96)
        .umap(color=['louvain_sub'])
    )