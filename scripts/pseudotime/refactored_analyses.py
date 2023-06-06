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
MYSTERY_CLUSTERS = ['Seurat_1', 'Seurat_3', 'Seurat_8', 'Seurat_0']
MYSTERY_IROOT = 'Seurat_4'
MYSTERY_COMPONENTS = '1,2'
MYSTERY_PATHS = [('p1', ["Seurat_4", "Seurat_7", "Seurat_2", "Seurat_5","Seurat_1-Seurat_3-Seurat_8-Seurat_0,0"])]
MYSTERY_NAME = 'p1'


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
    tls_adata_ext = tls_adata.tls.copy()
    
    
    # -- TLS 120h data processing?
    tls_adata_120h = (
        tls_adata_ext
        .query_timepoints(['Organoid_120h'])
        .paga(groups='louvain')
        .louvain(restrict_to=['louvain', MYSTERY_CLUSTERS], resolution=0, key_added='louvain2')
        .paga(groups='louvain2')
        .assign_uns(key='iroot', value_func=lambda ad: np.flatnonzero(ad.obs['louvain2'] == MYSTERY_IROOT)[0])
        .dpt() # Adds observation 'dpt_pseudotime'
        .assign_obs(key='iroot', value_func=lambda ad: ad.obs['dpt_pseudotime'])
        .scale()
    )


    # -- TLS 120h 'diffmap'?
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        components=MYSTERY_COMPONENTS,
        color='dpt_pseudotime'
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_120h.diffmap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'some_diffmap.pdf')
    
    
    # -- TLS 120h UMAP?
    # TODO:
    _FROM_CONFIG = dict(
        show=False,
        legend_loc='on data',
        color=['dpt_pseudotime']
    )
    fig, ax = plt.subplots(figsize=(4, 4))
    tls_adata_120h.umap(ax=ax, **_FROM_CONFIG)
    # plt.savefig(figures_dir / 'some_umap.pdf')
    
    
    # -- Somitic trajectory?
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
                       name=MYSTERY_NAME,
        )
    )
    # plt.savefig(figures_dir / 'Figure_5A_PAGApaths.pdf')


    # -- Neural trajectory?
    # TODO:
    _FROM_CONFIG = dict()
    (
        tls_adata_120h
        .pagapath_hmap(MYSTERY_PATHS[0][1],
                       data_assets.somitic_trajectory_genes,
                       plot_kwargs=_FROM_CONFIG,
                       height=6,
                       save_df=True,
                       name=MYSTERY_NAME,
        )
    )
    # plt.savefig(figures_dir / 'some_fig.pdf')
