#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plot somitic trajectory.

Author: Thomas R. Silvers (silvers@molgen.mpg.de)

Usage:
    python TLS/scripts/pseudotime/plot_somitic_trajectory.py --from-cache
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from TLS.configs.config_manager import config


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str,
                        default=config.data['TLS AnnData']['path'],
                        help='Path to the data file')
    parser.add_argument('--plot_param_key', type=str,
                        default='somitic_trajectory_umap',
                        help='Path to the plot parameters file')
    parser.add_argument('--from-cache',
                        default=True,
                        action=argparse.BooleanOptionalAction)

    args = parser.parse_args()
    
    return vars(args)


def plot_somitic_trajectory():
    args = parse_args()
    tls_adata = sc.read(args['data'], sparse=True, cache=args['from_cache'])
    
    timepoint_plot_params = config.plots[args['plot_param_key']]
    timepoint_plot_params = timepoint_plot_params.items()    
    n_timepoints = len(timepoint_plot_params)
    
    fig = plt.figure(figsize=(8, 3))
    gs = fig.add_gridspec(1, n_timepoints)

    for i, (key, params) in enumerate(timepoint_plot_params):
        params['plot_kwargs'] = dict(legend_fontsize=8)
        params['ax'] = fig.add_subplot(gs[i])
        tls_adata.tls.umap_timepoint(**params)

        params['ax'].set_title(key)
        if i > 0:
            params['ax'].set_ylabel(None)
        if i < n_timepoints - 1:
            params['ax'].get_legend().remove()
        else:
            legend = params['ax'].get_legend()
            legend.set_title('louvain cluster,\nfrom Seurat', prop={'size':9})
            for text in legend.texts:
                text.set_text(text.get_text().replace('Seurat_', ''))

    fig.suptitle('TLS differentiation: scRNA-seq, louvain clustering', y=1.1);

    return fig


if __name__ == '__main__':
    somitic_trajectory_fig = plot_somitic_trajectory()

    figures_dir = Path(config.plotting['output'])
    plt.savefig(figures_dir / 'Figure_5A_PAGApaths.pdf')