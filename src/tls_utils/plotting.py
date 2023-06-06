from pathlib import Path

from matplotlib.axes import Axes
import scanpy as sc


def plot_pagapath(
    adata: sc.AnnData,
    nodes: list,
    genes: list,
    ax: Axes,
    *args,
    annotations: list = ['distance'], # TODO: Include rename of `dpt_pseudotime` to `distance` in docs
    name: str = None,
    save_results: bool = False,
    save_path: str = None,
    **kwargs,
):
    # Set custom defaults
    if 'show_node_names' not in kwargs:
        kwargs['show_node_names'] = False

    if 'ytick_fontsize' not in kwargs:
        kwargs['ytick_fontsize'] = 12
    
    if 'left_margin' not in kwargs:
        kwargs['left_margin'] = 0.5
    
    if 'show_colorbar' not in kwargs:
        kwargs['show_colorbar'] = True

    if 'normalize_to_zero_one' not in kwargs:
        kwargs['normalize_to_zero_one'] = True

    # Plotting
    __, df = sc.pl.paga_path(
        adata, nodes, genes,
        *args,
        annotations=annotations,
        return_data=True,
        ax=ax,
        show=False,
        **kwargs
    )
    
    # Aesthetics
    ax.set_frame_on(False)
    
    # Save results
    if save_results:
        if save_path is None:
            if name is None:
                name = '_'.join(nodes)
            save_path = f'./reports/results/pagapath_{name}.csv'
            
        if not isinstance(save_path, Path):
            save_path = Path(save_path)
            
        if save_path.exists():
            raise FileExistsError(f'File {save_path} already exists.')

        if not save_path.parent.exists():
            save_path.parent.mkdir(parents=True)

        df.to_csv(save_path)
