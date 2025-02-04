from abc import ABC
import datetime
from pathlib import Path
import types
from typing import Dict, List, Union

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import scanpy as sc


def register_anndata_accessor(name):
    def decorator(cls):
        if hasattr(sc.AnnData, name):
            # raise AttributeError(f"Accessor {name} already exists in scanpy's AnnData.")
            pass
        setattr(sc.AnnData, name, property(lambda self: cls(self)))
        return cls
    return decorator


class BaseAnnDataMixin(ABC):

    def build_methods(self, module, exclude=None):
        for name, func in module.__dict__.items():
            if (
                callable(func)
                and not name.startswith("_")
                and (exclude is None or name not in exclude)
            ):
                setattr(self, name, types.MethodType(self._make_func(func), self))

    def _make_func(self, func):
        def wrapper(*__, **kwargs):
            func(self._obj, **kwargs)
            return self
        return wrapper


class ScanpyPreprocessMixin(BaseAnnDataMixin):
    
    _scpp_module = sc.pp
    _scpp_exclude = []
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.build_methods(module=self._scpp_module, exclude=self._scpp_exclude)


class ScanpyToolsMixin(BaseAnnDataMixin):
    
    _sctl_module = sc.tl
    _sctl_exclude = ['pca', 'tsne', 'umap', 'diffmap', 'draw_graph']
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.build_methods(module=self._sctl_module, exclude=self._sctl_exclude)


class ScanpyPlottingMixin(BaseAnnDataMixin):
    
    _scpl_module = sc.pl
    _scpl_exclude = ['paga', 'louvain', 'dpt']
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.build_methods(module=self._scpl_module, exclude=self._scpl_exclude)

        
@register_anndata_accessor("tls")
class TLSAnnDataAccessor(
    ScanpyPreprocessMixin,
    ScanpyToolsMixin,
    ScanpyPlottingMixin,
):
    def __init__(self, anndata_obj: sc.AnnData):
        self._obj = anndata_obj
        super().__init__()
        self._get_timepoints()

    @property
    def timepoints(self) -> np.ndarray:
        return self._obj.obs['donor'].cat.categories.values

    def _get_timepoints(self) -> None:
        # TODO: Temporary...
        self.timepoints_cat = self.timepoints
        self.timepoints_dt = (
            pd.Series(self.timepoints_cat)
            .str.split('_', expand=True).get(1)
            .str.split('h', expand=True).get(0)
            .apply(lambda x: datetime.timedelta(hours=float(x)))
            .values
        )
        
    def set_raw(self, raw) -> None:
        # TODO: Unclear what this is doing...
        self._obj.raw = raw
        return self

    # --- Query methods ---
    # TODO: Refactor exclude and query methods into a more general method
    
    def exclude_clusters(self, method: str, cluster_ids: List[str]) -> None:
        self._obj = self._obj[~self._obj.obs[method].isin(cluster_ids)].copy()
        return self

    def query_clusters(self, method: str, cluster_ids: List[str]) -> None:
        return self._query(cluster_ids, method)

    def query_timepoints(self, timepoints: List[str]) -> None:
        _timepoint_key = 'donor'
        return self._query(timepoints, _timepoint_key)

    def _query(self, values: Union[List, None], key: str) -> None:
        if values is None or values == [None] or len(values) == 0:
            return self

        else:
            if key not in self._obj.obs.keys():
                raise ValueError(f'Invalid key: {key}')
            
            for value in values:
                if value not in self._obj.obs[key].values:
                    raise ValueError(f'Invalid {key}: {value}')

            self._obj = self._obj[np.in1d(self._obj.obs[key], values)].copy()

            return self
    
    # --- Assign methods ---
    
    def assign_uns(self, **kwargs) -> None:
        return self._assign('uns', **kwargs)

    def assign_obs(self, **kwargs) -> None:
        return self._assign('obs', **kwargs)

    def _assign(self, ad_attr, key, value_func) -> None:
        getattr(self._obj, ad_attr)[key] = value_func(self._obj)
        return self

    # --- Methods that combine mix-in methods ---

    def umap_timepoint(
        self,
        timepoint: str = None,
        exclude: List[str] = None,
        ax: Axes = None,
        plot_kwargs: Dict = {},
    ) -> None:
        return (
            self
            .query_timepoints(timepoints=[timepoint])
            .exclude_clusters(method='louvain', cluster_ids=exclude)
            .umap(color=['louvain'], palette=sc.pl.palettes.vega_20, show=False, ax=ax, **plot_kwargs)
        )

    def pagapath_hmap(
        self,
        nodes: List[str],
        gene_set: List[str],
        plot_kwargs: Dict = {},
        height: float = 6,
        save_df: bool = False,
        name: str = None,
        save_path: str = None,
    ) -> Figure:
        """
        Plot a heatmap of gene expression along a PAGA path.
        
        Notes
        -----
        - If `save_df` is True, the resulting DataFrame will be saved to `self.uns[name]`.
        - Because of how `scanpy.pl.paga_path` is implemented, makes most sense to
            create a figure and axes object inside of this method.
        """
        # -- Configure data --
        # Withouth assigning `distance`, `paga_path` throws a type error (??)
        self.assign_obs(key='distance', value_func=lambda ad: ad.obs['dpt_pseudotime'])
        
        # -- Configure plot --
        default_plot_kwargs = dict(
            show_node_names=False,
            ytick_fontsize=12,
            left_margin=0.5,
            show_colorbar=True,
            normalize_to_zero_one=True,
        )
        default_plot_kwargs.update(plot_kwargs)
        fig, ax = plt.subplots()
        paga_path_kwargs = dict(
            nodes=nodes,
            keys=gene_set,
            annotations=['distance'],
            return_data=save_df,
            ax=ax,
            show=False
        )

        # -- Plot --
        if save_df is False:
            self.paga_path(**paga_path_kwargs, **default_plot_kwargs)

        # Registered extension method doens't work to return a DataFrame
        else:
            __, df = sc.pl.paga_path(self._obj, **paga_path_kwargs, **default_plot_kwargs)

            # Save results
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

        # Aesthetics
        ax.set_frame_on(False)
        fig.set_size_inches(height, len(gene_set)/(height*(2/3)))

        return fig

    # ~*~

    @classmethod
    def from_file(cls, path) -> None:
        return cls(sc.read(path, sparse=True))

    def copy(self) -> None:
        return TLSAnnDataAccessor(self._obj.copy())

    def __repr__(self):
        return (f'TLSAccessor object with '
                f'{self._obj.shape[0]} cells and '
                f'{self._obj.shape[1]} variables')
    
    # ~~~~ Methods that are not implemented ~~~~
    
    def query():
        raise NotImplementedError

    def assign():
        raise NotImplementedError