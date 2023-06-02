from abc import ABC
import datetime
import types
from typing import Dict, List, Union

from matplotlib.axes import Axes
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