__all__ = ['ignore_warnings']

def ignore_warnings() -> None:
    import warnings

    from anndata import ImplicitModificationWarning
    from numba import NumbaDeprecationWarning

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    
    warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
    warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)