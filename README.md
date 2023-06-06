# Mouse embryonic stem cells self-organize into trunk-like structures with neural tube and somites

## Note on code repository

The original code repository is preserved on the `master` branch or available in the [original repository](https://github.com/HeleneKretzmer/TLS).

In the `legacy-refactor` branch, the following major changes have been made:

- :white_check_mark: Directory restructuring
- :white_check_mark: Develop `src` code and separate it from `scripts` code
- :white_check_mark: Refactored `scripts/Pseudotime/Pseudotime_and_Subclustering.py`
- :white_check_mark: Data and plot specs can be set in config files
- :white_check_mark: A dashboard at `reports/app`
- :white_square_button: Dashboard Improvements
- :white_square_button: Legacy code refactoring
- :white_square_button: Cache data intermediates, make public
- :white_square_button: Write tests


A more detailed roadmap can be found in `reports/app/pages/4_🎙️_Feedback.py`


### Registering custom `AnnData` accessors

A key strength and weakness of `Anndata` objects, and in turn the entire `scanpy` ecosystem, is that analyses modify data inplace by default without a transparent 'paper trail' of analyses. While memory-efficient and user friendly,
- resulting code is often difficult to follow
- the attribute space of each object becomes cluttered
- uncertainty which analyses or data generated object

Libraries can use the decorator `anndata_extensions.register_anndata_accessor()` to add additional “namespaces” to `AnnData` objects. All of these follow a similar convention: you decorate a class, providing the name of attribute to add. The class’s __init__ method gets the object being decorated. For example:

```python
@anndata_extensions.register_anndata_accessor('devo')
class DevoAccessor:

    _tp_key = 'donor'

    def __init__(self, anndata_obj: sc.AnnData):
        self._obj = anndata_obj
        self._get_timepoints()

    @property
    def timepoints(self) -> np.ndarray:
        return self._obj.obs[self._tp_key].cat.categories.values

    def _get_timepoints(self) -> None:
        self.timepoints_dt = (
            pd.Series(self.timepoints)
            .str.split('_', expand=True).get(1)
            .str.split('h', expand=True).get(0)
            .apply(lambda x: datetime.timedelta(hours=float(x)))
            .values
        )

    def query_timepoints(self, timepoints: List[str]) -> None:
        return self._query(timepoints, self._tp_key)

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
```
Now users can access your methods using the `devo` namespace:
```python
adata = sc.read(path_to_h5ad, sparse=True, cache=True)
adata.devo.query_timepoints(['E9.5'])
```

An accessor that extends `AnnData` objects with methods from `scanpy` modules is part of the namespace already (`tls`) when this repository is used as a python package (see below).

To demonstrate its usage, we can replicate a good chunk of [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html). Note the number of parameters and desctructive data operations that are performed in the original code. While `adata.raw = adata` can 'freeze' analysis for use later, creating a clean slate, this approach does not scale and results in similarly opaque `AnnData` objects.

Original code (abridged):
```python
adata.var_names_make_unique()
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')
sc.pl.pca_variance_ratio(adata, log=True)
```

Revised code (abridged):
```python
data_conf = dict(
    n_top=20,
    min_genes=200,
    min_cells=3,
    mt_varname='MT-',
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True,
    min_n_genes_by_counts=2500,
    min_pct_counts_mt=5,
    target_sum=1e4,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    regress=['total_counts', 'pct_counts_mt'],
    max_value=10
)

# Prepare data
# -- First plotting object
pbmc3k_part_one_adata = (
    adata
    .tls
    .var_names_make_unique()
    .highest_expr_genes(n_top=data_conf['n_top'])
    .filter_cells(min_genes=data_conf['min_genes'])
    .filter_genes(min_cells=data_conf['min_cells'])
    .query_var('mt', lambda x: x.var_names.str.startswith(data_conf['mt_varname']))
    .calculate_qc_metrics(**data_conf)
)

# -- Second plotting object
pbmc3k_part_two_adata = (
    pbmc3k_part_one_adata
    # .copy() # Optional
    .query_obs('n_genes_by_counts', lambda x: x < data_conf['min_n_genes_by_counts'])
    .query_obs('pct_counts_mt', lambda x: x < data_conf['min_pct_counts_mt'])
    .normalize_total(target_sum=data_conf['target_sum'])
    .log1p()
    .highly_variable_genes(**data_conf)
    .raw
    .filter_highly_variable()
    .regress_out(data_conf['regress'])
    .scale(max_value=data_conf['max_value'])
)


# Make plots
plot_conf = dict() # Can specify kwargs and aesthetics here
fig, axs = plt.subplots(2, 3)

# -- Part one
# ----- violin plot
pbmc3k_part_one_adata.violin(['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                             jitter=0.4, multi_panel=True, ax=axs[0, 0])
# ----- scatter plot (pct_counts_mt)
pbmc3k_part_one_adata.scatter(x='total_counts', y='n_genes_by_counts', ax=axs[0, 1])
# ----- scatter plot (n_genes_by_counts)
pbmc3k_part_one_adata.scatter(x='total_counts', y='n_genes_by_counts', ax=axs[0, 2])

# -- Part two
# ----- pca plot one
pbmc3k_part_two_adata.pca(svd_solver='arpack', ax=axs[1, 0])
# ----- pca plot two
pbmc3k_part_two_adata.pca(color='CST3', ax=axs[1, 1])
# ----- pca plot var ratio
pbmc3k_part_two_adata.pca_variance_ratio(log=True, ax=axs[1, 2])

```



## Abstract

> Post-implantation embryogenesis is a highly dynamic process comprising multiple lineage decisions and morphogenetic changes inaccessible to deep analysis in vivo. Here, we show that pluripotent mouse embryonic stem cells (mESCs) form aggregates that upon embedding in an extra-cellular matrix compound induce the formation of highly organized "Trunk-Like-Structures (TLS)" comprising the neural tube and somites. Comparative single-cell RNA-seq analysis confirms that this process is highly analogous to mouse development and follows the same step-wise gene-regulatory program. Tbx6 knockout TLS develop additional neural tubes mirroring the embryonic mutant phenotype, and chemical modulation can induce excess somite formation. Trunk-like-structures thus reveal an unprecedented level of self-organization and provide a powerful platform for investigating post-implantation embryogenesis in a dish.

## Code availability

This repository contains the collection of R and python snippets and scripts used to perform the transcriptomics analysis presented in Veenvliet and Bolondi et al. The main processing and data integration steps are show, based on which general visualization and plots of e.g. single genes were done.
A more detailed description and references for used tools can be found in the [corresponding publication](https://www.biorxiv.org/content/10.1101/2020.03.04.974949v1).

### Setup virtual environment


#### Using python's `venv`

1. Ensure you have both Python and R installed on your system, as the setup involves creating a Python virtual environment and installing R packages.

2. Create a new virtual environment in a directory named `tls-env`:
    
    ```bash
    python src/venv_setup.py tls-env
    ```

3. Activate the virtual environment. On Unix or MacOS, run `source tls-env/bin/activate`.

    On Windows, run `.\tls-env\Scripts\activate`.

4. Install the Python project requirements with:

    ```bash
    pip install -r requirements.txt
    ```

5. The R packages listed in the `requirements.R` file will be automatically installed during the creation of the virtual environment. Make sure the Rscript command is available in your system's PATH.

    Please note that the `requirements.R` file should be in the same directory as `venv_setup.py` and should contain the names of the required R packages, one per line, with an indication whether it's from CRAN or Bioconductor.

    The `venv_setup.py` script will install these packages automatically when setting up the virtual environment.




#### Using `conda` or `mamba`

```bash
mamba create -n tls-env --no-default-packages -y
mamba activate tls-env
mamba env update -n tls-env --file environment.yml
```

## Project structure and usage

This project is structured as a Python package, which allows for better organization and reusability of the code. To use the project as a package:

1. Add the parent directory of `TLS` to your `PYTHONPATH` environment variable. This allows Python to find the `TLS` package. You can do this either temporarily for a single terminal session:

    ```bash
    export PYTHONPATH="/path/to/parent/of:$PYTHONPATH"
    ```

    or permanently for all terminal sessions by adding the line above to your shell's configuration file (e.g., `~/.bashrc`, `~/.zshrc`, `~/.bash_profile`, depending on your shell).

2. Ensure that every subdirectory in `TLS` that you want to use as a Python package contains an `__init__.py` file (it can be empty). This tells Python that the directory should be treated as a package.

3. You can then import the code in `TLS` as you would import any other Python package. For example:

    ```python
    from TLS.configs import config
    ```

Please remember to replace `/path/to/parent/of` with the actual path to the parent directory of `TLS`.

Note: If you are using Jupyter notebooks and encounter issues with imports, you may need to restart your notebook's kernel to pick up changes in your environment variables.

#### Using the interactive app

A suite of interactive data exploration tools are available as an app. [Networked computers at the MPIMG](http://mariux64.molgen.mpg.de/) can reach the app at [http://141.14.18.2:8501](http://141.14.18.2:8501).

The `streamlit` app can also be run locally from `reports/app` by cloning the `legacy-refactor` branch; however, the app requires intermediate data files that aren't yet included in the repository and must be provided by the user.


## Data availability

All scRNAseq datasets can be found in [the Gene Expression Omnibus (GEO) under accession code GSE141175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141175).