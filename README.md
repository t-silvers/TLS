# Mouse embryonic stem cells self-organize into trunk-like structures with neural tube and somites

## Abstract

> Post-implantation embryogenesis is a highly dynamic process comprising multiple lineage decisions and morphogenetic changes inaccessible to deep analysis in vivo. Here, we show that pluripotent mouse embryonic stem cells (mESCs) form aggregates that upon embedding in an extra-cellular matrix compound induce the formation of highly organized "Trunk-Like-Structures (TLS)" comprising the neural tube and somites. Comparative single-cell RNA-seq analysis confirms that this process is highly analogous to mouse development and follows the same step-wise gene-regulatory program. Tbx6 knockout TLS develop additional neural tubes mirroring the embryonic mutant phenotype, and chemical modulation can induce excess somite formation. Trunk-like-structures thus reveal an unprecedented level of self-organization and provide a powerful platform for investigating post-implantation embryogenesis in a dish.

## Code Availability

This repository contains the collection of R and python snippets and scripts used to perform the transcriptomics analysis presented in Veenvliet and Bolondi et al. The main processing and data integration steps are show, based on which general visualization and plots of e.g. single genes were done.
A more detailled description and references for used tools can be found in the [corresponding publication](https://www.biorxiv.org/content/10.1101/2020.03.04.974949v1).

### Setup Virtual Environment


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

## Project Structure and Usage

This project is structured as a Python package, which allows for better organization and reusability of the code. To use the project as a package:

1. Add the parent directory of `TLS` to your `PYTHONPATH` environment variable. This allows Python to find the `TLS` package. You can do this either temporarily for a single terminal session:

    ```bash
    export PYTHONPATH="/path/to/parent/of:$PYTHONPATH"
    ```

    or permanently for all terminal sessions by adding the line above to your shell's configuration file (e.g., `~/.bashrc`, `~/.zshrc`, `~/.bash_profile`, depending on your shell).

2. Ensure that every subdirectory in `TLS` that you want to use as a Python package contains an `__init__.py` file (it can be empty). This tells Python that the directory should be treated as a package.

3. You can then import the code in `TLS` as you would import any other Python package. For example:

    ```python
    from myproj.configs import config
    ```

Please remember to replace `/path/to/parent/of` with the actual path to the parent directory of `TLS`.

Note: If you are using Jupyter notebooks and encounter issues with imports, you may need to restart your notebook's kernel to pick up changes in your environment variables.

#### Using the interactive app

An interactive data exploration tools is available as a `streamlit` app at [app.com](app.com).


## Data Availability

All scRNAseq datasets can be found in [the Gene Expression Omnibus (GEO) under accession code GSE141175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141175).