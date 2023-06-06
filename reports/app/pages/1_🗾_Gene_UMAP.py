import streamlit as st

import matplotlib.pyplot as plt
import numpy as np # TODO: Unclear why this is needed
import scanpy as sc

from TLS.configs.config_manager import config
from TLS.reports.app.app_utils import page_footer
from TLS.src.tls_utils import ignore_warnings


@st.cache_resource
def get_tls_adata():
    # TODO: Factor out to src
    return sc.read(config.data['TLS AnnData']['path'], sparse=True, cache=True)

@st.cache_resource
def make_gene_umap(gene):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
    (
        tls_adata
        .tls
        .umap(ax=ax, color=gene)
    )

    return fig


if __name__ == '__main__':
    st.set_page_config(page_title='Gene Expression in a 2D Embedding', page_icon='🗾')
    ignore_warnings()

    page_footer(
        page_title='Gene Expression in a 2D Embedding',
        page_text="""
        2D Embeddings (e.g., UMAP) colored by expression of indicated genes.
        """,
        abbr_title='Expression Embedded',
    )
    
    tls_adata = get_tls_adata()
    
    user_input = st.text_input("Enter a gene name here (e.g., `Cdx2`) and press Enter:")
    if user_input:
        st.write(f'You entered: {user_input}')
        try:
            gene_umap_fig = make_gene_umap(user_input)
            st.pyplot(gene_umap_fig)
        
        except KeyError:
            st.write(f'Gene `{user_input}` not found in the data.')