from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np # TODO: Unclear why this is needed
from PIL import Image
import scanpy as sc
import streamlit as st

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
        # Not sure why this works ... can't pass the AnnData directly because it's unhashable
        tls_adata
        .tls
        .umap(ax=ax, color=gene)
    )
    cbar = _get_cbar(fig)
    cbar.set_title('Expression,\nscaled', fontsize=9)

    return fig

def _get_cbar(fig):
    cbar = [ax for ax in fig.get_axes() if ax.get_label() == '<colorbar>']
    try:
        cbar = cbar[0]
        return cbar
    except:
        # No colorbar
        pass

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
        st.write(f'You entered: _:blue[{user_input}]_')
        try:
            gene_umap_fig = make_gene_umap(user_input)
            cell_type_ref_img = Image.open('./reports/figures/cell_type_clusters_umap.png')
        
            expr_col, ref_col = st.columns(2)
            with expr_col:
                st.markdown('Gene expression:')
                st.pyplot(gene_umap_fig)

            with ref_col:
                st.markdown('Cell type reference:')
                st.image(cell_type_ref_img, width=400)

        except KeyError:
            st.write(f'Gene `{user_input}` not found in the data.')