import matplotlib.pyplot as plt
import numpy as np # TODO: Unclear why this is needed
import scanpy as sc
import streamlit as st

from TLS.configs.config_manager import config
from TLS.data.assets.assets_manager import data_assets
from TLS.reports.app.app_utils import page_footer
from TLS.src.tls_utils import ignore_warnings


if __name__ == '__main__':
    st.set_page_config(page_title='Pseudotime-Ordered Gene Expression', page_icon='🖼️')
    ignore_warnings()
    page_footer(
        page_title='Pseudotime-Ordered Gene Expression',
        page_text="""
        Clustered, scaled gene expression values from scRNA-seq for
        a given gene set from the 120 hr time point. Cells are ordered by pseudotime,
        estimated using the diffusion pseudotime method of [Haghverdi16].
        """
    )
    
    user_input = st.text_input("Enter a gene set name here and press Enter:")
    st.markdown(
    """
    **Gene set options**:
    - `g2_genes`
    - `hoxa_genes`
    - `hoxb_genes`
    - `hoxc_genes`
    - `hoxd_genes`
    - `neural_trajectory_genes`
    - `somitic_trajectory_genes`
    """
    )

    # Prepare 120 hr data
    tls_adata = sc.read(config.data['TLS AnnData']['path'], sparse=True)    
    tls_adata_120h = (
        tls_adata
        .tls
        .query_timepoints(['Organoid_120h'])
        .paga(groups='louvain')
        .louvain(restrict_to=['louvain', ['Seurat_1', 'Seurat_3', 'Seurat_8', 'Seurat_0']], resolution=0, key_added='louvain2')
        .paga(groups='louvain2')
        .assign_uns(key='iroot', value_func=lambda ad: np.flatnonzero(ad.obs['louvain2'] == 'Seurat_4')[0])
        .dpt() # Adds observation 'dpt_pseudotime'
        .assign_obs(key='iroot', value_func=lambda ad: ad.obs['dpt_pseudotime'])
        .scale()
    )

    if user_input:
        st.write(f'You entered: {user_input}')

        # Plot data
        fig, ax = plt.subplots(1, 1, figsize=(4, 4),
                               # Must be false for `paga_path` to work
                               constrained_layout=False)
        try:
            gene_set = getattr(data_assets, user_input)
        except AttributeError:
            st.write(f'Gene set `{user_input}` not found in the data.')
            raise
        (
            tls_adata_120h
            # Withouth assigning distance, `paga_path` throws an error because of a type error (??)
            .assign_obs(key='distance', value_func=lambda ad: ad.obs['dpt_pseudotime'])
            .paga_path(
                nodes=["Seurat_4", "Seurat_7", "Seurat_2", "Seurat_5", "Seurat_1-Seurat_3-Seurat_8-Seurat_0,0"],
                keys=gene_set,
                show_node_names=False,
                ytick_fontsize=12,
                left_margin=0.5,
                annotations=['distance'],
                show_colorbar=True,
                normalize_to_zero_one=True,
                return_data=True,
                ax=ax,
                show=False,
            )
        )
        ax.set_frame_on(False)
        fig.set_size_inches(6, len(gene_set)/4)
    
        # Show in app
        st.pyplot(fig)