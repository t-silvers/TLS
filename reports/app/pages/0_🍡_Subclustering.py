import streamlit as st

from TLS.reports.app.app_utils import page_footer
from TLS.scripts.pseudotime.make_subclustering_plots import make_subclustering_plots
from TLS.src.tls_utils import ignore_warnings


if __name__ == '__main__':
    st.set_page_config(page_title='Subclustering TLS Cells', page_icon='🍡')
    ignore_warnings()
    
    subclustering_fig = make_subclustering_plots()

    page_footer(
        page_title='Subclustering TLS Cells',
        page_text="""
            With the Seurat package, we can cluster cells based on their
            gene expression profiles using the Louvain method.
        """,
        abbr_title='scRNA-seq Clustering',
        fig=subclustering_fig,
    )