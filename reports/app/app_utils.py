from matplotlib.figure import Figure
import streamlit as st


def page_footer(
    page_title: str,
    page_text: str = None,
    abbr_title: str = None,
    fig: Figure = None,
) -> None:
    # Check inputs
    if abbr_title is None:
        abbr_title = page_title
    
    # Render text
    st.markdown(f'# {abbr_title}')
    st.sidebar.header(abbr_title)
    st.write(page_text)

    # Render figure    
    if fig is not None:
        st.pyplot(fig=fig)