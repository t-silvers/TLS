import streamlit as st

from TLS.reports.app.app_utils import page_footer


@st.cache_data
def get_data():
    # TODO: Allow concurrent access to the file
    return []


if __name__ == '__main__':
    st.set_page_config(page_title='Feedback', page_icon='🎙️')
    # APP
    page_footer(
        page_title='Feedback on the TLS App',
        page_text='Provide feedback on the TLS app here.',
        abbr_title='Feedback',
    )

    # Declare a form and call methods directly on the returned object
    form = st.form(key='feedback_form')
    user_feedback = form.text_input(label='Feedback:')
    submit_button = form.form_submit_button(label='Submit')
    
    if submit_button:
        with open('user_feedback.txt', 'a') as f:
            f.write(user_feedback + '\n')
        
        st.write('Feedback submitted!')
        st.markdown(f'_:blue[{user_feedback}]_')
        
    st.markdown(
    """
    **Known issues**:
    - Label axes and add legend(s) to plots
    - Replace clustering method and labels with meaningful ones (i.e., cell type)
    - `3_📊_DataFrame_Demo` is a placeholder page
    """
    )

    st.markdown(
    """
    **Roadmap**:
    - Allow exploration of the data in `3_📊_DataFrame_Demo`
    - Add a page for more of the plots in the paper
    - Accept user input for the gene set
    - Accept user input for the timepoint(s)
    - Accept user input for the cluster method and cluster(s) rendered
    - Accept user input for plot aesthetics (palette, size, etc.)
    - Allow concurrency for the feedback form
    - More flexibility (drop-down menu?) for gene input
    """
    )

    st.markdown('_:red[Report any other issues you find here]_')