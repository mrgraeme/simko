import streamlit as st

st.set_page_config(
    page_title="SimKO",
    page_icon="🧫",
)

st.write("# Welcome to SimKO! 🧫")

# st.sidebar.success("Select a demo above.")

st.markdown(
    """
    SimKO is a tool to interrogate and visualise protein abundance change in celllines.
    
    **👈 Select a tool from the sidebar** 
    ### Want to learn more?
    - Check out [pkg.rcrds.com](https://pkg.rcrds.com)
"""
)
