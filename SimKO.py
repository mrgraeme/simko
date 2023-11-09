import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


st.set_page_config(
    page_title="SimKO",
    page_icon="🥼",
    layout="wide"
)

# protein_list = ['ARID1A', 'PBRM1', 'BRAF']
# get_classes_by_mean_abundance(protein_list, abundance)


abundance = get_abundance_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

st.write("# Welcome to SimKO! 🥼")

# st.sidebar.success("Select a demo above.")

st.markdown(
    """
    SimKO is a tool to interrogate and visualise protein abundance change in celllines.
    
    **👈 Select a tool from the sidebar** 
    ### Want to learn more?
    - Check out [pkg.rcrds.com](https://pkg.rcrds.com)
"""
)
