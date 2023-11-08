import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(
    page_title="SimKO - Abundance data",
    page_icon="🥼",
)

@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('data/abundance.csv')
    abundance = abundance.set_index('protein')
    return abundance

abundance = get_abundance_data()
cmap = plt.cm.get_cmap('RdYlBu_r')


protein_list = st.multiselect(
    'Proteins for KO',
     abundance.index, placeholder='Add proteins to analyse')

cl_list = st.multiselect(
    'CLs for KO',
     abundance.columns, placeholder='Add CLs to view')


if protein_list:
    abundance_filter = abundance.loc[abundance.index.isin(protein_list)]

    if cl_list:
        abundance_view = abundance_filter.filter(cl_list)
    else:
        abundance_view = abundance_filter

    st.dataframe(abundance_view.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")


else:
    st.write('Please select at least one protein.')


