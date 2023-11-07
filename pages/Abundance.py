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
    abundance = pd.read_csv('data/CCLE_imputed_all.csv')
    abundance = abundance.rename(columns={'Gene':'Protein'})
    return abundance

# protein_list = ['ARID1A', 'PBRM1', 'BRAF']
# get_classes_by_mean_abundance(protein_list, abundance)


abundance = get_abundance_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

# x = st.slider('x', min_value=10, max_value=50)  # 👈 this is a widget
# st.write(x, 'squared is', x * x)
# y=st.slider('y', min_value=0.1, max_value=1.0)  # 👈 this is a widget

# chart_data = chart_data[0:x]
# chart_data = chart_data * y

protein_list = st.multiselect(
    'Proteins for KO',
     abundance['Protein'], placeholder='Add proteins to analyse')


cl_list = st.multiselect(
    'CLs for KO',
     abundance.columns, placeholder='Add CLs to view')


if protein_list:
    abundance_filter = abundance.loc[abundance['Protein'].isin(protein_list)]

    if cl_list:
        abundance_view = abundance_filter.filter(['Protein'] + cl_list)
    else:
        abundance_view = abundance_filter

    st.dataframe(abundance_view.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None))


else:
    st.write('Please select at least one protein.')


