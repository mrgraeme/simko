import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(
    page_title="SimKO - Filter data",
    page_icon="🥼",
)

@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('data/abundance.csv')
    abundance = abundance.set_index('protein')
    return abundance

@st.cache_data
def get_expression_data():
    expression = pd.read_csv('./data/expression.csv')
    expression = expression.set_index('protein')
    return expression

@st.cache_data
def get_mutation_data():
    mutation = pd.read_csv('./data/mutation.csv')
    mutation = mutation.set_index('protein')
    return mutation


abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
cmap = plt.cm.get_cmap('RdYlBu_r')


protein_list = st.multiselect(
    'Proteins for KO',
     abundance.index, placeholder='Add proteins to analyse')

cl_list = st.multiselect(
    'CLs for KO',
     abundance.columns, placeholder='Add CLs to view')


if protein_list:
    abundance_filter = abundance.loc[abundance.index.isin(protein_list)]
    expression_filter = expression.loc[expression.index.isin(protein_list)]
    mutation_filter = mutation.loc[mutation.index.isin(protein_list)]

    if cl_list:
        abundance_filter = abundance_filter.filter(cl_list)
        expression_filter = expression_filter.filter(cl_list)
        mutation_filter = mutation_filter.filter(cl_list)

    tab1, tab2, tab3 = st.tabs(["abundance", "expression", "mutation"])

    with tab1:
        st.dataframe(abundance_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")
    with tab2:
        st.dataframe(expression_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")
    with tab3:
        st.dataframe(mutation_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")






else:
    st.write('Please select at least one protein.')


