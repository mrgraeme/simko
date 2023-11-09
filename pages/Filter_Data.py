import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(
    page_title="SimKO - Filter data",
    page_icon="🥼",
    layout="wide"
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

cls = list(abundance.columns)
tissues = set([str(i.split('_', 1)[1:][0]) for i in abundance.columns])

st.write("### View and filter protein data! 🥼")



lineage_list = st.multiselect(
    'Filter for tissue',
     tissues, placeholder='Add tissues to view')

if lineage_list:
    cls = [s for s in cls if any(xs in s for xs in lineage_list)]

cl_list = st.multiselect(
    'cell-lines to view',
     cls, placeholder='Add cell-lines to view')

if cl_list:
    cls = cl_list

protein_list = st.multiselect(
    'Proteins to view (Tip: Filter for cell-lines / lineage first to make this run faster!)',
     abundance.index, placeholder='Add proteins to view')



if protein_list:

    tab1, tab2, tab3 = st.tabs(["abundance", "expression", "mutation"])

    with tab1:
        abundance_filter = abundance.loc[protein_list].filter(cls)
        st.dataframe(abundance_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")

    with tab2:
        expression_filter = expression.loc[expression.index.isin(protein_list)].filter(cls)
        st.dataframe(expression_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")
    
    with tab3:
        mutation_filter = mutation.loc[mutation.index.isin(protein_list)].filter(cls)
        st.dataframe(mutation_filter.style.background_gradient(cmap=cmap,vmin=(-6),vmax=6,axis=None).format("{:.3f}"), ) # .format("{:.2%}")






else:
    st.write('Please select at least one protein.')


