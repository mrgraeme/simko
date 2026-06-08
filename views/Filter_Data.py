import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('data/abundance.csv')
    return abundance.set_index('protein')

@st.cache_data
def get_expression_data():
    expression = pd.read_csv('./data/expression.csv')
    return expression.set_index('protein')

@st.cache_data
def get_mutation_data():
    mutation = pd.read_csv('./data/full_mutation.csv')
    return mutation.set_index('protein')

# Load raw matrices
abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

st.write("### View and filter protein data 🔬")

# 1. Inherit Global Parameters from Configuration State

ko_targets = st.session_state.get('saved_ko_proteins', [])
active_proteins = st.session_state.get('saved_protein_list', []) 
active_tissues = st.session_state.get('saved_tissues', [])
active_cell_lines = st.session_state.get('saved_cell_lines', [])


# 2. Determine active cohort columns based on global scope
all_cols = list(abundance.columns)
if active_cell_lines:
    active_cls = active_cell_lines
elif active_tissues:
    active_cls = [s for s in all_cols if any(xs in s for xs in active_tissues)]
else:
    active_cls = all_cols

# 3. Render Inspector Tables
if active_proteins:
    st.caption(f"Showing matrix view for {len(active_proteins)} focus proteins across {len(active_cls)} matched cell lines.")
    
    tab1, tab2, tab3 = st.tabs(["Abundance", "Expression", "Mutation"])
    
    with tab1:
        abundance_filter = abundance.loc[abundance.index.isin(active_proteins)].filter(active_cls)
        st.dataframe(abundance_filter.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))
    
    with tab2:
        expression_filter = expression.loc[expression.index.isin(active_proteins)].filter(active_cls)
        st.dataframe(expression_filter.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))
    
    with tab3:
        mutation_filter = mutation.loc[mutation.index.isin(active_proteins)].filter(active_cls)
        st.dataframe(mutation_filter.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))
else:
    st.info("💡 No focus proteins found in global state. Go to **Global Settings** to select or paste a list of proteins to display here.")