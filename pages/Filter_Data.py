import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



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
    mutation = pd.read_csv('./data/full_mutation.csv')
    mutation = mutation.set_index('protein')
    return mutation

abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

cls = list(abundance.columns)
tissues = set([str(i.split('_', 1)[1:][0]) for i in abundance.columns])

st.write("### View and filter protein data 🔬")

lineage_list = st.multiselect(
    'Filter for tissue',
     tissues, placeholder='Add tissues to view')

if lineage_list:
    cls = [s for s in cls if any(xs in s for xs in lineage_list)]

cl_list = st.multiselect(
    'Filter for cell-lines',
     cls, placeholder='Add cell-lines to view')

if cl_list:
    cls = cl_list

# --- NEW PASTE LOGIC START ---
st.write("Protein Quick Paste")
with st.expander("Paste a list of proteins"):
    # Text area for raw input
    pasted_proteins = st.text_area("Paste proteins (comma or newline separated)", 
                                  help="Example: P53, EGFR, BRAC1")
    
    if st.button("Populate Protein Selection"):
        if pasted_proteins:
            # Parse the input: split by commas or newlines and strip whitespace
            import re
            input_list = re.split(r'[,\n]+', pasted_proteins)
            input_list = [p.strip() for p in input_list if p.strip()]
            
            # Filter to ensure the pasted proteins actually exist in your index
            valid_proteins = [p for p in input_list if p in abundance.index]
            
            # Update session state
            current_selection = list(st.session_state['protein_selector'])
            updated_selection = list(set(current_selection + valid_proteins))
            
            st.session_state['protein_selector'] = updated_selection
            
            if len(valid_proteins) < len(input_list):
                missing = set(input_list) - set(valid_proteins)
                st.warning(f"Populated {len(valid_proteins)} proteins. Ignored {len(missing)} unknown IDs.")
            else:
                st.success(f"Added {len(valid_proteins)} proteins!")


# Initialize the session state key if it doesn't exist
if 'selected_proteins' not in st.session_state:
    st.session_state['selected_proteins'] = []


protein_list = st.multiselect(
    'Proteins to view',
    options=abundance.index,
    default=st.session_state['selected_proteins'], # This links it to your button
    key='protein_selector', # Giving it a unique key
    placeholder='Add proteins to view'
)

st.session_state['selected_proteins'] = protein_list



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


