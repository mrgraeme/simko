import streamlit as st
import pandas as pd
import numpy as np
import re

# 1. Initialize persistent "Shadow Keys" that Streamlit won't delete
if 'saved_ko_proteins' not in st.session_state:
    st.session_state['saved_ko_proteins'] = []
if 'saved_tissues' not in st.session_state:
    st.session_state['saved_tissues'] = []
if 'saved_cell_lines' not in st.session_state:
    st.session_state['saved_cell_lines'] = []
if 'saved_protein_list' not in st.session_state:
    st.session_state['saved_protein_list'] = []

@st.cache_data
def get_abundance_options():
    abundance = pd.read_csv('data/abundance.csv')
    abundance = abundance.set_index('protein')
    cls = list(abundance.columns)
    tissues = sorted(list(set([str(i.split('_', 1)[1:][0]) for i in abundance.columns])))
    return list(abundance.index), cls, tissues

protein_options, all_cell_lines, tissue_options = get_abundance_options()

st.title("Global Analysis Control Panel ⚙️")
st.markdown("Configure your dataset scopes here. These settings automatically persist across all analysis tabs.")

# --- CALLBACK FUNCTIONS TO PERSIST DATA ---
def update_ko_proteins():
    st.session_state['saved_ko_proteins'] = st.session_state['temp_ko_proteins']

def update_tissues():
    st.session_state['saved_tissues'] = st.session_state['temp_tissues']
    # Clear cell lines if they no longer match selected tissues
    if st.session_state['temp_tissues']:
        valid_cls = [cl for cl in all_cell_lines if any(ts in cl for ts in st.session_state['temp_tissues'])]
        st.session_state['saved_cell_lines'] = [cl for cl in st.session_state['saved_cell_lines'] if cl in valid_cls]

def update_cell_lines():
    st.session_state['saved_cell_lines'] = st.session_state['temp_cell_lines']

def update_protein_list():
    st.session_state['saved_protein_list'] = st.session_state['temp_protein_list']


# --- SECTION 1: PROTEINS FOR KO ---
st.subheader("1. Knockout Targets")
st.multiselect(
    'Proteins for KO',
    options=protein_options,
    default=st.session_state['saved_ko_proteins'],
    key='temp_ko_proteins',
    on_change=update_ko_proteins,
    placeholder='Select targeted KO proteins'
)

# --- SECTION 2: COHORT FILTERING (TISSUE & CELL LINE) ---
st.subheader("2. Cohort Filtering")
col1, col2 = st.columns(2)

with col1:
    st.multiselect(
        'Tissue Filter',
        options=tissue_options,
        default=st.session_state['saved_tissues'],
        key='temp_tissues',
        on_change=update_tissues,
        placeholder='Filter by broader tissue groups'
    )

# Dynamically constrain available cell lines based on chosen tissues
available_cell_lines = all_cell_lines
if st.session_state['saved_tissues']:
    available_cell_lines = [cl for cl in all_cell_lines if any(ts in cl for ts in st.session_state['saved_tissues'])]

with col2:
    st.multiselect(
        'Cell Line Filter',
        options=available_cell_lines,
        default=st.session_state['saved_cell_lines'],
        key='temp_cell_lines',
        on_change=update_cell_lines,
        placeholder='Defaults to all available'
    )

# --- SECTION 3: COMPREHENSIVE PROTEIN LIST & PASTE ---
st.subheader("3. Downstream Protein Focus List")

with st.expander("💡 Quick Paste Proteins (Comma or Newline Separated)"):
    pasted_input = st.text_area("Paste tracking list IDs here:", help="Example: P53, EGFR, BRCA1")
    
if st.button("Append to Protein List"):
        if pasted_input:
            # 1. Parse the text into clean strings
            parsed_list = re.split(r'[,\n]+', pasted_input)
            parsed_list = [p.strip() for p in parsed_list if p.strip()]
            valid_pasted = [p for p in parsed_list if p in protein_options]
            
            # 2. FIX: Deduplicate while explicitly preserving the structural sequence order
            current_list = st.session_state['saved_protein_list']
            for p in valid_pasted:
                if p not in current_list:
                    current_list.append(p)
            
            # 3. Synchronize both state keys to avoid the linear lifecycle trap
            st.session_state['saved_protein_list'] = current_list
            st.session_state['temp_protein_list'] = current_list
            
            if len(valid_pasted) < len(parsed_list):
                missing = set(parsed_list) - set(valid_pasted)
                st.warning(f"Added {len(valid_pasted)} proteins. Ignored {len(missing)} invalid symbols.")
            else:
                st.success(f"Successfully appended {len(valid_pasted)} proteins to list (order preserved)!")
            
            # 4. Force immediate screen redraw
            st.rerun()

st.multiselect(
    'Selected Protein List',
    options=protein_options,
    default=st.session_state['saved_protein_list'],
    key='temp_protein_list',
    on_change=update_protein_list,
    placeholder='Search or select specific proteins to isolate'
)

st.success("Configuration saved! Your settings are locked in and ready across all tabs.")