import streamlit as st
import pandas as pd
import numpy as np
import re

from simko_utils import compute_tissue_variance_explained, get_tissue_of

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
            parsed_list = re.split(r'[,\n]+', pasted_input)
            parsed_list = [p.strip() for p in parsed_list if p.strip()]
            valid_pasted = [p for p in parsed_list if p in protein_options]

            current_list = st.session_state['saved_protein_list']
            for p in valid_pasted:
                if p not in current_list:
                    current_list.append(p)

            st.session_state['saved_protein_list'] = current_list
            st.session_state['temp_protein_list'] = current_list

            if len(valid_pasted) < len(parsed_list):
                missing = set(parsed_list) - set(valid_pasted)
                st.warning(f"Added {len(valid_pasted)} proteins. Ignored {len(missing)} invalid symbols.")
            else:
                st.success(f"Successfully appended {len(valid_pasted)} proteins to list (order preserved)!")

            st.rerun()

st.multiselect(
    'Selected Protein List',
    options=protein_options,
    default=st.session_state['saved_protein_list'],
    key='temp_protein_list',
    on_change=update_protein_list,
    placeholder='Search or select specific proteins to isolate'
)

st.caption(
    "This list also defines the 'chosen proteins' scope used for tracked-list statistics on the "
    "Abundance/Expression/Mutation pages (see Statistical Analysis Settings below) -- because it's "
    "set here, independently of any run's results, it's a valid pre-specified candidate list rather "
    "than a post-hoc selection."
)

# --- SECTION 3b: TISSUE-VARIANCE ANNOTATION ---
# A gene's own abundance can vary substantially by cancer type independently
# of whether it has any CRISPR fitness phenotype -- validation work found
# some KO targets (e.g. proliferation-linked genes) carry a lot of
# tissue-linked baseline variation, which can make pan-cancer results for
# that specific target harder to interpret without the tissue filter.
# Eta-squared (proportion of this gene's cell-line-to-cell-line variance
# explained by tissue-of-origin) is a fast, always-available warning signal
# for this -- not a p-value, just a 0-1 effect size to flag by eye.
@st.cache_data
def get_tissue_variance_table():
    abundance = pd.read_csv('data/abundance.csv').set_index('protein')
    tissue_of = get_tissue_of(abundance.columns)
    eta_sq = compute_tissue_variance_explained(abundance, tissue_of)
    return eta_sq

with st.expander("Tissue-variance check for your selected proteins", expanded=False):
    st.markdown(
        "How much of each selected protein's abundance varies by cancer type (tissue-of-origin), "
        "independent of essentiality or any KO-simulation result. **High values mean pan-cancer "
        "results for that protein are more likely to be tissue-confounded -- consider the tissue "
        "filter above, or interpret with extra caution.** Rough guide: >0.14 is a conventionally "
        "'large' effect size, 0.06-0.14 medium, below 0.06 small -- but there's no hard cutoff; use "
        "judgement."
    )
    eta_table = get_tissue_variance_table()
    proteins_to_check = list(dict.fromkeys(
        st.session_state['saved_ko_proteins'] + st.session_state['saved_protein_list']
    ))
    if proteins_to_check:
        rows = []
        for p in proteins_to_check:
            role = []
            if p in st.session_state['saved_ko_proteins']:
                role.append('KO target')
            if p in st.session_state['saved_protein_list']:
                role.append('tracked')
            value = eta_table.get(p, None)
            rows.append({
                'Protein': p,
                'Role': ' + '.join(role),
                'Tissue variance (η²)': round(value, 3) if value is not None and not pd.isna(value) else 'not covered in Abundance panel',
            })
        st.dataframe(pd.DataFrame(rows).sort_values('Protein'), hide_index=True, use_container_width=True)
    else:
        st.caption("Select KO target(s) and/or a protein list above to see their tissue-variance values here.")

# --- SECTION 4: STATISTICAL ANALYSIS SETTINGS ---
st.subheader("4. Statistical Analysis Settings")
st.markdown(
    "These apply across the Abundance, Expression and Mutation pages, which no longer show their own "
    "per-page statistics controls -- set everything here once."
)

with st.expander("Cohort combination mode", expanded=False):
    st.radio(
        "When multiple KO targets are selected:",
        options=[
            "Strict — cell line must qualify on every target (recommended for double-knockout / genetic-interaction analyses)",
            "Compensatory — a cell line can qualify via one target alone (partial/average across targets)",
        ],
        index=0,
        key='saved_strict_multi_target_choice',
        help=(
            "Strict maps to rank-percentile MAX for Abundance/Expression (a cell line only counts as "
            "low if it's low in every target) and to 'all targets mutated' for Mutation. Compensatory "
            "maps to averaged abundance/expression, or 'any target mutated', respectively. Has no "
            "effect with a single target selected."
        ),
    )
st.session_state['saved_strict_multi_target'] = st.session_state.get(
    'saved_strict_multi_target_choice', ''
).startswith('Strict')

with st.expander("'Most changed' table — genome-wide scan", expanded=True):
    st.markdown(
        "**Effect-size floor.** At typical cohort sizes here, genome-wide FDR correction leaves few or "
        "no hits even for well-established effects -- see Design and Implementation / Discussion for "
        "why. Filtering the genome-wide scan by a minimum |differential| before correction is the "
        "tool's current default approach; it is NOT formally independent of the p-value the way a "
        "covariate-based filter would be (Bourgon, Gentleman & Huber, PNAS 2010), so treat resulting "
        "q-values as indicative rather than a strict FDR guarantee if reporting them."
    )
    st.radio(
        "Restrict the genome-wide correction to:",
        options=[
            "Features with |diff| ≥ threshold (default — not formally independent of p, see above)",
            "Features with overall variability ≥ threshold (independent filter, DESeq2-style — formally justified, but usually removes far fewer features)",
            "All genome-wide features (no restriction — most conservative, may show few or no hits)",
        ],
        index=0,
        key='saved_correction_scope_choice',
    )
    c1, c2, c3 = st.columns(3)
    st.session_state['saved_correction_threshold_ab'] = c1.number_input(
        "Abundance threshold", min_value=0.0, value=st.session_state.get('saved_correction_threshold_ab', 0.25), step=0.05, format="%.3f"
    )
    st.session_state['saved_correction_threshold_ex'] = c2.number_input(
        "Expression threshold", min_value=0.0, value=st.session_state.get('saved_correction_threshold_ex', 0.25), step=0.05, format="%.3f"
    )
    st.session_state['saved_correction_threshold_dep'] = c3.number_input(
        "Dependency threshold", min_value=0.0, value=st.session_state.get('saved_correction_threshold_dep', 0.25), step=0.05, format="%.3f"
    )

    st.markdown("**Display filter for the 'most changed' table itself** (applied on top of the above).")
    st.session_state['saved_diff_floor'] = st.number_input(
        "Minimum |differential| to display", min_value=0.0,
        value=st.session_state.get('saved_diff_floor', 0.25), step=0.05, format="%.3f",
    )
    st.checkbox(
        "Also require a significance threshold",
        value=st.session_state.get('saved_sig_filter_enabled', False),
        key='saved_sig_filter_enabled',
        help="Off by default -- with the effect-size floor above already doing the filtering, an "
             "additional significance requirement is optional, not the primary filter.",
    )
    if st.session_state['saved_sig_filter_enabled']:
        cc1, cc2 = st.columns(2)
        sig_choice = cc1.radio(
            "Significance column", options=['Nominal p (uncorrected)', 'Adjusted q'], index=0,
            key='saved_sig_col_choice', horizontal=True,
        )
        st.session_state['saved_sig_col'] = 'p' if sig_choice.startswith('Nominal') else 'p_adj'
        st.session_state['saved_sig_threshold'] = cc2.number_input(
            "Threshold", min_value=0.0, max_value=1.0,
            value=st.session_state.get('saved_sig_threshold', 0.01), step=0.001, format="%.4f",
        )
    st.radio(
        "Order the table by",
        options=['Differential magnitude (largest changes first)', 'Significance (lowest p/q first)'],
        index=0, key='saved_sort_by_choice', horizontal=True,
    )
    st.session_state['saved_sort_by'] = (
        'diff' if st.session_state.get('saved_sort_by_choice', '').startswith('Differential') else 'significance'
    )

with st.expander("'Chosen proteins' table — tracked list", expanded=False):
    st.markdown(
        "The Abundance/Expression/Mutation pages' 'chosen proteins' section (your Selected Protein "
        "List from Section 3, plus KO targets) now computes its own, separate correction restricted to "
        "exactly that list -- m = however many of those genes are present in each matrix, rather than "
        "the whole genome. Because that list is specified here, independently of any analysis run, "
        "this correction carries the full formal FDR guarantee (unlike the diff-prefiltered genome-wide "
        "scan above). No further settings needed -- this always uses your current Selected Protein List."
    )

st.success("Configuration saved! Your settings are locked in and ready across all tabs.")