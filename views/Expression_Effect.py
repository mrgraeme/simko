import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns

@st.cache_data
def get_abundance_data():
    return pd.read_csv('./data/abundance.csv').set_index('protein')

@st.cache_data
def get_expression_data():
    return pd.read_csv('./data/expression.csv').set_index('protein')

@st.cache_data
def get_mutation_data():
    return pd.read_csv('./data/full_mutation.csv').set_index('protein')

# Math and analytical utilities
def get_classes_by_mean_expression(protein_list, expression, n):
    protein_list_expression = expression.loc[expression.index.isin(protein_list)].T
    protein_list_expression['mean'] = protein_list_expression.mean(axis=1)
    protein_list_expression = protein_list_expression.sort_values("mean", ascending=False)
    
    median = protein_list_expression.tail(round((protein_list_expression.shape[0]/2) + n/2)).head(n)
    median['class'] = 'median'
    low = protein_list_expression.tail(n)
    low['class'] = 'low'
    
    return pd.concat([median, low])

def ttest_from_sample_stats(row, n_cls=20):
    pooled_sd = np.sqrt((((n_cls-1)*(row['low_std']**2)) + ((n_cls-1)*(row['median_std']**2))) / (n_cls + n_cls - 2))
    if pooled_sd == 0:
        return 1.0
    t_stat = (row['low'] - row['median']) / (pooled_sd * np.sqrt(1/n_cls + 1/n_cls))
    return 2 * (1 - t.cdf(abs(t_stat), (n_cls + n_cls - 2)))

def get_differentials(class_df, data_df, n):
    median_class = list(class_df.loc[class_df['class']=='median'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    
    diff_df = pd.DataFrame()
    diff_df['median'] = data_df.filter(median_class).mean(axis=1)
    diff_df['median_std'] = data_df.filter(median_class).std(axis=1)
    diff_df['low'] = data_df.filter(low_class).mean(axis=1)
    diff_df['low_std'] = data_df.filter(low_class).std(axis=1)
    
    diff_df['diff'] = diff_df['low'] - diff_df['median']
    diff_df['p'] = diff_df.apply(ttest_from_sample_stats, n_cls=n, axis=1)
    return diff_df.drop(columns=['low_std', 'median_std']).sort_values('diff', ascending=True)

def get_differentials_boxplot(class_df, data_df, protein_list, n):
    # 1. Grab the global protein order straight from session state configuration
    # This ensures it exactly mirrors the order you picked or pasted them in.
    global_order = st.session_state.get('saved_protein_list', [])
    
    # Filter global order to only include proteins actually present in the data slice
    plot_order = [p for p in global_order if p in protein_list]
    
    # Fallback to current slice if global list isn't populated yet
    if not plot_order:
        plot_order = protein_list

    # 2. Reshape and merge the matrix data
    data_df = data_df.reset_index()
    data_df = data_df.loc[data_df['protein'].isin(plot_order)]
    data_df = data_df.melt(id_vars='protein')
    
    class_df = class_df[['class']].reset_index()
    box_data = data_df.merge(class_df, how='left', left_on='variable', right_on='index').dropna()
    
    # 3. Define standard ggplot2-inspired aesthetic colors
    # Hex approximations for ggplot2 defaults: Cornflower Blue-ish (#619CFF) & Rose Pink (#F8766D)
    ggplot2_palette = {
        'median': '#619CFF',
        'low': '#F8766D'
    }
    
    fig = plt.figure(figsize=(10, 4))
    
    # 4. Enforce order and palette strictly
    sns.boxplot(
        data=box_data, 
        x="protein", 
        y="value", 
        hue='class', 
        order=plot_order,                 # Force your configuration page sorting sequence
        hue_order=['median', 'low'],      # Lock category hierarchy
        palette=ggplot2_palette           # Apply constant color map
    )
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    return fig

def process_mutations(class_df, data_df):
    median_class = list(class_df.loc[class_df['class']=='median'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    diff_df = pd.DataFrame()
    diff_df['median'] = data_df.filter(median_class).sum(axis=1)
    diff_df['low'] = data_df.filter(low_class).sum(axis=1)
    diff_df['diff'] = diff_df['low'] - diff_df['median']
    return diff_df.sort_values('diff', ascending=True)

def get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, protein_list):
    ab_sub = diff_abund_df.loc[diff_abund_df.index.isin(protein_list)].copy()
    ab_sub.columns = ['Abundance - ' + n for n in ab_sub.columns]
    
    ex_sub = diff_exp_df.loc[diff_exp_df.index.isin(protein_list)].copy()
    ex_sub.columns = ['Expression - ' + n for n in ex_sub.columns]
    
    mu_sub = diff_mut_df.loc[diff_mut_df.index.isin(protein_list)].copy()
    mu_sub.columns = ['Mutation - ' + n for n in mu_sub.columns]
    
    return pd.concat([ab_sub, ex_sub, mu_sub], axis=1)

# Core Execution Pipeline
abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

st.write("### Explore expression change effects 🫶")

# 1. Grab global properties
ko_targets = st.session_state.get('saved_ko_proteins', [])
focus_proteins = st.session_state.get('saved_protein_list', [])
active_tissues = st.session_state.get('saved_tissues', [])
active_cell_lines = st.session_state.get('saved_cell_lines', [])

# 2. Slice to cohort
all_cols = list(abundance.columns)
if active_cell_lines:
    active_cls = active_cell_lines
elif active_tissues:
    active_cls = [s for s in all_cols if any(xs in s for xs in active_tissues)]
else:
    active_cls = all_cols

abundance = abundance.filter(active_cls)
expression = expression.filter(active_cls)
mutation = mutation.filter(active_cls)

if ko_targets:
    st.caption(f"Analyzing downstream markers based on {ko_targets} background expressions.")
    n = ((abundance.shape[1]-1) // 3) if abundance.shape[1] < 60 else 20
    class_df = get_classes_by_mean_expression(ko_targets, expression, n)
    
    tab1, tab2 = st.tabs(["Expression Slices", "Cohort Partition Matrix"])
    with tab1:
        fig = plt.figure(figsize=(10, 3))
        sns.heatmap(class_df[['mean']].sort_values('mean').T.astype(float).round(1), square=True, cmap="vlag", annot=True, annot_kws={'size': 6}, cbar=False)
        st.pyplot(fig)
    with tab2:
        st.table(class_df.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None))

    st.markdown("**Mean abundance, expression, and mutation count updates**")

    protein_select = st.multiselect(
        'Select additional proteins differences to view',
        options=abundance.index, placeholder='Add additional markers'
    )

    show_proteins = list(set(ko_targets + focus_proteins + protein_select))

    diff_abund_df = get_differentials(class_df, abundance, n)
    diff_exp_df = get_differentials(class_df, expression, n)
    diff_mut_df = process_mutations(class_df, mutation)

    if show_proteins:
        diff_summary = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, show_proteins)
        fig_box = get_differentials_boxplot(class_df, abundance, show_proteins, n)
        st.pyplot(fig_box)
        st.dataframe(diff_summary.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))

    st.markdown("**Proteins most significantly affected by background expression shifts**")
    n_outlayers = st.slider('Number proteins by median differential', value=5, min_value=0, max_value=50)
 
    significant_hits = diff_abund_df[diff_abund_df['p'] < 0.01]
    diff_proteins = list(significant_hits.head(n_outlayers).index) + list(significant_hits.tail(n_outlayers).index)
    diff_top = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_proteins)

    tab_fc, tab_ab = st.tabs(["Foldchange Highlights", "Raw Abundances"])
    with tab_fc:
        st.dataframe(diff_top.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))
        st.download_button('Download Foldchange Data', diff_abund_df.to_csv(index=True).encode('utf-8'), f"expression_foldchange_{'_'.join(ko_targets)}.csv", "text/csv")
        
    with tab_ab:
        abundance_out = abundance.loc[abundance.index.isin(diff_abund_df.index)]
        st.dataframe(abundance_out[class_df.index].head(100).style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))
        st.download_button('Download Raw Abundance Matrix', abundance_out[class_df.index].to_csv(index=True).encode('utf-8'), f"expression_raw_abundance_{'_'.join(ko_targets)}.csv", "text/csv")
else:
    st.warning("⚠️ No Knockout targets set. Go to **Global Settings** to pick your target proteins.")