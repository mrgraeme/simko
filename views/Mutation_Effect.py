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

@st.cache_data
def get_dependency_data():
    return pd.read_csv('./data/dependency.csv').set_index('protein')

# Math & Analytics Utilities
def get_classes_by_mutation(protein_list, mutation, n):
    mutation_filter = mutation.loc[mutation.index.isin(protein_list)].T
    mutation_filter['mutation_count'] = mutation_filter.sum(axis=1)
    mutated = mutation_filter.loc[mutation_filter['mutation_count'] > 0].sort_values('mutation_count', ascending=False).copy()
    non_mutated = mutation_filter.loc[mutation_filter['mutation_count'] == 0].copy()
    
    mutated['class'] = 'mutated'
    non_mutated['class'] = 'non-mutated'
    return pd.concat([mutated, non_mutated])

def ttest_from_sample_stats(row, n_cls=20):
    pooled_sd = np.sqrt((((n_cls-1)*(row['mutated_std']**2)) + ((n_cls-1)*(row['non-mutated_std']**2))) / (n_cls + n_cls - 2))
    if pooled_sd == 0:
        return 1.0
    t_stat = (row['mutated'] - row['non-mutated']) / (pooled_sd * np.sqrt(1/n_cls + 1/n_cls))
    return 2 * (1 - t.cdf(abs(t_stat), (n_cls + n_cls - 2)))

def get_differentials(class_df, data_df, n):
    non_mutated_class = list(class_df.loc[class_df['class']=='non-mutated'].index)
    mutated_class = list(class_df.loc[class_df['class']=='mutated'].index)
    
    diff_df = pd.DataFrame()
    diff_df['non-mutated'] = data_df.filter(non_mutated_class).mean(axis=1)
    diff_df['non-mutated_std'] = data_df.filter(non_mutated_class).std(axis=1)
    diff_df['mutated'] = data_df.filter(mutated_class).mean(axis=1)
    diff_df['mutated_std'] = data_df.filter(mutated_class).std(axis=1)
    
    diff_df['diff'] = diff_df['mutated'] - diff_df['non-mutated']
    diff_df['p'] = diff_df.apply(ttest_from_sample_stats, n_cls=n, axis=1)
    return diff_df.drop(columns=['mutated_std', 'non-mutated_std']).sort_values('diff', ascending=True)

def get_differentials_boxplot(class_df, data_df, protein_list, n):
    data_df = data_df.reset_index()
    data_df = data_df.loc[data_df['protein'].isin(protein_list)]
    data_df = data_df.melt(id_vars='protein')
    class_df = class_df[['class']].reset_index()
    box_data = data_df.merge(class_df, how='left', left_on='variable', right_on='index').dropna()
    
    fig = plt.figure(figsize=(10, 4))
    sns.boxplot(data=box_data, x="protein", y="value", hue='class', palette='pastel')
    plt.xticks(rotation=45)
    plt.tight_layout()
    return fig

def process_mutations(class_df, data_df):
    # Mapping back safely using mutated/non-mutated groups
    m_class = list(class_df.loc[class_df['class']=='mutated'].index)
    nm_class = list(class_df.loc[class_df['class']=='non-mutated'].index)
    diff_df = pd.DataFrame()
    diff_df['non-mutated'] = data_df.filter(nm_class).sum(axis=1)
    diff_df['mutated'] = data_df.filter(m_class).sum(axis=1)
    diff_df['diff'] = diff_df['mutated'] - diff_df['non-mutated']
    return diff_df.sort_values('diff', ascending=True)

def get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_dependency_df, protein_list):
    ab_sub = diff_abund_df.loc[diff_abund_df.index.isin(protein_list)].copy()
    ab_sub.columns = ['Abundance - ' + n for n in ab_sub.columns]
    
    ex_sub = diff_exp_df.loc[diff_exp_df.index.isin(protein_list)].copy()
    ex_sub.columns = ['Expression - ' + n for n in ex_sub.columns]
    
    mu_sub = diff_mut_df.loc[diff_mut_df.index.isin(protein_list)].copy()
    mu_sub.columns = ['Mutation - ' + n for n in mu_sub.columns]
    
    dep_sub = diff_dependency_df.loc[diff_dependency_df.index.isin(protein_list)].copy()
    dep_sub.columns = ['Dependency - ' + n for n in dep_sub.columns]
    
    return pd.concat([ab_sub, ex_sub, mu_sub, dep_sub], axis=1)

# Pipeline Setup
abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
dependency = get_dependency_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

st.write("### Explore mutation effects 👾")

# 1. Pull settings from configuration state
ko_targets = st.session_state.get('saved_ko_proteins', [])
focus_proteins = st.session_state.get('saved_protein_list', [])
active_tissues = st.session_state.get('saved_tissues', [])
active_cell_lines = st.session_state.get('saved_cell_lines', [])

# 2. Re-index active cohort slice
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
dependency = dependency.filter(active_cls)

if ko_targets:
    st.caption(f"Analyzing variants for targets {ko_targets} across {len(active_cls)} matched cell lines.")
    n = ((abundance.shape[1]-1) // 3) if abundance.shape[1] < 60 else 20
    class_df = get_classes_by_mutation(ko_targets, mutation, n)

    if len(class_df) == 0:
        st.error("There is not enough mutation data available to slice cohorts for this target selection.")
    else:
        tab1, tab2 = st.tabs(["Mutation Heatmap", "Cohort Matrix"])
        with tab1:
            fig = plt.figure(figsize=(10, 3))
            sns.heatmap(class_df[['mutation_count']].sort_values('mutation_count').T.astype(float).round(1), square=True, cmap="vlag", annot=True, annot_kws={'size': 6}, cbar=False)
            st.pyplot(fig)
        with tab2:
            st.table(class_df.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None))

        st.markdown("**Mean abundance, expression, and mutation changes across cohorts**")

        protein_select = st.multiselect(
            'Select additional proteins differences to view',
            options=abundance.index, placeholder='Supplement tracking pool'
        )

        show_proteins = list(set(ko_targets + focus_proteins + protein_select))

        diff_abund_df = get_differentials(class_df, abundance, n)
        diff_exp_df = get_differentials(class_df, expression, n)
        diff_dependency_df = get_differentials(class_df, dependency, n)
        diff_mut_df = process_mutations(class_df, mutation)

        if show_proteins:
            diff_summary = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_dependency_df, show_proteins)
            fig_box = get_differentials_boxplot(class_df, abundance, show_proteins, n)
            st.pyplot(fig_box)
            st.dataframe(diff_summary.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))

        st.markdown("**Proteins most significantly affected by mutation state**")
        n_outlayers = st.slider('Number proteins by median differential', value=5, min_value=0, max_value=50)
    
        significant_hits = diff_abund_df[diff_abund_df['p'] < 0.05]
        diff_proteins = list(significant_hits.head(n_outlayers).index) + list(significant_hits.tail(n_outlayers).index)

        diff_top = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_dependency_df, diff_proteins)
        st.dataframe(diff_top.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None).format("{:.3f}"))

        st.download_button('Download Differential Summary Table', diff_abund_df.to_csv(index=True).encode('utf-8'), f"mutation_differentials_{'_'.join(ko_targets)}.csv", "text/csv")
else:
    st.warning("⚠️ No Knockout targets set. Go to **Global Settings** to pick your targets.")