import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(
    page_title="SimKO - Mutation Effect",
    page_icon="🥼",
    layout="wide"
)
@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('./data/abundance.csv')
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


def get_classes_by_mutation(protein_list, mutation, n):
    mutation_filter = mutation.loc[mutation.index.isin(protein_list)].T
    mutation_filter['mutation_count'] = mutation_filter.sum(axis=1)
    mutated = mutation_filter.loc[mutation_filter['mutation_count'] > 0]
    non_mutated = mutation_filter.loc[mutation_filter['mutation_count'] == 0]
    if len(mutated) < len(non_mutated):
        non_mutated = non_mutated.sample(n=len(mutated))
    elif len(mutated) > len(non_mutated):
        mutated = mutated.sample(n=len(non_mutated))
    mutated['class'] = 'mutated'
    non_mutated['class'] = 'non-mutated'
    protein_classes = pd.concat([mutated, non_mutated])
    return pd.DataFrame(protein_classes)

def ttest_from_sample_stats(row, n_cls = 20):
    pooled_sd = np.sqrt((((n_cls-1)*(row['mutated_std']**2)) + ((n_cls-1)*(row['non-mutated_std']**2))) / (n_cls + n_cls - 2))
    # print(pooled_sd)
    t_stat = (row['mutated'] - row['non-mutated']) / (pooled_sd * np.sqrt(1/n_cls + 1/n_cls))
    # print(t_stat)
    p_value = 2*(1 - t.cdf(abs(t_stat), (n_cls + n_cls - 2)))
    # print(p_value)
    return p_value

def get_differentials(class_df, data_df, n):
    non_mutated_class = list(class_df.loc[class_df['class']=='non-mutated'].index)
    mutated_class = list(class_df.loc[class_df['class']=='mutated'].index)
    diff_df = pd.DataFrame()
    diff_df['non-mutated'] = data_df.filter(non_mutated_class).mean(axis=1)
    diff_df['non-mutated_std'] = data_df.filter(non_mutated_class).std(axis=1)
    diff_df['mutated'] = data_df.filter(mutated_class).mean(axis=1)
    diff_df['mutated_std'] = data_df.filter(mutated_class).std(axis=1)
    # diff_df['diff'] = (diff_df['low'] - diff_df['median'])
    diff_df['diff'] = diff_df['mutated'] - diff_df['non-mutated']
    diff_df['p'] = diff_df.apply(ttest_from_sample_stats, n_cls=n, axis=1)
    diff_df = diff_df.drop(columns=['mutated_std', 'non-mutated_std'])
    diff_df = diff_df.sort_values('diff', ascending=True)
    return diff_df

def get_differentials_boxplot(class_df, data_df, protein_list, n):
    # median_class = list(class_df.loc[class_df['class']=='median'].index)
    # low_class = list(class_df.loc[class_df['class']=='low'].index)
    data_df = data_df.reset_index()
    data_df = data_df.loc[data_df['protein'].isin(protein_list)]
    data_df = data_df.melt(id_vars='protein')
    class_df = class_df[['class']].reset_index()
    box_data = data_df.merge(class_df, how='left', left_on='variable', right_on='index').dropna()
    fig = plt.figure(figsize=(10, 4))
    sns.boxplot(data=box_data, x="protein", y="value", hue='class')
    return fig


def process_mutations(class_df, data_df):
    median_class = list(class_df.loc[class_df['class']=='median'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    diff_df = pd.DataFrame()
    diff_df['median'] = data_df.filter(median_class).sum(axis=1)
    diff_df['low'] = data_df.filter(low_class).sum(axis=1)
    diff_df['diff'] = (diff_df['low'] - diff_df['median'])
    diff_df = diff_df.sort_values('diff', ascending=True)
    return diff_df

def get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, protein_list):
    abundance_summary = diff_abund_df.loc[diff_abund_df.index.isin(protein_list)]
    abundance_summary.columns = ['Abundance - ' + n for n in abundance_summary.columns]
    expression_summary = diff_exp_df.loc[diff_exp_df.index.isin(protein_list)]
    expression_summary.columns = ['Expression - ' + n for n in expression_summary.columns]
    mutation_summary = diff_mut_df.loc[diff_mut_df.index.isin(protein_list)]
    mutation_summary.columns = ['Mutation - ' + n for n in mutation_summary.columns]
    diff_summary = pd.concat([abundance_summary, expression_summary, mutation_summary], axis=1)
    return diff_summary


abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

cls = list(abundance.columns)
tissues = set([str(i.split('_', 1)[1:][0]) for i in abundance.columns])

st.write("### Simulate protein KO data! 🥼")

protein_list = st.multiselect(
    'Proteins for KO',
     abundance.index, placeholder='Add proteins to analyse')

tissue_list = st.multiselect(
    'Filter for Tissue',
     tissues, placeholder='Specify specific tissues')

if tissue_list:
    cls = [s for s in cls if any(xs in s for xs in tissue_list)]
    st.write('Filtered to', len(cls), 'cell lines')
    abundance = abundance.filter(cls)
    expression = expression.filter(cls)
    mutation = mutation.filter(cls)

    # st.dataframe(abundance.head(4).style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None).format("{:.3f}"),)  ### For debugging - delete

if protein_list:
    n = ((abundance.shape[1]-1) // 3) if abundance.shape[1] < 60 else 20
    class_df = get_classes_by_mutation(protein_list, mutation, n)
    
    tab1, tab2 = st.tabs(["figure", "data"])

    with tab1:
        fig = plt.figure(figsize=(10, 4))
        sns.heatmap(class_df[['mutation_count']].sort_values('mutation_count').T.astype(float).round(1), square=True, cmap="vlag", annot=True,annot_kws={'size': 5.5}, cbar=False)
        st.pyplot(fig)

    with tab2:
        st.table(class_df.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None))



    st.markdown(
        """
        **Mean abundance, expression and mutation count change across median and low classes** 
        """
    )


    protein_select = st.multiselect(
        'Select additional proteins differences to view',
        abundance.index, placeholder='Add additional proteins to view')

    if protein_select:
        show_proteins = protein_list + protein_select
    else:
        show_proteins = protein_list

    diff_abund_df = get_differentials(class_df, abundance, n)
    diff_exp_df = get_differentials(class_df, expression, n)
    diff_mut_df = process_mutations(class_df, mutation)

    diff_summary = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, show_proteins)

    fig = get_differentials_boxplot(class_df, abundance, show_proteins, n)
    st.pyplot(fig)


    st.dataframe(diff_summary.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None).format("{:.3f}"),)

    


    st.markdown(
        """
        **Proteins most effected by simulated KO of given proteins** 
        """
    )
    
    n_outlayers = st.slider('Number proteins by median differential', value = 5, min_value=0, max_value=50, )  # 👈 this is a widget
 
    diff_proteins = list(diff_abund_df.loc[diff_abund_df['p'] < 0.01].head(n_outlayers).index) + list(diff_abund_df.loc[diff_abund_df['p'] < 0.01].tail(n_outlayers).index)

    diff_top = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_proteins)
    st.dataframe(diff_top.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None).format("{:.3f}"),)

    st.download_button('Download All Abundance', 
                       diff_abund_df.to_csv(index=True).encode('utf-8'), 
                       "abundance_%s.csv" % ('_'.join(protein_list)),
                       "text/csv",
                       'download-csv')
    


else:
    st.write('Please select at least one protein.')


