import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(
    page_title="SimKO - Simulate KOs",
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


def get_classes_by_mean_abundance(protein_list, abundance, n=20):
    protein_list_abundance = abundance.loc[abundance.index.isin(protein_list)]
    protein_list_abundance = protein_list_abundance.T
    # protein_list_abundance.columns = protein_list_abundance.iloc[0]
    # protein_list_abundance = protein_list_abundance.iloc[1:]
    protein_list_abundance['mean'] = protein_list_abundance.mean(axis=1)
    protein_list_abundance = protein_list_abundance.sort_values("mean", ascending = False)
    high = protein_list_abundance.head(n)
    high['class'] = 'high'
    low = protein_list_abundance.tail(n)
    low['class'] = 'low'
    protein_classes = pd.concat([high, low])
    return pd.DataFrame(protein_classes)

def get_differentials(class_df, data_df):
    high_class = list(class_df.loc[class_df['class']=='high'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    diff_df = pd.DataFrame()
    diff_df['high'] = data_df.filter(high_class).mean(axis=1)
    diff_df['low'] = data_df.filter(low_class).mean(axis=1)
    diff_df['diff'] = (diff_df['low'] - diff_df['high'])
    diff_df = diff_df.sort_values('diff', ascending=True)
    return diff_df

def process_mutations(class_df, data_df):
    high_class = list(class_df.loc[class_df['class']=='high'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    diff_df = pd.DataFrame()
    diff_df['high'] = data_df.filter(high_class).sum(axis=1)
    diff_df['low'] = data_df.filter(low_class).sum(axis=1)
    diff_df['diff'] = (diff_df['low'] - diff_df['high'])
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

# protein_list = ['ARID1A', 'PBRM1', 'BRAF']
# get_classes_by_mean_abundance(protein_list, abundance)



# x = st.slider('x', min_value=10, max_value=50)  # 👈 this is a widget
# st.write(x, 'squared is', x * x)
# y=st.slider('y', min_value=0.1, max_value=1.0)  # 👈 this is a widget

# chart_data = chart_data[0:x]
# chart_data = chart_data * y

st.write("### Simulate protein KO data! 🥼")


protein_list = st.multiselect(
    'Proteins for KO',
     abundance.index, placeholder='Add proteins to analyse')


if protein_list:
    class_df = get_classes_by_mean_abundance(protein_list, abundance)
    
    tab1, tab2 = st.tabs(["figure", "data"])

    with tab1:
        fig = plt.figure(figsize=(10, 4))
        sns.heatmap(class_df[['mean']].sort_values('mean').T.astype(float).round(1), square=True, cmap="vlag", annot=True,annot_kws={'size': 5.5}, cbar=False)
        st.pyplot(fig)

    with tab2:
        st.table(class_df.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None))

    st.markdown(
        """
        **Mean abundance, expression and mutation count change across high and low classes** 
        """
    )


    protein_select = st.multiselect(
        'Select additional proteins differences to view',
        abundance.index, placeholder='Add additional proteins to view')

    if protein_select:
        show_proteins = protein_list + protein_select
    else:
        show_proteins = protein_list

    diff_abund_df = get_differentials(class_df, abundance)
    diff_exp_df = get_differentials(class_df, expression)
    diff_mut_df = process_mutations(class_df, mutation)

    diff_summary = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, show_proteins)
    st.dataframe(diff_summary.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None).format("{:.3f}"),)


    st.markdown(
        """
        **Proteins most effected by simulated KO of given proteins** 
        """
    )
    
    n = st.slider('Number proteins by high differential', value = 5, min_value=0, max_value=50, )  # 👈 this is a widget
 
    
    diff_proteins = list(diff_abund_df.head(n).index) + list(diff_abund_df.tail(n).index)

    diff_top = get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, diff_proteins)
    st.dataframe(diff_top.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None).format("{:.3f}"),)


    


else:
    st.write('Please select at least one protein.')

