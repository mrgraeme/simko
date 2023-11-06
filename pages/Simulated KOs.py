import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns

st.set_page_config(
    page_title="SimKO - SimKO",
    page_icon="🧫",
)

@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('~/work/data/global/ccle/proteomics/CCLE_imputed_all.csv')
    abundance = abundance.rename(columns={'Gene':'Protein'})
    return abundance

def get_classes_by_mean_abundance(protein_list, abundance, n=10):
    protein_list_abundance = abundance.loc[abundance['Protein'].isin(protein_list)]
    protein_list_abundance = protein_list_abundance.T
    protein_list_abundance.columns = protein_list_abundance.iloc[0]
    protein_list_abundance = protein_list_abundance.iloc[1:]
    protein_list_abundance['mean'] = protein_list_abundance.mean(axis=1)
    protein_list_abundance = protein_list_abundance.sort_values("mean", ascending = False)
    high = protein_list_abundance.head(n)
    high['class'] = 'high'
    low = protein_list_abundance.tail(n)
    low['class'] = 'low'
    protein_classes = pd.concat([high, low])
    return pd.DataFrame(protein_classes)

def get_differentials(class_df, abundance):
    abundance = abundance.set_index('Protein', append=True)
    high_class = list(class_df.loc[class_df['class']=='high'].index)
    low_class = list(class_df.loc[class_df['class']=='low'].index)
    abundance_diff = pd.DataFrame()
    abundance_diff['high'] = abundance[high_class].mean(axis=1)
    abundance_diff['low'] = abundance[low_class].mean(axis=1)
    abundance_diff = abundance_diff.reset_index().iloc[:,1:]
    abundance_diff['diff'] = (abundance_diff['low'] - abundance_diff['high'])
    abundance_diff = abundance_diff.sort_values('diff', ascending=True)
    return abundance_diff 



# protein_list = ['ARID1A', 'PBRM1', 'BRAF']
# get_classes_by_mean_abundance(protein_list, abundance)


abundance = get_abundance_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

# x = st.slider('x', min_value=10, max_value=50)  # 👈 this is a widget
# st.write(x, 'squared is', x * x)
# y=st.slider('y', min_value=0.1, max_value=1.0)  # 👈 this is a widget

# chart_data = chart_data[0:x]
# chart_data = chart_data * y

protein_list = st.multiselect(
    'Proteins for KO',
     abundance['Protein'], placeholder='Add proteins to analyse')


cmap = plt.cm.get_cmap('RdYlBu_r')
if protein_list:
    class_df = get_classes_by_mean_abundance(protein_list, abundance)
    
    tab1, tab2 = st.tabs(["figure", "data"])

    with tab1:
        fig = plt.figure(figsize=(10, 4))
        sns.heatmap(class_df[['mean']].sort_values('mean').T.astype(float).round(2), square=True, cmap="vlag", annot=True,annot_kws={'size': 6}, cbar=False)
        st.pyplot(fig)

    with tab2:
        st.table(class_df.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None))


  
    
    n = st.slider('Number of top / bottom proteins to show by differential', value = 5, min_value=0, max_value=50, )  # 👈 this is a widget
    protein_select = st.multiselect(
    'Select additional proteins differences to view',
     abundance['Protein'], placeholder='Add additional proteins to view')



    diff_df = get_differentials(class_df, abundance)
    diff_proteins = list(diff_df.head(n)['Protein']) + list(diff_df.tail(n)['Protein'])

    if protein_select:
        diff_proteins = diff_proteins + protein_select

    diff_df = diff_df.loc[diff_df['Protein'].isin(diff_proteins)]
    st.table(diff_df.style.background_gradient(cmap = cmap, vmin=(-6), vmax=6, axis=None))


    


else:
    st.write('Please select at least one protein.')


