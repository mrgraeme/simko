import streamlit as st
import pandas as pd
import numpy as np

name = st.text_input("Your name", key="name")

# You can access the value at any point with:
st.session_state.name




@st.cache_data
def get_data():
    chart_data = pd.DataFrame(
        np.random.randn(50, 3),
        columns=['a', 'b', 'c'])
    return chart_data


chart_data = get_data()

x = st.slider('x', min_value=10, max_value=50)  # 👈 this is a widget
st.write(x, 'squared is', x * x)
y=st.slider('y', min_value=0.1, max_value=1.0)  # 👈 this is a widget

chart_data = chart_data[0:x]
chart_data = chart_data * y


option = st.multiselect(
    'Which column do you want to see?',
     chart_data.columns)

if option:
    chart_data = chart_data[option]

st.line_chart(chart_data )


st.dataframe(chart_data.style.highlight_max(axis=0))


map_data = pd.DataFrame(
    np.random.randn(100, 2) / [2, 2] + [51.5, -0.2],
    columns=['lat', 'lon'])

st.map(map_data)

