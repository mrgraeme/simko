import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import streamlit as st


# import simko.simko as simko
# import simko.simko_stream as simko_stream

st.set_page_config(
    layout='wide',
    page_title="IPQC",
    page_icon="🥼",
)

# st.title("Welcome to IPQC!")



pages = {
    "Filter Data": [
        st.Page("pages/Filter_Data.py", title="Filter Data", icon="🔬"),
    ],
    "Effect Analysis": [
        st.Page("pages/Simulated_KO_Effect.py", title="Abundance Effects", icon="📤"),
        st.Page("pages/Expression_Effect.py", title="Expression Effect", icon="🫶"),
        st.Page("pages/Mutation_Effect.py", title="Mutation Effect", icon="👾"),
    ],
    # "Analysis": [
    #     st.Page(
    #         "pages/Functional_Categories.py",
    #         title="Functional Categories",
    #         icon="📊"
    #     ),
    # ],
}

# Create the navigation menu
pg = st.navigation(pages)
pg.run()

# st.write("# Welcome to SimKO! 🥼")

# # st.sidebar.success("Select a demo above.")

# st.markdown(
#     """
#     SimKO is a tool to interrogate and visualise protein abundance change in celllines.
    
#     **👈 Select a tool from the sidebar** 
#     ### Want to learn more?
#     - Check out [pkg.rcrds.com](https://pkg.rcrds.com)
# """
# )
