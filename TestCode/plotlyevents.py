import streamlit as st
from streamlit_plotly_events import plotly_events
import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from streamlit import runtime


import sys
sys.path.extend(['..'])
import MolDisplay

@st.cache (show_spinner=False)
def read_data ():
    df = pd.read_csv ('~/OneDrive/Consulting/Projects/Triana/22_10_05_Macrocycle/MACROCYCLE_2/Samples/MACROCYCLE_2.full.all.props.csv')
    return np.array(df['ExactMW']), np.array (df['SP3'])

def runplot ():
    df, np1 = read_data ()
    fig = px.scatter(df, x='ExactMW', y='SP3', render_mode='webgl')
    selected_points = plotly_events(fig, click_event=True, hover_event=False , key = 'test3')
    st.text (str (selected_points))

    if len(selected_points) > 0:
        idx = selected_points[0] ['pointIndex']
        smiles = np1[idx,12]
        st.text (smiles)
        st.image (MolDisplay.ShowMols([smiles], cols = 1, subImgSize = (200,200), outtype = 'bytesio' ))

def Run2 ():
    np0, np1 = read_data()

    # Create figure
    fig = go.Figure()

    fig.add_trace(
        go.Scattergl(
            x = np0,
            y = np1,
            mode = 'markers',
            marker = dict(
                line = dict(
                    width = 1,
                    color = 'DarkSlateGrey')
            )
        )
    )

    selected_points = plotly_events(fig, click_event=True, hover_event=False, key='test3')
    st.text(str(selected_points))

    if len(selected_points) > 0:
        idx = selected_points[0]['pointIndex']
        # smiles = np1[idx, 12]
        # st.text(smiles)
        #st.image(MolDisplay.ShowMols([smiles], cols=1, subImgSize=(200, 200), outtype='bytesio'))


if __name__=="__main__":
    if runtime.exists():
        Run2 ()