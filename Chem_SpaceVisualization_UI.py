import Chem_Similarity_Diversity
import pandas as pd
import streamlit as st
import json
import NamingStandardizer
import matplotlib.pyplot as plt
import umap
import pickle
import seaborn as sns
import Enumeration_CombinatorialLibrary as cle
import numpy as np
import matplotlib
import pathlib
import Chem_SpaceVisualization

class Chem_SpaceVisualizationUI:
    initpath = '../CIxTools.init.json'
    paramslist = ['chemspace_outpath', 'chemspace_strucct', 'chemspace_outprefix', 'chemspace_infilelist',
                  'chemspace_pklfile', 'chemspace_geninfilelist', 'chemspace_genoutpath',  'chemspace_genoutprefix',  'chemspace_genfrac']
    ChSV = Chem_SpaceVisualization.Chem_SpaceVisualization ()
    def __init__ (self):
        f = open(self.initpath)
        initjson = json.load(f)
        f.close()
        for p in self.paramslist:
            if p not in st.session_state or st.session_state[p] == None or st.session_state [p] == '' :
                if p in initjson:
                    st.session_state[p] = initjson [p]
                else:
                    st.session_state[p] = ''
    def body (self):
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
            Umap Chemspace Visualization</h1>""", unsafe_allow_html=True)

        with st.expander(label='Generate Map Model'):
            gen_infilelist = st.text_area(label='chemspace_geninfilelist', key='chemspace_geninfilelist', on_change=self.SaveToInit).strip().split('\n')
            gen_outpath = st.text_input(label='chemspace_genoutpath', key='chemspace_genoutpath', on_change=self.SaveToInit)
            gen_outprefix = st.text_input(label='chemspace_genoutprefix', key='chemspace_genoutprefix', on_change=self.SaveToInit)
            gen_frac = st.text_input(label='chemspace_genfrac', key='chemspace_genfrac', on_change=self.SaveToInit)
            if st.button(label='run umap generation', key='RunUMAP'):
                print (float(gen_frac))
                pklfile, imgfile, csvfile, plotlyfig = self.ChSV.CreateMultiFile_UMap(gen_outpath, gen_outprefix, gen_infilelist, frac= float(gen_frac))
                st.text(pklfile)
                st.image(imgfile)
                if plotlyfig is not None:
                    st.plotly_chart(plotlyfig)
        with st.expander(label='Run Visualization', expanded= True):
            outpath = st.text_input (label ='chemspace_outpath', key='chemspace_outpath', on_change=self.SaveToInit)
            outprefix = st.text_input (label='chemspace_outprefix', key='chemspace_outprefix',on_change=self.SaveToInit)
            strucct = st.text_input(label='chemspace_strucct', key='chemspace_strucct',  on_change=self.SaveToInit)
            infilelist = st.text_area (label='chemspace_infilelist', key='chemspace_infilelist', on_change=self.SaveToInit).strip ().split ('\n')
            pkfile = st.text_input(label='chemspace_pklfile', key='chemspace_pklfile', on_change=self.SaveToInit )
            underlying_mapfile = pkfile.replace ('.pkl', '.csv')
            if st.button(label='run umap vs model', key='RunUMAPvModel'):
                with st.spinner('Building Map...'):
                    fig_file, plotlyfig = self.ChSV.ChemicalSpace_UMap(outpath, outprefix, infilelist, ['red', 'blue', 'aqua', 'purple'],
                                                       pkfile, int(strucct), underlying_mapfile)
                    st.image(fig_file)
                    if plotlyfig is not None:
                        st.plotly_chart(plotlyfig)
                    st.success('Done!')

    def SaveToInit(self):

        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        for p in self.paramslist:
            if p in st.session_state:
                data[p] = st.session_state[p]
        with open(self.initpath, "w") as jsonFile:
            print (json.dumps(data, indent=4))
            jsonFile.write(json.dumps(data, indent=4))

    def RunUI (self):
        self.body ()

if __name__ == "__main__":
    csvis = Chem_SpaceVisualization ()
    csvis.RunUI()