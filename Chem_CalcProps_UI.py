import Chem_CalcProps
import streamlit as st
import json

class Chem_CalcPropsUI:
    initpath = '../CIxTools.init.json'
    paramslist = ['props_smilescol', 'props_bbidcols', 'props_infile', 'props_filepath']

    def __init__(self):
        f = open(self.initpath)
        initjson = json.load(f)
        f.close()
        for p in self.paramslist:
            if p not in st.session_state or st.session_state[p] == None or st.session_state[p] == '':
                if p in initjson:
                    st.session_state[p] = initjson[p]
                else:
                    st.session_state[p] = ''

    def head(self):
        st.markdown("""
            <h1 style='text-align: center; margin-bottom: -35px;'>
            Calculate Properties
            </h1>
        """, unsafe_allow_html=True
                   )
        return

    def body (self):
        filename = st.text_input (label = 'Props Input File', key = 'props_infile',  on_change=self.SaveToInit)
        smiles_col = st.text_input(label='SMILES Column', key='props_smilescol',  on_change=self.SaveToInit)
        bbidcols = st.text_input(label='BBID Cols', key='props_bbidcols',  on_change=self.SaveToInit).split(',')
        with st.spinner('Calculating...'):
            RunClusters = st.button(label='Run', key='RunProps')
            if RunClusters == True:
                outfile = filename.replace ('.csv', '.props.csv')
                if outfile == filename:
                    st.text ('Incompatible input file, .csv file is required')
                    return
                Chem_CalcProps.props.addPropsToFile(filename, outfile,  smiles_col, bbidcols)
        with st.expander ('Property Graphs'):
            infile = st.text_input ('property file path', key='props_filepath', on_change=self.SaveToInit)
            gen_graphs = st.button ('Generate Graphs')
            if gen_graphs:
                with st.spinner ('Generating Plots...'):
                    res  = Chem_CalcProps.props.GeneratePlots(None, infile )
                    print (res)
                    st.image (res, width = 600)
        return

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
        self.head ()
        self.body ()
