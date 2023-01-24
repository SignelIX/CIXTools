import pandas as pd
import streamlit as st
import json
import sys
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import  GridOptionsBuilder
from streamlit import runtime
import Chem_Similarity_Diversity

class Chem_Similarity_DiversityUI:
    initpath = '../CIxTools.init.json'
    paramslist = ['simdiv_inputfile', 'simdiv_inputsmiles', 'bbsim_inputfile', 'bbsim_inputsmiles', 'bbsim_rxnschemefile','bbsim_schemename',
                  'bbsim_rxtntnum', 'bbsim_retct', 'bbsim_filtcol', 'bbsim_filtvals', 'bbsim_addtlcols', 'bbsim_savetofldr']
    sim = Chem_Similarity_Diversity.Similarity ()

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
    def body (self, param):
        if param == 'Sim Search':
            self.simsearchbody ()
        if param == 'BB Similarity':
            self.bbsimbody ()

    def simsearchbody (self):
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
            Similarity Search</h1>""", unsafe_allow_html=True)
        inputfile = st.text_input (label = 'input file', key = 'simdiv_inputfile', on_change=self.SaveToInit)
        smiles = st.text_input (label = 'input smiles', key = 'simdiv_inputsmiles', on_change=self.SaveToInit)
        sep = st.selectbox (label= 'Delimiter' , options= ['COMMA', 'TAB'], index = 0)
        if st.button (label = 'Run Similarity (smiles v. file)'):
            with st.spinner('Running Similarity'):
                outpath = inputfile.replace('.csv', '') + '.simrun.csv'
                if sep == 'TAB':
                    sep = '\t'
                else:
                    sep  =','
                self.sim.SimSearchLargeFile([smiles], inputfile, outpath, sep= sep)

    def bbsimbody(self):
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
               BB Similarity Search</h1>""", unsafe_allow_html=True)
        inputfile = st.text_input(label='input file', key='bbsim_inputfile', on_change=self.SaveToInit)
        smiles = st.text_input(label='input smiles (; separated', key='bbsim_inputsmiles', on_change=self.SaveToInit)
        rxnschemefile = st.text_input(label='rxnschemefile', key='bbsim_rxnschemefile', on_change=self.SaveToInit)
        schemename = st.text_input(label='schemename', key='bbsim_schemename', on_change=self.SaveToInit)
        rxtntnum = st.text_input(label='rxtntnum (0 based)', key='bbsim_rxtntnum', on_change=self.SaveToInit)
        retct = st.text_input(label='retct', key='bbsim_retct', on_change=self.SaveToInit)
        sep = st.selectbox(label='Delimiter', options=['COMMA', 'TAB'], index=0)
        filtercol = st.text_input(label='filtcol', key='bbsim_filtcol', on_change=self.SaveToInit)
        filtvals = st.text_input(label='filtvals (; separated)', key='bbsim_filtvals', on_change=self.SaveToInit)
        addtlcols = st.text_input(label='addtlcols(; separated)', key='bbsim_addtlcols', on_change=self.SaveToInit)
        savetofldr = st.text_input(label='savetofldr', key='bbsim_savetofldr', on_change=self.SaveToInit)

        if st.button(label='Run BB Similarity (smiles v. file)'):
            with st.spinner('Running BB Similarity'):
                outpath = inputfile.replace('.csv', '') + '.simrun.csv'
                if sep == 'TAB':
                    sep = '\t'
                else:
                    sep = ','
                comp_bbdf = pd.read_csv (inputfile)
                filtvalsplit = filtvals.split(';')
                addtlcolssplit = addtlcols.split(';')
                comp_bbdf=comp_bbdf[comp_bbdf[filtercol].isin( filtvalsplit)]
                smilessplit = smiles.split (';')
                outdf = None
                for s in smilessplit:
                    # if outdf is not None:
                    #     excludelist = list(outdf['BB_ID'])
                    # else:
                    excludelist = None

                    res = self.sim.Find_MostSimilarBB(s, comp_bbdf, rxnschemefile, schemename, int(rxtntnum), int(retct), excludelist=excludelist, addtlcols= addtlcolssplit)
                    cols = ['BB_ID', 'SIMILARITY',  'SMILES']

                    for c in addtlcolssplit:
                        cols.append(c)
                    print (cols)
                    resdf = pd.DataFrame (res, columns = cols)
                    resdf ['search_smiles'] = s
                    if outdf is None:
                        outdf = resdf
                    else:
                        outdf = outdf.append(resdf)


                st.session_state['aggriddata'] = outdf

                gb = GridOptionsBuilder()
                gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, editable=True,enableRangeSelection=True, )
                gb.configure_column("BB_ID", headerName="BB ID", width=100)
                gb.configure_column("SIMILARITY", headerName="Similarity", type=["numericColumn", "numberColumnFilter"],width=50)
                gb.configure_column("SMILES", headerName="SMILES", width=50)
                gb.configure_column("search_smiles", headerName="search_smiles", width=50)
                for c in addtlcolssplit:
                    gb.configure_column(c, headerName=c, width=50)
                gridOptions = gb.build()
                st.session_state['gridOptions'] = gridOptions
                AgGrid(st.session_state['aggriddata'], enable_enterprise_modules=True,gridOptions=gridOptions )
                if savetofldr is not None:
                    outdf.to_csv(savetofldr + '/' + schemename + '_' + str(rxtntnum) + '.csv', index = False)
        else:
            if 'aggriddata' in  st.session_state:
                AgGrid(st.session_state['aggriddata'], enable_enterprise_modules=True,gridOptions=st.session_state['gridOptions'],)


    def SaveToInit(self):

        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        for p in self.paramslist:
            if p in st.session_state:
                data[p] = st.session_state[p]
        with open(self.initpath, "w") as jsonFile:
            jsonFile.write(json.dumps(data, indent=4))

    def RunUI (self, param):
        self.body (param)

if __name__ == "__main__":
    if runtime.exists():
        if len (sys.argv) == 0:
            csvis = Chem_Similarity_DiversityUI ('Sim Search')
            csvis.RunUI()

