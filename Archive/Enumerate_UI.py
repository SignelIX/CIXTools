import json
import pandas as pd
from numpy import random
import MolDisplay
import pathlib
import streamlit as st
import re
import sys
import os
import argparse
import yaml
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import  GridOptionsBuilder
from st_aggrid import GridUpdateMode, DataReturnMode
import numexpr
from streamlit import runtime
import Enumerate
import numpy as np

class EnumerationUI:
    enum = Enumerate.Enumerate()
    initpath =  str(pathlib.Path (__file__).parent.parent)  + '/CIxTools.init.json'
    smiles_colnames = ['SMILES', 'smiles']
    bbid_colnames = ['BB_ID', 'id']
    paramslist = ['enumerate_rxnscheme', 'enumerate_remdups',  'enumerate_lastschemepath', 'enumerate_specstr', 'enumerate_schemename', 'enumerate_rndcount']
    schemelist = []
    schemename = None
    specstr = None
    lspath = None
    rxnschemefile = None
    bbdfs = None
    rxtnts = None
    smilescol='__SMILES'
    ls = None
    Enumerate = False

    SchemeDef_exp = None
    Setup_col = None
    Activity_col = None
    Def_col = None
    Export_exp = None
    grid_exp = None
    Enumeration_exp = None
    ccont = None
    ccontA = None

    def __init__ (self):
        try:
            st.set_page_config(layout="wide")
        except:
            print ('page config already set')
        f = open(self.initpath)
        initjson = json.load(f)
        f.close()
        for p in self.paramslist:
            if p not in st.session_state or st.session_state[p] == None or st.session_state [p] == '' :
                if p in initjson:
                    st.session_state[p] = initjson [p]
                    if p == 'enumerate_rxnscheme':
                        self.Init_RxnSchemeFile()
                else:
                    st.session_state[p] = ''

    def head(self):
        st.markdown("""
            <h1 style='text-align: center; margin-bottom: -35px;'>
            Enumeration Tester
            </h1>
        """, unsafe_allow_html=True
                    )

    def UpdateBBDfs(self, inpath, do_reload):

        flist = pathlib.Path(inpath).glob('*.csv')

        bbfiles =  {}
        dfs = {}

        for f in flist:
            c = str(f)
            result = re.search(r'\.BB[0-9]+\.',  c)
            if result is  not None:
                bbfiles [result.group(0)[3:4]] = c

        if not do_reload:
            for k in bbfiles.keys ():
                if 'df' + str(k) not in st.session_state:
                    do_reload = True
        if do_reload:
            dfs = self.enum.load_BBlists (bbfiles.values (), self.bbid_colnames, self.smiles_colnames)
            for k in dfs.keys ():
                st.session_state['df' + str(k)] = dfs[k]
        else:
            for k in dfs.keys():
                dfs[k] = st.session_state['df' + str(k)]

        dflist = sorted(dfs.items())
        dflist = [x[1] for x in dflist]
        return dflist

    def SaveScheme (self):
        with open(st.session_state['enumerate_rxnscheme'], "r") as jsonFile:
            data = json.load(jsonFile)
        schemejson = json.loads (st.session_state ['enumerate_lastschemedef'])
        data[st.session_state['enumerate_schemename']] = schemejson
        with open(st.session_state['enumerate_rxnscheme'], "w") as jsonFile:
            json.dump(data, jsonFile, indent=4)

    def SetScheme (self):
        st.session_state['rxtnts'] = [None, None, None]
        for k in st.session_state.keys ():
            if k.startswith( 'bb') and k.endswith('idx'):
                st.session_state[k] = 0
            if k.startswith( 'bb') and k.endswith('txt'):
                st.session_state[k] = ''
        if os.path.exists (st.session_state['enumerate_rxnscheme']):
            with open(st.session_state['enumerate_rxnscheme'], "r") as jsonFile:
                data = json.load(jsonFile)
                if st.session_state['enumerate_schemename'] in data:
                    st.session_state ['enumerate_lastschemedef'] = json.dumps (data [st.session_state['enumerate_schemename']], indent=4)
                    self.SaveToInit()
        else:
            st.write ('Rxn Scheme File Not Found')

    def SaveToInit(self):
        with open(self.initpath, "r") as jsonFile:
            try:
                data = json.load(jsonFile)
            except Exception as e :
                st.write(self.initpath)
                st.write (str(e))
                return
        for p in self.paramslist:
            if p in st.session_state:
                data[p] = st.session_state[p]
        with open(self.initpath, "w") as jsonFile:
            jsonFile.write(json.dumps(data, indent=4))

    def Init_RxnSchemeFile (self):
        self.enum.named_reactions = {}
        if not os.path.exists (st.session_state['enumerate_rxnscheme'],):
            st.write ('Rxn Scheme Path Not Found')
            return
        self.enum.named_reactions.update(self.enum.ReadRxnScheme(st.session_state['enumerate_rxnscheme'], 'Named_Reactions', FullInfo=True)[0])
        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        for rxl in data ['rxnlists']:
            if '/' not in rxl:
                rxl = str (pathlib.Path (__file__).parent.parent) + '/' + rxl
            with open(rxl, "r") as jsonFile:
                data = json.load(jsonFile)
            self.enum.named_reactions.update (data)
        self.SaveToInit ()

    def Setup_Column (self):

        if ('schemename' in st.session_state):
            self.ls = st.session_state['schemename']
        else:
            self.ls = ''
        with self.Setup_col:
            with st.expander('Setup', expanded=True):
                self.rxnschemefile = st.text_input(label='Rxn Scheme File', on_change=self.Init_RxnSchemeFile,
                                              key='enumerate_rxnscheme')
                rxnschemerefresh = st.button(label='Re-initialize to Rxn Scheme File')
                if rxnschemerefresh == True:
                    self.Init_RxnSchemeFile()
                if os.path.exists(self.rxnschemefile) == False:
                    st.text(self.rxnschemefile + ' does not exist. Please adjust path')
                    return

                self.lspath = st.text_input(label='Scheme BB Path', on_change=self.SaveToInit, key='enumerate_lastschemepath')
                if self.lspath == '' or self.lspath is None:
                    st.write ('bbpath not specified')
                    return
                else:
                    if not self.lspath.endswith ('/'):
                        self.lspath += '/'
                f = open(self.rxnschemefile, 'r')
                schemejson = json.load(f)

                if 'enumerate_schemename' not in st.session_state:
                    st.session_state['enumerate_schemename'] = self.ls
                    bbpath = self.lspath + '/BBLists'
                    self.specstr = ''
                else:
                    if 'enumerate_specstr' in st.session_state:
                        self.specstr = st.session_state['enumerate_specstr']
                        slash_specstr = '/' + str(self.specstr)
                    else:
                        self.specstr = ''
                        slash_specstr = ''
                    bbpath = self.lspath + str(st.session_state['enumerate_schemename']) + slash_specstr + '/BBLists'

                for s in schemejson:
                    self.schemelist.append(s)
                self.schemelist.sort()

                self.bbdfs = self.UpdateBBDfs(bbpath, False)

                if self.ls == '' or self.ls not in self.schemelist:
                    lsidx = 0
                    if self.schemelist [0] == 'Named_Reactions' and len(self.schemelist) > 1:
                        lsidx = 1
                else:
                    lsidx = self.schemelist.index(self.ls)
                self.schemename = st.selectbox(label='Scheme', options=self.schemelist, key='schemename', index=lsidx)
                self.specstr = st.text_input(key='spec', label='specstr')

    def Activity_Column (self):
        self.SchemeDef1()
        self.Grid_Exp()
        self.Enumerate_Exp()
        self.SchemeDef2()
        self.Export_Expander()
        self.Enumerate_Clicked()

    def SchemeDef1(self):
        with self.SchemeDef_exp:
            st.write ('Test')
            self.ccont = st.container()
            self.ccontA = st.container()
            with self.ccontA:
                newname = st.text_input(label='New Scheme Name')
                if st.button(label='Add New Scheme'):
                    st.session_state['enumerate_schemename'] = newname
                    self.ls = newname
                    st.session_state[
                        'enumerate_lastschemedef'] = '{\n"altnames":[],\n"steps":[\n{"stepname":"",\n"Reactants":["",""],\n"Rxns":{\n"default":[]\n}\n}\n],"scaffold_dummy_structures": [], "filters":{"names":{}, "BBfilters":{}}}'
                    self.SaveScheme()
                    self.SetScheme()
                    self.schemelist.append(newname)
                    self.bbdfs = []



            if 'spec' not in st.session_state and self.specstr is not None and self.specstr != '':
                st.session_state['spec'] = self.specstr


            if self.schemename != st.session_state['enumerate_schemename'] or self.specstr != st.session_state[
                'enumerate_specstr']:
                addspec = ''
                if self.specstr != '' and self.specstr is not None:
                    addspec = '/' + self.specstr
                if self.lspath is not None and self.lspath != '':
                    self.bbdfs = self.UpdateBBDfs(
                        st.session_state['enumerate_lastschemepath'] + self.schemename + addspec + '/BBLists', True)
                st.session_state['enumerate_schemename'] = self.schemename
                st.session_state['enumerate_specstr'] = self.specstr
                st.session_state['aggriddata'] = None
                self.SaveToInit()
                self.SetScheme()
            else:
                if 'enumerate_lastschemedef' not in st.session_state:
                    st.session_state['enumerate_schemename'] = self.schemename
                    st.session_state['enumerate_specstr'] = self.specstr
                    self.SetScheme()
    def SchemeDef2(self):
        with self.SchemeDef_exp:
            ccol1, ccol2 = st.columns(2)
            st.text('Current Scheme: ' + str (st.session_state['enumerate_schemename']))
            with ccol1:
                if st.button(label='revert'):
                    self.SetScheme()
            with ccol2:
                if st.button(label='save scheme'):
                    self.SaveScheme()
                    self.Enumerate = True
            with self.ccont:
                dcol1, dcol2 = st.columns([2, 1])
                with dcol2:
                    if st.button(label='Add Step'):
                        print(st.session_state['enumerate_lastschemedef'])
                        defjson = json.loads(st.session_state['enumerate_lastschemedef'])

                        stepdict = {}
                        stepdict["Reactants"] = []
                        stepdict["Rxns"] = {'default': ''}
                        defjson['steps'].append(stepdict)
                        st.session_state['enumerate_lastschemedef'] = json.dumps(defjson, indent=4)
                with dcol1:
                    st.text_area(height=450, label='Scheme Definition', key='enumerate_lastschemedef')

    def Enumerate_Exp (self):
        self.Enumerate = False

        with self.Enumeration_exp :
            col1, col2, col3 = st.columns(3)
            with col1:
                if st.button('random'):
                    for dx in range(0, len(self.bbdfs)):
                        if self.bbdfs[dx] is not None:
                            st.session_state['bb' + str(dx) + 'idx'] = self.bbdfs[dx].index[
                                random.randint(0, len(self.bbdfs[dx]))]
                    self.Enumerate = True
            with col2:
                if st.button('enumerate'):
                    self.Enumerate = True
            with col3:
                if st.button('clear'):
                    for dx in range(0, len(self.bbdfs)):
                        st.session_state['bb' + str(dx) + 'idx'] = 0
                        st.session_state['bb' + str(dx) + 'txt'] = ''
        if self.bbdfs is not None:
            for n in range(0, len(self.bbdfs)):
                if 'bb' + str(n) + 'idx' not in st.session_state:
                    st.session_state['bb' + str(n) + 'idx'] = 0
                if 'bb' + str(n) + 'txt' not in st.session_state:
                    st.session_state['bb' + str(n) + 'txt'] = ''

            with self.Enumeration_exp :
                self.rxtnts = [None] * len(self.bbdfs)
                for n in range(0, len(self.bbdfs)):
                    df = self.bbdfs[n]
                    self.rxtnts[n] = st.text_input(key='bb' + str(n) + 'txt', label='BB' + str(n + 1) + ' (override)')
                    if df is not None:
                        if self.rxtnts[n] == '':
                            self.rxtnts[n] = st.selectbox(label='BB' + str(n + 1), options=self.bbdfs[n][self.smilescol]
                                                     , key='bb' + str(n),
                                                     index=int(st.session_state['bb' + str(n) + 'idx']))
                if self.rxtnts != st.session_state ['rxtnts'] :
                    st.session_state ['rxtnts'] = self.rxtnts
                    self.Enumerate = True
    def Enumerate_Clicked (self):
        if self.Enumerate == True:
            print ('ENUMERATING')
            with self.Enumeration_exp:
                if self.rxtnts is None or len(self.rxtnts) == 0:
                    st.text('Reactants not correctly loaded')
                else:
                    try:
                        self.Init_RxnSchemeFile ()
                        res, intermeds, rxninfo= self.enum.TestReactionScheme(self.schemename, self.rxtnts, st.session_state ['enumerate_lastschemedef'], True)
                        if res is None or res == 'FAIL' or res.startswith ('FAIL') or res == 'NOSCHEMEFAIL':
                            print ('HERE')
                            st.text ('Reaction Failure: ' + res)
                            with self.EnumSteps_cont:
                                st.image(MolDisplay.ShowMols([''], cols=1, subImgSize=(400, 400),
                                                         outtype='b64_datauri'))
                        else:
                            try:
                                with self.EnumSteps_cont:
                                    st.image(MolDisplay.ShowMols([res], cols=1, subImgSize=(400, 400),
                                                             outtype='b64_datauri'))
                            except Exception as e2:
                                st.text ('EXCEPTION:' + str(e2))
                    except Exception as e:
                        st.text ('Error: bad scheme definition')
                        for r in self.rxtnts:
                            st.text (r)
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        st.text (str(exc_type) + ' ' +  fname + ' ' + str( exc_tb.tb_lineno))
                        st.text (e)
                        intermeds = None

            with self.EnumSteps_cont:
                coly1, coly2 = st.columns (2)
                with coly1:
                    with st.expander(label='Reagents', expanded = True):
                        st.image(MolDisplay.ShowMols(self.rxtnts, cols=1, subImgSize=(400, 200), outtype='b64_datauri'))
                with coly2:
                    with st.expander(label='Reaction Products', expanded = True):
                        if intermeds is not None and len (intermeds) > 0:
                            st.image(MolDisplay.ShowMols(intermeds, cols=1, subImgSize=(400, 200), outtype='b64_datauri'))
                displist = []
                for x in range (0, len(self.rxtnts)):
                    if intermeds is not None:
                        if  len (intermeds) - 1 < x:
                            if len (rxninfo) - 1 < x:
                                displist.append([self.rxtnts[x], None, None])
                            else:
                                displist.append([self.rxtnts[x], None, rxninfo[x]])
                        else:
                            displist.append ([self.rxtnts [x], intermeds[x], rxninfo [x]])
                dispdf = pd.DataFrame (displist, columns = ['Reactant', 'Intermediate','RxnInfo'])
                with st.expander (label='Rxn Info Grid'):
                    MolDisplay.ShowMols_StreamLit_Grid (dispdf,['Reactant', 'Intermediate'], rowheight = 200 )

    def Grid_Exp (self):
        with self.grid_exp:
            getr100 = st.button(label='get random 100')
            if getr100:
                self.Init_RxnSchemeFile()
                with st.spinner('Enumerating'):
                    resdf = self.enum.EnumFromBBFiles(self.schemename , self.specstr, self.specstr, self.lspath, self.specstr,
                                                      100, self.rxnschemefile,
                                                      SMILEScolnames=self.smiles_colnames,
                                                      BBcolnames=self.bbid_colnames,
                                                      rem_dups=False, returndf=True)
                    st.write (len(resdf))
                    cols =  ['full_smiles']
                    for c in resdf.columns [0:len(resdf.columns) -1]:
                        cols.append (c)
                    resdf = resdf.loc [:, cols]
                    st.session_state['aggriddata'] = resdf

            if 'aggriddata' in  st.session_state and st.session_state ['aggriddata'] is not None:
                gb = GridOptionsBuilder.from_dataframe(st.session_state['aggriddata'])

                gb.configure_selection('single')
                gridOptions = gb.build()
                selected = AgGrid(st.session_state['aggriddata'], update_mode='SELECTION_CHANGED',
                                          gridOptions= gridOptions,  data_return_mode=DataReturnMode.AS_INPUT, reload_data=True, custom_css={ "#gridToolBar": { "padding-bottom": "0px !important", } })
                if len(selected['selected_rows']) > 0:
                    for n in range(0, len(self.bbdfs) ):
                        st.session_state['bb' + str(n) + 'txt'] = selected['selected_rows'][0]['bb' + str (n+1) + '_smiles' ]
                    self.Enumerate = True

    def Export_Expander (self):
        with self.Export_exp:
            with st.spinner ('Exporting'):
                expval = st.button('Export Random Selection')
                countval = st.text_input(label='Count', key='enumerate_rndcount', on_change=self.SaveToInit )
                remdup_val = False
                if 'enumerate_remdups' in st.session_state:
                    if st.session_state ['enumerate_remdups'] == 'True':
                        remdup_val = True
                remdups = st.checkbox (label='Remove Duplicate Products', value = remdup_val )
                addintermeds = st.checkbox(label='Add Intermediates', value=remdup_val)
                if expval:
                    try:
                        ct = int (countval)
                    except:
                        ct = 5000
                    self.enum.EnumFromBBFiles(self.schemename, self.specstr, self.specstr, self.lspath,  self.specstr, ct, self.rxnschemefile, SMILEScolnames=self.smiles_colnames, BBcolnames=self.bbid_colnames, rem_dups=remdups, retIntermeds=addintermeds)

    def Layout (self):
        self.Setup_col, self.Activity_col, Output_col = st.columns([2, 2,2])
        self.Setup_Column()
        with self.Setup_col:
            self.SchemeDef_exp = st.expander(label='Scheme Definition Tools', expanded=True)
        with self.Activity_col:
            self.Enumeration_exp = st.expander(label='Enumeration', expanded=True)
            self.grid_exp = self.Scroll_Expanders (title= 'Test structure grid', expanded=True)

        with Output_col:
            EnumStepsTab, ExportTab =  st.tabs (['Single Enumeration', 'Export'])
            with EnumStepsTab:
                self.EnumSteps_cont = st.container()
            with ExportTab:
                self.Export_exp = st.expander(label='Export', expanded=True)
        self.Activity_Column()

    def Scroll_Expanders (self, title, expanded):
        css = '''
              <style>
                  [data-testid="stExpander"] div:has(>.streamlit-expanderContent) {
                      overflow: scroll;
                      height: 400px;
                  }
              </style>
              '''
        st.markdown(css, unsafe_allow_html=True)
        exp = st.expander(label=title, expanded = expanded)



        return exp

    def body(self):
        if 'enumerate_specstr' not in st.session_state:
            st.session_state['enumerate_specstr'] = 'Empty'
        self.Layout()

    def RunUI(self):
        self.head()
        self.body()

if __name__=="__main__":
    if runtime.exists ():
        enum = EnumerationUI ()
        enum.RunUI ()


