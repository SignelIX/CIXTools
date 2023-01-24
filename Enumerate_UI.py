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

class EnumerationUI:
    enum = Enumerate.Enumerate()
    initpath = '../CIxTools.init.json'
    smiles_colnames = None
    bbid_colnames = None
    paramslist = ['enumerate_rxnscheme', 'enumerate_remdups',  'enumerate_lastschemepath', 'enumerate_specstr', 'enumerate_schemename', 'enumerate_rndcount']

    def __init__ (self):
        f = open(self.initpath)
        initjson = json.load(f)
        f.close()
        for p in self.paramslist:
            if p not in st.session_state or st.session_state[p] == None or st.session_state [p] == '' :
                if p in initjson:
                    st.session_state[p] = initjson [p]
                    if p == 'enumerate_rxnscheme':
                        self.enum.named_reactions = self.enum.ReadRxnScheme(st.session_state['enumerate_rxnscheme'],'Named_Reactions', FullInfo=True)[0]
                        print ('NAMED Rxns:', self.enum.named_reactions )
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
        for k in st.session_state.keys ():
            if k.startswith( 'bb') and k.endswith('idx'):
                st.session_state[k] = 0
        with open(st.session_state['enumerate_rxnscheme'], "r") as jsonFile:
            data = json.load(jsonFile)
            if st.session_state['enumerate_schemename'] in data:
                st.session_state ['enumerate_lastschemedef'] = json.dumps (data [st.session_state['enumerate_schemename']], indent=4)
                self.SaveToInit()

    def SaveToInit(self):
        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        for p in self.paramslist:
            if p in st.session_state:
                data[p] = st.session_state[p]
        with open(self.initpath, "w") as jsonFile:
            jsonFile.write(json.dumps(data, indent=4))

    def Init_RxnSchemeFile (self):
        self.enum.named_reactions= self.enum.ReadRxnScheme(st.session_state ['enumerate_rxnscheme'], 'Named_Reactions', FullInfo = True )[0]
        self.SaveToInit ()

    def body(self):
        smilescol='SMILES'

        if 'enumerate_specstr' not in st.session_state:
            st.session_state['enumerate_specstr'] = 'Empty'

        rxnschemefile = st.text_input(  label='Rxn Scheme File', on_change=self.Init_RxnSchemeFile, key='enumerate_rxnscheme')
        rxnschemerefresh = st.button (label = 'Re-initialize to Rxn Scheme File')
        if rxnschemerefresh == True:
            self.Init_RxnSchemeFile()
        if os.path.exists(rxnschemefile) == False:
            st.text (rxnschemefile + ' does not exist. Please adjust path')
            return

        lspath = st.text_input( label='Scheme Path', on_change=self.SaveToInit, key = 'enumerate_lastschemepath')

        if ('schemename' in st.session_state):
            ls = st.session_state['schemename']
        else:
            ls = ''

        f = open(rxnschemefile, 'r')
        schemejson = json.load(f)
        schemelist = []

        if 'enumerate_schemename' not in st.session_state:
            st.session_state['enumerate_schemename'] = ls
            bbpath = lspath +  '/BBLists'
            specstr = ''
        else:
            if 'enumerate_specstr' in st.session_state:
                specstr = st.session_state['enumerate_specstr']
                slash_specstr = '/' + specstr
            else:
                specstr = ''
                slash_specstr = ''
            bbpath = lspath + st.session_state['enumerate_schemename'] + slash_specstr +  '/BBLists'

        for s in schemejson:
            schemelist.append(s)
        schemelist.sort ()

        dfs = self.UpdateBBDfs(bbpath, False)
        with st.expander (label='Scheme Definition Tools'):
            cont1 = st.container ()
        cont2 = st.container()
        cont3= st.container ()
        cont4 = st.container ()
        cont5 = st.container()


        Enumerate = False
        with cont1:
            ccont = st.container()
            ccontA = st.container ()
            with ccontA:
                newname = st.text_input(label = 'New Scheme Name' )
                if st.button (label = 'Add New Scheme'):
                    st.session_state['enumerate_schemename']= newname
                    ls = newname
                    st.session_state ['enumerate_lastschemedef'] = '{\n"altnames":[],\n"steps":[\n{"stepname":"",\n"Reactants":["",""],\n"Rxns":{\n"default":[]\n}\n}\n],"scaffold_dummy_structures": [], "filters":{"names":{}, "BBfilters":{}}}'
                    self.SaveScheme()
                    self.SetScheme()
                    schemelist.append (newname)
                    dfs = []



        with cont3:
            colx1, colx2 = st.columns(2)
            with colx1:
                if ls == '' or ls not in schemelist:
                    lsidx = 0
                else:
                    lsidx = schemelist.index(ls)
                if  'spec' not in st.session_state and specstr is not None and specstr != '':
                    st.session_state['spec'] = specstr


                schemename = st.selectbox(label='Scheme', options=schemelist, key= 'schemename', index=lsidx)
                specstr = st.text_input(key='spec', label='specstr')


                if schemename != st.session_state['enumerate_schemename'] or specstr != st.session_state['enumerate_specstr']:
                    addspec = ''
                    if specstr != '' and specstr is not None:
                        addspec = '/' + specstr
                    dfs  = self.UpdateBBDfs( st.session_state['enumerate_lastschemepath'] + schemename + addspec +'/BBLists', True)
                    st.session_state['enumerate_schemename'] = schemename
                    st.session_state['enumerate_specstr'] = specstr
                    st.session_state['aggriddata'] = None
                    self.SaveToInit()
                    self.SetScheme()
                else:
                    if 'enumerate_lastschemedef' not in st.session_state:
                        st.session_state['enumerate_schemename'] = schemename
                        st.session_state['enumerate_specstr'] = specstr
                        self.SetScheme()

        with cont2:
            col1, col2, col3 = st.columns(3)
            with col1:
                if st.button('random'):
                    for dx in range(0, len(dfs)):
                        if dfs[dx] is not None:
                            st.session_state['bb' + str(dx) + 'idx'] = dfs[dx].index[
                                random.randint(0, len(dfs[dx]))]
                    Enumerate = True
            with col2:
                if st.button('enumerate'):
                    Enumerate = True

        for n in range(0, len(dfs)):
            if 'bb' + str(n) + 'idx' not in st.session_state:
                st.session_state['bb' + str(n) + 'idx'] = 0
            if 'bb' + str(n) + 'txt' not in st.session_state:
                st.session_state['bb' + str(n) + 'txt'] = ''
        with cont5:
            with st.expander(label='Test structure grid', expanded=True):
                getr100 = st.button(label='get random 100')
                if getr100:
                    with st.spinner('Enumerating'):
                        print(schemename, specstr, specstr, lspath, specstr, rxnschemefile)
                        resdf = self.enum.EnumFromBBFiles(schemename, specstr, specstr, lspath, specstr,
                                                          100, rxnschemefile,
                                                          SMILEScolnames=self.smiles_colnames,
                                                          BBcolnames=self.bbid_colnames,
                                                          rem_dups=False, returndf=True)

                        cols =  ['full_smiles']
                        for c in resdf.columns [0:len(resdf.columns) -1]:
                            cols.append (c)
                        resdf = resdf.loc [:, cols]
                        st.session_state['aggriddata'] = resdf

                if 'aggriddata' in  st.session_state and st.session_state ['aggriddata'] is not None:
                    gb = GridOptionsBuilder.from_dataframe(st.session_state['aggriddata'])

                    gb.configure_selection('single')
                    gridOptions = gb.build()
                    selected = AgGrid(st.session_state['aggriddata'], height=250, update_mode='SELECTION_CHANGED',
                                              gridOptions= gridOptions,  data_return_mode=DataReturnMode.AS_INPUT, reload_data=True)
                    if len(selected['selected_rows']) > 0:
                        for n in range(0, len(dfs) ):
                            st.session_state['bb' + str(n) + 'txt'] = selected['selected_rows'][0]['bb' + str (n+1) + '_smiles' ]
                        Enumerate = True

        with cont3:
            with colx1:
                rxtnts = [None] * len(dfs)

                for n in range (0, len(dfs)) :
                    df = dfs[n]
                    rxtnts[n] = st.text_input(key='bb' + str(n) + 'txt', label='BB' + str(n + 1) + ' (override)')
                    if df is not None:
                        if rxtnts[n] == '':
                             rxtnts[n] = st.selectbox(label='BB' + str (n + 1), options=dfs[n][smilescol]
                                 , key='bb' + str (n ),
                                   index=int(st.session_state['bb' + str(n) + 'idx']))
        with cont1:
            ccol1, ccol2 = st.columns(2)
            st.text ('Current Scheme: '  + st.session_state['enumerate_schemename'])
            with ccol1:
                if st.button (label='revert'):
                    self.SetScheme()
            with ccol2:
                if st.button (label='save scheme'):
                    self.SaveScheme()
            with ccont:
                dcol1, dcol2 = st.columns([2,1] )
                with dcol2:
                    if st.button (label = 'Add Step'):
                        print (st.session_state ['enumerate_lastschemedef'])
                        defjson = json.loads(st.session_state ['enumerate_lastschemedef'] )
                        stepdict = {}
                        stepdict ["Reactants"] = []
                        stepdict["Rxns"] = {'default':''}
                        defjson ['steps'].append (stepdict)
                        st.session_state['enumerate_lastschemedef'] = json.dumps (defjson, indent=4)
                with dcol1:
                    st.text_area(height=200, label='Scheme Definition', key='enumerate_lastschemedef')

        with cont2:
            with st.expander (label='Export'):
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
                        self.enum.EnumFromBBFiles(schemename, specstr, specstr, lspath,  specstr, ct, rxnschemefile, SMILEScolnames=self.smiles_colnames, BBcolnames=self.bbid_colnames, rem_dups=remdups, retIntermeds=addintermeds)


        if Enumerate == True:
            with cont3:
                with colx2:
                    if rxtnts is None or len(rxtnts) == 0:
                        st.text('Reactants not correctly loaded')
                    else:
                        try:
                            res, intermeds, rxninfo= self.enum.TestReactionScheme(schemename, rxtnts, st.session_state ['enumerate_lastschemedef'], True)
                            if res is None or res == 'FAIL' or res.startswith ('FAIL') or res == 'NOSCHEMEFAIL':
                                st.text ('Reaction Failure: ' + res)
                            else:
                                try:
                                    st.pyplot(MolDisplay.ShowMol(res))
                                except Exception as e2:
                                    st.text ('EXCEPTION:',e2)


                        except Exception as e:
                            st.text ('Error: bad scheme definition')
                            for r in rxtnts:
                                st.text (r)
                            exc_type, exc_obj, exc_tb = sys.exc_info()
                            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                            st.text (str(exc_type) + ' ' +  fname + ' ' + str( exc_tb.tb_lineno))
                            st.text (e)
            with cont4:
                coly1, coly2 = st.columns (2)
                with coly1:
                    with st.expander(label='Reagents'):
                        st.image(MolDisplay.ShowMols(rxtnts, cols=1, subImgSize=(400, 200), outtype='b64_datauri'))
                with coly2:
                    with st.expander(label='Reaction Products'):
                        st.image(MolDisplay.ShowMols(intermeds, cols=1, subImgSize=(400, 200), outtype='b64_datauri'))
                displist = []
                print (len(rxtnts), len(intermeds), len(rxninfo))
                for x in range (0, len(rxtnts)):
                    if len (intermeds) - 1 < x:
                        displist.append([rxtnts[x], None, rxninfo[x]])
                    else:
                        displist.append ([rxtnts [x], intermeds[x], rxninfo [x]])
                dispdf = pd.DataFrame (displist, columns = ['Reactant', 'Intermediate','RxnInfo'])
                with st.expander (label='Rxn Info Grid'):
                    MolDisplay.ShowMols_StreamLit_Grid (dispdf,['Reactant', 'Intermediate'], rowheight = 200 )


    def RunUI(self):
        self.head()
        self.body()

if __name__=="__main__":
    if runtime.exists ():
        enum = EnumerationUI ()
        enum.RunUI ()


