import pandas as pd
from rdkit import  Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import streamlit as st
import json
import matplotlib.pyplot as plt
import math

class Chem_CalcProps:
    def addPropsToFile (self, infile, outfilename,  smiles_col = 'SMILES', bbidcols = ['bb1', 'bb2','bb3']):
        def taskfcn(row):
            rowvals = []
            if row [smiles_col] == 'FAIL':
                for b in bbidcols:
                    rowvals.append (row [b])
                for b in bbidcols:
                    rowvals.append ( row[b + '_smiles'])
                rowvals.append (row['full_smiles']),
                rowvals.append ([None, None, None, None, None, None, None, None, None, None])
            else:
                mol = Chem.MolFromSmiles(row[smiles_col])
                TPSA = rdMolDescriptors.CalcTPSA(mol)
                RBs = rdMolDescriptors.CalcNumRotatableBonds(mol)
                HBD = rdMolDescriptors.CalcNumHBD(mol)
                HBA = rdMolDescriptors.CalcNumHBA(mol)
                AMW = Chem.Descriptors.MolWt(mol)
                HAC = mol.GetNumHeavyAtoms()
                SlogP = Chem.Crippen.MolLogP(mol)
                SP3 = Chem.Lipinski.FractionCSP3(mol)
                ExactMW = rdMolDescriptors.CalcExactMolWt(mol)
                CSMILES = Chem.CanonSmiles(row[smiles_col])

            res = [ round(SlogP, 2), round(TPSA, 2), round (AMW,2), RBs, CSMILES, HBD, HBA, HAC, round(SP3, 2), round(ExactMW, 2)]

            return res

        df = pd.read_csv(infile)
        NUM_WORKERS = 16
        ddf = dd.from_pandas(df, npartitions=NUM_WORKERS)
        pbar = ProgressBar()
        pbar.register()
        meta_dict = {0:float, 1:float, 2:float, 3:int, 4:str, 5:int, 6:int, 7:int, 8:float, 9:float}
        res = ddf.apply(taskfcn, axis=1, result_type='expand', args=(), meta=meta_dict).compute()
        res.columns = ['SlogP', 'TPSA', 'AMW', 'RBs', 'CSMILES', 'HBD', 'HBA', 'HAC', 'SP3', 'ExactMW']
        df = df.merge(res, left_index=True, right_index=True)
        df.to_csv(outfilename, index=False)
        pbar.unregister()
        self.GeneratePlots (df, outfilename)

    def GeneratePlots (self, df, infilename, props = {'SlogP':(-1, 10, 12), 'TPSA':(70, 300, 24 ), 'RBs':(0, 15, 16), 'HBD':(0,10, 11), 'HBA':(0, 20, 21), 'HAC':(0, 120, 25), 'SP3':(0, 1, 21), 'ExactMW':(250, 1250,21 ) }):
        if df is None:
            df = pd.read_csv (infilename)
        plotlen = len (props)
        gridwidth = 2

        plt.tight_layout ()
        plt.figure(figsize=(8, 12), dpi=80)
        plt.subplots_adjust(left=0.15, right=.9, top=0.9, bottom=0.1)
        grid = plt.GridSpec(math.ceil (plotlen/gridwidth), gridwidth, wspace=.8, hspace=.8)

        ct = 0
        for p in props.keys ():
            subplt = plt.subplot( grid[int(ct/gridwidth), ct % gridwidth])
            step = (props[p][1]-props[p][0])/(props[p][2]-1)
            numsteps = (props[p][1]-props[p][0])/step
            ticks =  [float(x) for x in range (0, math.ceil(numsteps))]
            ticks = [step * t + props[p][0]  for t in ticks  ]
            subplt.hist(df[p], edgecolor = 'black', color='blue', bins = ticks)
            ct += 1
            subplt.set (xlabel=p, ylabel = 'Count')
            plt.title(p + " Histogram")
        plt.savefig(infilename + '.props_hist.png')

        return infilename + '.props_hist.png'

    def addPropsToFile2 (self, infile, outfile, proplist):
        df = pd.read_csv(infile)
        for p in proplist:
            if p not in df.columns:
                df [p] = None
        for ix, row in df.iterrows ():
            mol = Chem.MolFromSmiles(row['csmiles'])
            if 'TPSA' in proplist:
                TPSA = rdMolDescriptors.CalcTPSA(mol)
                df.loc[ix, 'TPSA'] = TPSA
            if 'RBs' in proplist:
                RBs = rdMolDescriptors.CalcNumRotatableBonds(mol)
                df.loc[ix, 'RBs'] = RBs
            if 'HBD' in proplist:
                HBD = rdMolDescriptors.CalcNumHBD(mol)
                df.loc[ix, 'HBD'] = HBD
            if 'HBA' in proplist:
                HBA = rdMolDescriptors.CalcNumHBA(mol)
                df.loc[ix, 'HBA'] = HBA
            if 'AMW' in proplist:
                AMW = Chem.Descriptors.MolWt(mol)
                df.loc [ix,'AMW'] = AMW
            if 'HAC' in proplist:
                HAC = mol.GetNumHeavyAtoms()
                df.loc[ix, 'HAC'] = HAC
            if 'SlogP' in proplist:
                SlogP = Chem.Crippen.MolLogP(mol)
                df.loc[ix, 'SlogP'] = SlogP
            if 'SP3' in proplist:
                SP3 = Chem.Lipinski.FractionCSP3(mol)
                df.loc[ix, 'SP3'] = SP3
            if 'ExactMW' in proplist:
                ExactMW = rdMolDescriptors.CalcExactMolWt(mol)
                df.loc[ix, 'ExactMW'] = ExactMW
            if 'CSMILES' in proplist:
                CSMILES = Chem.CanonSmiles(row['full_smiles'])
                df.loc[ix, 'CSMILES'] = CSMILES
        df.to_csv (outfile)

props = Chem_CalcProps ()

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
                props.addPropsToFile(filename, outfile,  smiles_col, bbidcols)
        with st.expander ('Property Graphs'):
            infile = st.text_input ('property file path', key='props_filepath', on_change=self.SaveToInit)
            gen_graphs = st.button ('Generate Graphs')
            if gen_graphs:
                with st.spinner ('Generating Plots...'):
                    res  = props.GeneratePlots(None, infile )
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
