import pandas as pd
from rdkit import  Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import streamlit as st

class CalcProps:
    def addPropsToFile (self, infile, outfile,  smiles_col = 'SMILES', bbidcols = ['bb1', 'bb2','bb3']):
        def taskfcn(row):

            rowvals = []
            if row [smiles_col] == 'FAIL':
                for b in bbcols:
                    rowvals.append (row [b])
                for b in bbcols:
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

            res = [ round(SlogP, 2), round(TPSA, 2), AMW, RBs, CSMILES, HBD, HBA, HAC, round(SP3, 2), round(ExactMW, 2)]

            return res

        df = pd.read_csv(infile).head () #DEBUG
        NUM_WORKERS = 16
        ddf = dd.from_pandas(df, npartitions=NUM_WORKERS)
        pbar = ProgressBar()
        pbar.register()
        meta_dict = {0:'SlogP', 1:'TPSA', 2:'AMW', 3:'RBs', 4:'CSMILES', 5:'HBD', 6:'HBA', 7:'HAC', 8:'SP3', 9:'ExactMW'}
        res = ddf.apply(taskfcn, axis=1, result_type='expand', args=(), meta=meta_dict).compute()
        res.columns = ['SlogP', 'TPSA', 'AMW', 'RBs', 'CSMILES', 'HBD', 'HBA', 'HAC', 'SP3', 'ExactMW']
        df = resdf.merge(res, left_index=True, right_index=True)
        df.to_csv(outfilename, index=False)
        pbar.unregister()

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

props = CalcProps ()
class CalcPropsUI:
    rootpath = ''
    initpath = '../CIxTools.init.json'
    bbdf = None

    def __init__(self):
        return

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
        smiles_col = st.text_input(label='SMILES Column', key='props_smilescol', value='full_smiles',  on_change=self.SaveToInit)
        bbidcols = st.text_input(label='BBID Cols', key='props_bbidcols', value='bb1,bb2,bb3',  on_change=self.SaveToInit).split(',')
        RunClusters = st.button(label='Run', key='RunProps')
        if RunClusters == True:
            outfile = filename.replace ('.csv', '.props.csv')
            if outfile == filename:
                st.text ('Incompatible input file, .csv file is required')
                return
            props.addPropsToFile(filename, outfile,  smiles_col, bbidcols)
        return

    def SaveToInit(self):
        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        if 'props_infile' in st.session_state:
            data["propsinfile"] = st.session_state.propsinfile
        if 'props_smilescol' in st.session_state:
            data["props_smilescol"] = st.session_state.schemepath
        if 'props_bbidcols' in st.session_state:
            data["props_bbidcols"] = st.session_state.props_bbidcols
        with open(self.initpath, "w") as jsonFile:
            json.dump(data, jsonFile)

    def RunUI (self):
        self.head ()
        self.body ()
