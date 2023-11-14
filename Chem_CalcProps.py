import pandas as pd
from rdkit import  Chem
from rdkit.Chem import Descriptors
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import json
import matplotlib.pyplot as plt
import math
from rdkit.Chem import AllChem

NUM_WORKERS = 16

class Chem_CalcProps:
    def addPropsToFile (self, infile, outfilename,  smiles_col = 'SMILES', bbidcols = ['bb1', 'bb2','bb3']):
        def taskfcn(row, smiles_col):
            rowvals = []

            if row [smiles_col] == 'FAIL':
                for b in bbidcols:
                    rowvals.append (row [b])
                for b in bbidcols:
                    rowvals.append ( row[b + '_smiles'])
                rowvals.append (row[smiles_col]),
                rowvals.append ([None, None, None, None, None, None, None, None, None, None])
                res = [None, None, None, None, None, None, None, None, None, None]
            else:
                try:
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
                except:
                    print ('FAIL: ', row[smiles_col])
                    res = [None, None, None, None, None, None, None, None, None, None]
            return res

        df = pd.read_csv(infile)

        ddf = dd.from_pandas(df, npartitions=NUM_WORKERS)
        pbar = ProgressBar()
        pbar.register()
        meta_dict = {0:float, 1:float, 2:float, 3:int, 4:str, 5:int, 6:int, 7:int, 8:float, 9:float}
        res = ddf.apply(taskfcn, axis=1, result_type='expand', args=([smiles_col]), meta=meta_dict).compute()
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
    def CalcPartialCharge (self, SMILES, atomidx):
        m = AllChem.MolFromSmiles(SMILES)
        AllChem.ComputeGasteigerCharges(m)
        res = m.GetAtomWithIdx(atomidx).GetProp('_GasteigerCharge')
        return res

    def Filters(self, m, filter_dict, stats):
        if stats is None:
            stats = {'MW': 0, 'tPSA': 0, 'logP': 0, 'rotB':0, 'hbd':0, 'hba':0}
        if filter_dict is None:
            return False, stats;
        if ('MW_high' in filter_dict) or ('MW_low' in filter_dict):
            mw = Descriptors.MolWt(m)
        if 'MW_high' in filter_dict:
            if mw > filter_dict['MW_high']:
                stats['MW'] += 1
                return True, stats
        if 'MW_low' in filter_dict:
            if mw < filter_dict['MW_low']:
                stats['MW'] += 1
                return True, stats
        if ('tPSA' in filter_dict):
            tpsa = Descriptors.TPSA(m)
            if tpsa > filter_dict['tPSA']:
                stats['tPSA'] += 1
                return True, stats
        if ('logP_high' in filter_dict) or ('logP_low' in filter_dict):
            logP = Descriptors.MolLogP(m)
        if ('logP_high' in filter_dict):
            logphigh = filter_dict['logP_high']
            if logP > logphigh:  # aLogP
                # print ('logphigh', logP)
                stats['logP'] += 1
                return True, stats
        if ('logP_low' in filter_dict):
            logplow = filter_dict['logP_low']
            if logP < logplow:  # aLogP
                stats['logP'] += 1
                return True, stats
        if ('RotB_low' in filter_dict) or ('RotB_high' in filter_dict):
            rotB = rdMolDescriptors.CalcNumRotatableBonds(m)
        if ('RotB_low' in filter_dict):
            rotb_low = filter_dict['RotB_low']
            if rotB < rotb_low:  # aLogP
                stats['rotB'] += 1
                return True, stats
        if ('RotB_high' in filter_dict):
            rotb_high = filter_dict['RotB_high']
            if rotB > rotb_high:  # aLogP
                stats['rotB'] += 1
                return True, stats
        if ('HBD_low' in filter_dict) or ('HBD_high' in filter_dict):
            hbd = rdMolDescriptors.CalcNumHBD(m)
        if ('HBD_low' in filter_dict):
            hbd_low = filter_dict['HBD_low']
            if hbd < hbd_low:
                stats['hbd'] += 1
                return True, stats
        if ('HBD_high' in filter_dict):
            hbd_high = filter_dict['HBD_high']
            if hbd > hbd_high:
                stats['hbd'] += 1
                return True, stats
        if ('HBA_low' in filter_dict) or ('HBA_high' in filter_dict):
            hba = rdMolDescriptors.CalcNumHBA(m)
        if ('HBA_low' in filter_dict):
            hba_low = filter_dict['HBA_low']
            if hba < hba_low:
                stats['hba'] += 1
                return True, stats
        if ('HBA_high' in filter_dict):
            hba_high = filter_dict['HBA_high']
            if hba > hba_high:
                stats['hba'] += 1
                return True, stats

        return False, stats


# props = Chem_CalcProps ()

if __name__ == "__main__":
    props = Chem_CalcProps()
    res = props.CalcPartialCharge ('O=CCCNC(=O)OCc1ccccc1', 1)
    print (res)