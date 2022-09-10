from rdkit import Chem
from rdkit.Chem import rdChemReactions
import json
import toml
import pandas as pd
from numpy import random
from multiprocessing.pool import ThreadPool as Pool
import time
import threading
import MolDisplay
import ChemUtilities
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import os
from rdkit.Chem import SaltRemover
import pathlib
from tqdm import tqdm
import streamlit as st

class Enumerate:
    rxndict = {}
    def React_Molecules (self, r1, r2, SM_rxn, showmols):
        if showmols == True:
            print ('Reaction Params:')
            print(r1, r2, SM_rxn)
        if isinstance (SM_rxn,str):
            SM_rxn = [SM_rxn]

        if r2 == 'None' or r2 is None or r2 == 'NOSMILES' or r2 == "SKIPCYC":
            if not isinstance(r1, Chem.rdchem.Mol):
                m1 = Chem.MolFromSmiles(r1)
            else:
                m1 = r1
            reacts = [m1]
        else:
            if not isinstance(r1, Chem.rdchem.Mol):
                m1 = Chem.MolFromSmiles(r1)
            else:
                m1 = r1
            if not isinstance(r2, Chem.rdchem.Mol):
                m2 = Chem.MolFromSmiles(r2)
            else:
                m2 = r2
            reacts = (m1, m2)

        for r in SM_rxn:
            rxn = rdChemReactions.ReactionFromSmarts(r)
            products = rxn.RunReactants(reacts)
            if len(products) == 0:
                product = None
            else:
                product = products [0][0]

            smi = None

            try:
                if product is  None:
                    smi = None
                else:
                    smi = Chem.MolToSmiles(product)
                    break
            except:
                smi = None

        if showmols == True:
            if smi is not None:
                MolDisplay.ShowMols([smi, r1, r2])
            else:
                MolDisplay.ShowMols([ r1, r2])

        return smi, len(products)

    def Show_Reactants (self, reactants):
        for r in reactants:
            MolDisplay.ShowMol (r)

    def ReadRxnScheme (self, fname, schemename, verbose = True, FullInfo = False):
        try:
            f = open (fname)
            data = json.load (f)
            f.close ()
        except Exception as E:
            if verbose:
                print ('json read error:', str(E))
            return None, None
        if schemename not in data:
            if verbose:
                print('schemename error:', schemename)
            return None, None
        if 'reactants' not in data[schemename]:
            reactants = None
        else:
            reactants = data[schemename]["reactants"]
        if FullInfo:
            return data[schemename], reactants
        else:
            return data[schemename]["steps"], reactants

    def Pull_ReactionDetails ( self, prior_product, in_reactants, step, rxtants):
        reactants = []

        def getReactants (rx, reactants):
            if step["Reactants"][rx] == 'None':
                reactants.append ('None')
            elif step["Reactants"][rx] == 'NOSMILES':
                reactants.append ('None')
            elif step["Reactants"][rx] == 'SKIPCYC':
                reactants.append('None')
            elif step["Reactants"][rx] == 'p':
                reactants.append(prior_product)
            elif rxtants != None and step["Reactants"][rx] in rxtants:
                idx = rxtants.index (step["Reactants"][rx])
                reactants.append(in_reactants[idx])
            elif ('r') in step["Reactants"][rx]:
                idx = int(step["Reactants"][rx][1:])
                reactants.append(in_reactants[idx])
            else:
                reactants.append(step["Reactants"][rx])
            return reactants

        reactants = getReactants(0, reactants)
        reactants = getReactants(1, reactants)
        if "default" in step["Rxns"] and len (step["Rxns"]) == 1:
            reaction = step["Rxns"]["default"]
        else:
            reaction = 'unassigned'
            for rxn in step["Rxns"]:
                if rxn != 'default':
                    if rxn == 'product':
                        reaction = 'product'
                        break
                    if rxn == 'NOSMILES':
                        if reactants [0] == 'NOSMILES' or reactants [1] == 'NOSMILES':
                            reaction = step["Rxns"][rxn]
                            break
                        else:
                            reaction = 'FAIL'
                            continue
                    if rxn == 'SKIPCYC':
                        if reactants[0] == 'SKIPCYC' or reactants[1] == 'SKIPCYC':
                            reaction = step["Rxns"][rxn]
                            break
                        else:
                            reaction = 'FAIL'
                            continue
                    if reactants [1] == 'None':
                        if not isinstance(reactants [0], Chem.rdchem.Mol):
                            m = Chem.MolFromSmiles(reactants [0])
                        else:
                            m = reactants [0]
                    else:
                        if not isinstance(reactants[1], Chem.rdchem.Mol):
                            m = Chem.MolFromSmiles(reactants[1])
                        else:
                            m = reactants[1]
                    pattern = Chem.MolFromSmarts(rxn)
                    if m is None:
                        reaction = 'FAIL'
                    else:
                        if m.HasSubstructMatch(pattern):
                            reaction =  step["Rxns"][rxn]
                            break
                        else:
                            reaction ='FAIL'
            if reaction == 'FAIL' and "default" in step["Rxns"] :
                reaction = step["Rxns"]["default"]
        return reactants, reaction

    def Read_ReactionDatabase (self, rxnfile):
        rxndb = toml.load (rxnfile)
        self.rxndict = {}
        for rtype in rxndb:
            subtypes = rxndb[rtype]
            for subtype in subtypes:
                self.rxndict [subtype]= subtypes [subtype]['smarts']

    def RunRxnScheme (self, in_reactants, schemefile, schemename, showmols, schemeinfo = None):

        if schemeinfo is not None:
            scheme, rxtants = [schemeinfo[0], schemeinfo [1]]
        else:
            scheme, rxtants = self.ReadRxnScheme(schemefile, schemename)

        if scheme is None:
            return 'NOSCHEMEFAIL'
        p=in_reactants [0]

        prod_ct = 1
        for step in scheme:
            stepname = step

            if type (scheme) == dict:
                step = scheme [step]
            reactants, reaction = self.Pull_ReactionDetails (p, in_reactants, step, rxtants)

            if reaction == 'FAIL':
                p = 'FAIL'
                break
            if reaction == 'terminate':
                break
            if reaction == 'product':
                continue
            else:
                p, outct = self.React_Molecules(reactants [0], reactants [1],  reaction,  showmols)
                prod_ct *= outct
            if showmols:
                print(stepname, p, reactants)
            if p is None:
                p = 'FAIL'
                break
        return p, prod_ct, [scheme, rxtants]

    def Sample_Library (self, BBlistpath, outfilepath, schemepath, scheme_name,  rndmsample_ct, ShowMols, saltstrip = True, CountOnly=False ):
        def recurse_library (prev_rvals, prev_idvals, cyc_num, cycct, enum_list, outfile, ct):
            for idx, row in cyc[cycles[cyc_num]].iterrows ():
                rvals = prev_rvals.copy ()
                rvals.append (row ['SMILES'])
                idvals = prev_idvals.copy ()
                idvals.append(str(row ['ID']))

                if cyc_num == cycct -1:
                    try:
                        enum_res, prod_ct = self.RunRxnScheme(rvals, schemepath, scheme_name, ShowMols)
                    except:
                        enum_res = 'FAIL'
                    outfile.write(enum_res)
                    for ir in range(0, len(rvals)):
                        outfile.write(',' + str(rvals[ir]) + ',' + str(idvals[ir]))
                    outfile.write('\n')
                    ct += 1

                    if ct %1000 ==0:
                        print (ct)
                else:
                    ct = recurse_library (rvals, idvals, cyc_num + 1, cycct, enum_list, outfile, ct)
            return ct

        def EnumFullMolecule (cycles, cycvals, r_list):
            rvals = []
            idvals = []
            resline = ''

            for c_i in range(0, len(cycles)):
                rnum = random.randint (0,len(cycvals [c_i]))
                rcx = cycvals [c_i].iloc [rnum]
                rvals.append (rcx['SMILES'])
                idvals.append (rcx ['ID'])
            string_ids = [str(int) for int in idvals]
            idstr = ','.join(string_ids)
            if not  idstr  in r_list:
                 resline = ''
                 r_list[idstr] = True
                 enum_res, prod_ct = self.RunRxnScheme(rvals, schemepath, scheme_name, ShowMols)
                 resline = enum_res
                 for ir in range(0, len(rvals)):
                      resline +=',' + rvals[ir] + ',' + str(idvals[ir])
            else:
                return None
            return resline

        class CompleteEnum ():
            tic = None
            block = ''
            lock = None
            outfile = None
            def __init__ (self, outf):
                self.tic = time.perf_counter()
                self.lock = threading.Lock()
                self.outfile = outf
            enumct = 0

            def CompleteEnumAsync (self, res):
                if res is not None:
                    self.lock.acquire()
                    self.enumct = self.enumct + 1
                    self.block += res
                    self.block += '\n'
                    if self.enumct % 1000 == 0:
                        print(self.enumct)
                        toc = time.perf_counter()
                        print(f"Time Count {toc - self.tic:0.4f} seconds")
                        self.WriteBlock ()
                    self.lock.release()

            def WriteBlock (self):
                self.outfile.write(self.block)
                self.block = ''

        df = pd.read_csv (BBlistpath)
        df = df.fillna(0)

        scheme, reactants = self.ReadRxnScheme(schemepath, scheme_name)
        if reactants is None:
            cycles = df.Cycle.unique ()
            cycles.sort ()
        else:
            cycles = reactants

        cyc = {}
        if saltstrip == True:
            df = ChemUtilities.SaltStripMolecules(df)

        for c in cycles:
            cyc[c] = df[df.Cycle.isin([c])]

        if CountOnly:
            sz = 1
            for cx in cyc:
                sz *= len (cyc[cx])
                print (cx, str(len (cyc[cx])))
            print ('libsize = ', sz)
            return ()

        ct = 0
        outfile = open(outfilepath, 'w')
        enum_list = []

        outfile = open(outfilepath, 'w')

        outfile.write('product_SMILES')

        for c_i in range(0, len(cycles)):
            outfile.write(',' + str(cycles [c_i]) + ',ID_' +  str(cycles [c_i]) )

        outfile.write('\n')

        if rndmsample_ct == -1:
            rvals = []
            idvals = []
            recurse_library (rvals, idvals, 0, len(cycles), enum_list, outfile, 0)
        else:
            cycvals = []
            for c_i in range (0,len(cycles)):
                cycvals.append (cyc[cycles[c_i]])

            r_list = {}
            block = ''
            self.enumct = 0
            CE = CompleteEnum (outfile)
            pool_size = 1
            pool = Pool(pool_size)
            for r in range (0, rndmsample_ct):
                pool.apply_async(EnumFullMolecule, args=(cycles, cycvals, r_list), callback=CE.CompleteEnumAsync)
            pool.close ()
            pool.join ()
            CE.WriteBlock()


        outfile.close ()

    def Deprotect (self, infile, deprotect_specfile, dp_outfile, smilescol, replace):
        if deprotect_specfile is None:
            return
        df = pd.read_csv (infile)
        dp_list = pd.DataFrame (columns=df.columns)
        for idx, row in df.iterrows ():
            m = Chem.MolFromSmiles(row[smilescol])
            try:
                res, deprotected = ChemUtilities.Deprotect(m, deprotect_specfile, True)
                if deprotected == True :
                    if replace == True:
                        l = row[smilescol]
                        df.at [idx, smilescol] =  Chem.MolToSmiles(res[0])
                    else:
                        r2 = row.copy ()
                        r2 [smilescol] =  Chem.MolToSmiles(res[0])
                        dp_list = dp_list.append(r2)
            except:
                deprotected = False
                res = [m]
        if replace:
            df.to_csv (dp_outfile)
        else:
            df =  df.append (dp_list)
            df.to_csv (dp_outfile)

    def enumerate_library_strux(self, libname, rxschemefile, infilelist, outpath, rndct=-1, bblistfile=None, SMILEScolnames = [], BBIDcolnames = []):
    # currently hard coded to 3 cycles, needs to be improved
        def taskfcn(row, libname, rxschemefile, showstrux, schemeinfo):
            rxtnts = [row['bb1_smiles'], row['bb2_smiles'], row['bb3_smiles']]
            try:
                res, prodct, schemeinfo = self.RunRxnScheme(rxtnts, rxschemefile, libname, showstrux, schemeinfo)
            except:
                print('FAIL--Error', rxtnts)
                return 'FAIL'
            if res == 'FAIL':
                print('FAIL in taskfcn')
                print(rxtnts)

            return res

        if type(infilelist) is list:
            cycdict = {}
            ct = 1
            for f in infilelist:
                cycdict['BB' + str(ct)] = pd.read_csv(f)
                changecoldict = {}
                for smc in SMILEScolnames:
                    if smc in cycdict['BB' + str(ct)].columns:
                        changecoldict[smc] = 'SMILES'
                for bbidc in BBIDcolnames:
                    if bbidc in cycdict['BB' + str(ct)].columns:
                        changecoldict[bbidc] = 'BB_ID'
                ct += 1

            if bblistfile is not None:
                picklistdf = pd.read_csv(bblistfile)
                b1df = picklistdf[picklistdf['Cycle'] == 'BB1']
                b1df = b1df.merge(cycdict['BB1'], on=['BB_ID'], how='left')
                b2df = picklistdf[picklistdf['Cycle'] == 'BB2']
                b2df = b2df.merge(cycdict['BB2'], on=['BB_ID'], how='left')
                b3df = picklistdf[picklistdf['Cycle'] == 'BB3']
                b3df = b3df.merge(cycdict['BB3'], on=['BB_ID'], how='left')
            else:
                b1df = cycdict['BB1']
                b2df = cycdict['BB2']
                b3df = cycdict['BB3']

            ct = 0
            fullct = len(b1df) * len(b2df) * len(b3df)
            print('molecules to enumerate:', rndct)
            if rndct > len(b1df) * len(b2df) * len(b3df):
                rndct = -1
            if rndct == -1:
                reslist = [[]] * fullct
                for i1, b1 in b1df.iterrows():
                    for i2, b2 in b2df.iterrows():
                        for i3, b3 in b3df.iterrows():
                            bb1 = b1['BB_ID']
                            bb1s = b1['SMILES']
                            bb2 = b2['BB_ID']
                            bb2s = b2['SMILES']
                            bb3 = b3['BB_ID']
                            bb3s = b3['SMILES']
                            reslist[ct] = [bb1, bb1s, bb2, bb2s, bb3, bb3s]
                            ct += 1
                            if ct % 10000 == 0:
                                print(ct, '/', fullct, end='\r')
            else:
                reslist = [[]] * rndct
                ct = 0

                while ct < rndct:
                    bb1df = cycdict['BB1'].sample()
                    b1 = bb1df.iloc[0]
                    bb1 = b1['BB_ID']
                    bb1s = b1['SMILES']
                    bb2df = cycdict['BB2'].sample()
                    b2 = bb2df.iloc[0]
                    bb2 = b2['BB_ID']
                    bb2s = b2['SMILES']
                    bb3df = cycdict['BB3'].sample()
                    b3 = bb3df.iloc[0]
                    bb3 = b3['BB_ID']
                    bb3s = b3['SMILES']
                    if [bb1, bb1s, bb2, bb2s, bb3, bb3s] not in reslist:
                        reslist[ct] = [bb1, bb1s, bb2, bb2s, bb3, bb3s]
                        ct += 1
                        if ct % 1000 == 0:
                            print(ct, '/', fullct, end='\r')
            resdf = pd.DataFrame(reslist, columns=['bb1', 'bb1_smiles', 'bb2', 'bb2_smiles', 'bb3', 'bb3_smiles'])
        else:
            resdf = infilelist

        print('starting enumeration')
        NUM_WORKERS = 16

        pbar = ProgressBar()
        pbar.register()
        ddf = dd.from_pandas(resdf, npartitions=NUM_WORKERS)
        schemeinfo = self.ReadRxnScheme(rxschemefile, libname, False)
        res = ddf.apply(taskfcn, axis=1, result_type='expand', args=(libname, rxschemefile, rndct == 1, schemeinfo),
                        meta=(0, str)).compute()
        pbar.unregister()

        df = resdf.merge(res, left_index=True, right_index=True)
        df = df.rename(columns={0: 'full_smiles'})
        outsuff = 'full'
        if rndct != -1:
            outsuff = str(rndct)
        if outpath is not None:
            df[df['full_smiles'] != 'FAIL'].to_csv(outpath + '.' + outsuff + '.enum.csv', index=False)
            df[df['full_smiles'] == 'FAIL'].to_csv(outpath + '.' + outsuff + '.fail.csv', index=False)
            df.to_csv(outpath + '.' + outsuff + '.all.csv', index=False)
            return outpath + '.' + outsuff + '.all.csv'
        else:
            return df

    def EnumFromBBFiles(self, libname, bbspec, outspec, inpath, foldername, num_strux, rxschemefile, picklistfile=None):
        if rxschemefile is None:
            rxschemefile = inpath + 'RxnSchemes.json'
        if bbspec != '':
            bbspec = '.' + bbspec
        infilelist = [
            inpath + foldername + '/BBLists/' + libname + bbspec + '.BB1.csv',
            inpath + foldername + '/BBLists/' + libname + bbspec + '.BB2.csv',
            inpath + foldername + '/BBLists/' + libname + bbspec + '.BB3.csv'
        ]
        print(infilelist)
        samplespath = inpath + foldername + '/Samples/'
        if not os.path.exists(samplespath):
            os.makedirs(samplespath)
        outpath = inpath + foldername + '/Samples/' + libname + outspec
        print('outpath:', outpath)
        print('outpath:', outpath)
        outfile = self.enumerate_library_strux(libname, rxschemefile, infilelist, outpath, num_strux, picklistfile)
        return outfile

    def FilterBBs(self, bbdict, filterfile):
        filtersdf = pd.read_csv(filterfile)
        patterndict = {}
        removeidxs = {}
        for ix, r in filtersdf.iterrows():
            pattern = r['smarts']
            patterndict[pattern] = Chem.MolFromSmarts(pattern)

        for dx in bbdict.keys():
            print(len(bbdict[dx]))
            removeidxs[dx] = []
            for idx, row in bbdict[dx].iterrows():
                m = Chem.MolFromSmiles(row['SMILES'])
                for v in patterndict.values():
                    if m.HasSubstructMatch(v) == True:
                        removeidxs[dx].append(idx)
                        continue
            bbdict[dx].drop(removeidxs[dx], axis=0, inplace=True)
            print(removeidxs[dx], len(removeidxs[dx]), len(removeidxs[dx]), len(bbdict[dx]), len(patterndict.values()))
        return bbdict

    def generate_BBlists(self, inpath, outpath, libname, rxnschemefile):
        def PassFilters(input_mol, filters_dict, filtname):
            filters = filters_dict[filtname]
            incl_list = filters['include']
            if 'exclude' in filters:
                excl_list = filters['exclude']
            else:
                excl_list = []
            for ifl in incl_list:
                if type(ifl) == list:
                    hassub = False
                    for orx in ifl:
                        pattern = Chem.MolFromSmarts(orx)
                        if (input_mol.HasSubstructMatch(pattern) == True):
                            hassub = True
                    if hassub == False:
                        return False
                else:
                    pattern = Chem.MolFromSmarts(ifl)
                    if (input_mol.HasSubstructMatch(pattern) == False):
                        return False
            for efl in excl_list:
                pattern = Chem.MolFromSmarts(efl)
                if (input_mol.HasSubstructMatch(pattern) == True):
                    return False
            return True

        f = open(rxnschemefile)
        data = json.load(f)
        f.close()
        names_dict = data[libname]['filters']['names']
        listnames = names_dict.keys()
        filters_dict = data[libname]['filters']['BB_filters']
        if 'GeneralFilters' in data:
            GenFilters_dict = data['GeneralFilters']
        else:
            GenFilters_dict = {}

        for n in names_dict:
            if n not in filters_dict:
                if n not in GenFilters_dict:
                    print('Error: Filter not found')
                else:
                    filters_dict[n] = GenFilters_dict[n]
        print(filters_dict)

        print('pull bbs')
        df = self.pull_BBs(inpath)
        print('pull bbs complete')
        for ln in listnames:
            df[ln] = None

        saltstrip = SaltRemover.SaltRemover()
        print('loop start')
        for idx, row in df.iterrows():
            try:
                input_mol = Chem.MolFromSmiles(row['SMILES'])
                input_mol = saltstrip.StripMol(input_mol)
                df.iloc[idx]['SMILES'] = Chem.MolToSmiles(input_mol)
            except:
                input_mol = None
            if input_mol is not None:
                for ln in listnames:
                    if PassFilters(input_mol, filters_dict, ln):
                        df.iloc[idx][ln] = 'Y'
                    else:
                        df.iloc[idx][ln] = 'N'
        print('loop done')
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        for ln in listnames:
            cdf = df[df[ln] == 'Y'].drop_duplicates(subset=['BB_ID'])
            if type(names_dict[ln]) == list:
                for ix in range(0, len(names_dict[ln])):
                    cdf.to_csv(outpath + '/' + libname + '.' + names_dict[ln][ix] + '.csv', index=False)
            else:
                cdf.to_csv(outpath + '/' + libname + '.' + names_dict[ln] + '.csv', index=False)

    def load_BBlists(self, bblists, BBIDcolnames= [], SMILEScolnames = []):
        cycs = []
        bbdict = {}

        if type (bblists) is dict:
            bbdict = bblists
        else:
            for l in bblists:
                cyclesplit = l.split('.')
                cyc = cyclesplit[len(cyclesplit) - 2]
                cycs.append(cyc)
                if type(l) is str:
                    bbdict[cyc] = pd.read_csv(l)

        for cyc in bbdict.keys():
            changecoldict = {}
            for smc in SMILEScolnames:
                if smc in bbdict[cyc].columns:
                    changecoldict[smc] = 'SMILES'
            for bbidc in BBIDcolnames:
                if bbidc in bbdict[cyc].columns:
                    changecoldict[bbidc] = 'BB_ID'
            if len(changecoldict) > 0:
                bbdict[cyc] = bbdict[cyc].rename(columns=changecoldict)
            print(bbdict[cyc].columns)
            bbdict[cyc] = bbdict[cyc].drop_duplicates(subset='BB_ID', keep="first")
        return bbdict

    def pull_BBs(self, inpath, idcol, smilescol):
        globlist = pathlib.Path(inpath).glob('*.csv')
        full_df = None
        for f in tqdm(globlist):
            print(f)
            if full_df is None:
                full_df = pd.read_csv(f)
                full_df = full_df[[idcol, smilescol]]
                full_df = full_df.rename(columns={idcol: 'BB_ID', smilescol: 'SMILES'})
            else:
                df = pd.read_csv(f)
                df = df.rename(columns={idcol: 'BB_ID', smilescol: 'SMILES'})
                try:
                    df = df[['BB_ID', 'SMILES']].dropna()
                except:
                    print('exiting:', f)
                    exit()
                full_df = full_df.append(df, ignore_index=True)
        return full_df

    def TestReactionScheme(self,schemename, rxtnts, rxnschemefile):
        enum = Enumerate.Enumerate()
        res = enum.RunRxnScheme(rxtnts, rxnschemefile, schemename, True)
        return res[0]

class EnumerationUI():
    enum = Enumerate()

    def head(self):
        st.markdown("""
            <h1 style='text-align: center; margin-bottom: -35px;'>
            Enumeration Tester
            </h1>
        """, unsafe_allow_html=True
                    )

    def UpdateBBDfs(self, inpath, do_reload):
        print(inpath)
        flist = pathlib.Path(inpath).glob('*.csv')
        bb1file = ''
        bb2file = ''
        bb3file = ''

        for f in flist:
            c = str(f)
            if '.BB1.' in c:
                bb1file = c
            if '.BB2.' in c:
                bb2file = c
            if '.BB3.' in c:
                bb3file = c

        df1 = None
        if 'df1' not in st.session_state or do_reload:
            if bb1file != '':
                df1 = pd.read_csv(bb1file)
                st.session_state['df1'] = df1
        else:
            df1 = st.session_state['df1']
        df2 = None
        if 'df2' not in st.session_state or do_reload:
            if bb2file != '':
                df2 = pd.read_csv(bb2file)
                st.session_state['df2'] = df2
        else:
            df2 = st.session_state['df2']
        df3 = None
        if 'df3' not in st.session_state or do_reload:
            if bb3file != '':
                df3 = pd.read_csv(bb3file)
                st.session_state['df3'] = df3
        else:
            df3 = st.session_state['df3']

        return df1, df2, df3

    def body(self, smilescol, rootpath):
        rxnschemefile = st.text_input(value=rootpath + 'RxnSchemes.json', label='rxscheme')

        col1, col2 = st.columns(2)
        lspath = st.text_input(value=rootpath, label='scheme path')
        ls = lspath.replace(rootpath, '')
        bbpath = lspath + '/BBLists'
        f = open(rxnschemefile, 'r')
        schemejson = json.load(f)
        schemelist = []

        if 'schemename' not in st.session_state:
            print(ls)
            st.session_state['schemename'] = ls

        for s in schemejson:
            schemelist.append(s)
        df1, df2, df3 = self.UpdateBBDfs(bbpath, False)

        if 'bb1idx' not in st.session_state:
            st.session_state['bb1idx'] = 0
        if 'bb2idx' not in st.session_state:
            st.session_state['bb2idx'] = 0
        if 'bb3idx' not in st.session_state:
            st.session_state['bb3idx'] = 0

        if st.button('random'):
            if df1 is not None:
                st.session_state['bb1idx'] = df1.index[random.randint(0, len(df1))]
            if df2 is not None:
                st.session_state['bb2idx'] = df2.index[random.randint(0, len(df2))]
            if df3 is not None:
                st.session_state['bb3idx'] = df3.index[random.randint(0, len(df3))]

        with col1:
            if ls not in schemelist:
                lsidx = 0
            else:
                lsidx = schemelist.index(ls)
            schemename = st.selectbox(label='Scheme', options=schemelist, key='scheme', index=lsidx)
            if schemename != st.session_state['schemename']:
                df1, df2, df3 = self.UpdateBBDfs(rootpath + schemename + '/BBLists', True)
                print(len(df2), len(df3), schemename)

            st.session_state['schemename'] = schemename

            if df1 is not None:
                rxtnt1 = st.text_input(key='bb1txt', label='bb1 (override)')
                if rxtnt1 == '':
                    rxtnt1 = st.selectbox(label='bb1', options=df1[smilescol], key='bb1',
                                          index=st.session_state['bb1idx'])
            if df2 is not None:
                rxtnt2 = st.text_input(key='bb2txt', label='bb2 (override)')
                if rxtnt2 == '':
                    rxtnt2 = st.selectbox(label='bb2', options=df2[smilescol], key='bb2',
                                          index=st.session_state['bb2idx'])
            if df3 is not None:
                rxtnt3 = st.text_input(key='bb3txt', label='bb3 (override)')
                if rxtnt3 == '':
                    rxtnt3 = st.selectbox(label='bb3', options=df3[smilescol], key='bb3',
                                          index=st.session_state['bb3idx'])

        if st.button('enumerate'):
            rxtnts = [rxtnt1, rxtnt2, rxtnt3]
            res = self.enum.TestReactionScheme(schemename, rxtnts, rxnschemefile)
            st.pyplot(MolDisplay.ShowMols(rxtnts))
            with col2:
                st.write(schemename)
                st.write(res)
                st.pyplot(MolDisplay.ShowMol(res))

        if st.button('Export Random Selection'):
            enumfile = self.enum.EnumFromBBFiles(schemename, '', '', lspath, '', 5000, rxnschemefile)

    def RunUI(self, smilescol, rootpath):
        self.head()
        self.body(smilescol, rootpath)

if __name__=="__main__":
    if st._is_running_with_streamlit:
        enum = EnumerationUI ()
        enum.RunUI ('SMILES', '')