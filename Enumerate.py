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
import matplotlib.pyplot as plt
import re
import operator


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

    def ReadRxnScheme (self, fname_json, schemename, verbose = True, FullInfo = False):
        try:
            data = {}
            data[schemename] = json.loads(fname_json)
        except Exception as e:
            try:
                f = open (fname_json)
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

    def RunRxnScheme (self, in_reactants, schemefile_jsontext, schemename, showmols, schemeinfo = None):

        if schemeinfo is not None:
            scheme, rxtants = [schemeinfo[0], schemeinfo [1]]
        else:
            scheme, rxtants = self.ReadRxnScheme(schemefile_jsontext, schemename)
        if scheme is None:
            return 'NOSCHEMEFAIL'
        p=in_reactants [0]
        prod_ct = 1
        intermeds = []
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
            intermeds.append (p)
        return p, prod_ct, [scheme, rxtants, intermeds]

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
        def taskfcn(row, libname, rxschemefile, showstrux, schemeinfo, cycct):
            rxtnts = []
            for ix in range (0, cycct):
                rxtnts.append (row['bb' + str (ix + 1) + '_smiles'])
            try:
                res, prodct, schemeinfo = self.RunRxnScheme(rxtnts, rxschemefile, libname, showstrux, schemeinfo)
            except:
                print('FAIL--Error', rxtnts)
                return 'FAIL'
            if res == 'FAIL':
                print('FAIL in taskfcn')
                print(rxtnts)

            return res

        def rec_bbpull(bdfs, level, cycct, bbslist, ct):
            for i, b in bdfs[level].iterrows():
                bb = b['BB_ID']
                bbs = b['SMILES']
                bbslist.append(bb)
                bbslist.append(bbs)
                if level < cycct - 1:
                    ct = rec_bbpull(bbdfs, level + 1, cycct, bbslist)
                else:
                    print(bbslist)
                    ct += 1
                    if ct % 10000 == 0:
                        print(ct, '/', fullct, end='\r')
                return ct

        infilelist.sort ()
        print (infilelist)
        cycct = len (infilelist)

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

            bdfs = [None] * cycct
            if bblistfile is not None:
                picklistdf = pd.read_csv(bblistfile)
                for ix in range (0, cycct):
                    bdfs [ix] = picklistdf[picklistdf['Cycle'] == 'BB' + str (ix + 1)]
                    bdfs = bdfs.merge(cycdict['BB' + str (ix +1)], on=['BB_ID'], how='left')
            else:
                for ix in range(0, cycct):
                    bdfs[ix] = cycdict['BB' + str(ix+1)]
            print ('BDFS', len (bdfs), cycdict.keys ())

            ct = 0
            fullct = 1
            for ix in range(0, cycct):
                fullct *= len(bdfs[ix])

            print('molecules to enumerate:', rndct)
            if rndct > fullct:
                rndct = -1



            if rndct == -1:
                reslist = [[]] * fullct
                bbslist= []
                rec_bbpull (bdfs, 0, cycct, bbslist, 0)
            else:
                reslist = [[]] * rndct
                ct = 0

                while ct < rndct:
                    bblist = []
                    for ix in range (0, cycct):
                        bbdf = cycdict['BB' + str (ix + 1)].sample()
                        b = bbdf.iloc[0]
                        bb = b['BB_ID']
                        bbs = b['SMILES']
                        bblist.append (bb)
                        bblist.append (bbs)
                    if bblist not in reslist:
                        reslist[ct] = bblist
                        ct += 1
                        if ct % 1000 == 0:
                            print(ct, '/', fullct, end='\r')
            hdrs = []
            for ix in range(0, cycct):
                hdrs.append ('bb' + str (ix + 1))
                hdrs.append('bb' + str(ix + 1) + '_smiles')
            resdf = pd.DataFrame(reslist, columns=hdrs)
        else:
            resdf = infilelist

        print('starting enumeration')
        NUM_WORKERS = 16

        pbar = ProgressBar()
        pbar.register()
        ddf = dd.from_pandas(resdf, npartitions=NUM_WORKERS)
        schemeinfo = self.ReadRxnScheme(rxschemefile, libname, False)
        res = ddf.apply(taskfcn, axis=1, result_type='expand', args=(libname, rxschemefile, rndct == 1, schemeinfo, cycct),
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
        print ('inpath: ', inpath,  libname)
        bbspecprefix = ''
        if bbspec != '' and bbspec is not None:
            bbspecprefix = bbspec + '.'
        flist = pathlib.Path(inpath + libname + '/BBLists').glob(bbspecprefix + '*.csv')
        infilelist = []

        for f in flist:
            c = str(f)
            result = re.search(r'\.BB[0-9]+\.', c)
            if result is not None:
                infilelist.append (c)
        infilelist.sort ()

        if rxschemefile is None:
            rxschemefile = inpath + 'RxnSchemes.json'

        samplespath = inpath + foldername + '/Samples/'
        if not os.path.exists(samplespath):
            os.makedirs(samplespath)
        outpath = inpath + foldername + '/Samples/' + libname + outspec
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

    def TestReactionScheme(self,schemename, rxtnts, rxnschemefile, retIntermeds = False):
        enum = Enumerate()
        res = enum.RunRxnScheme(rxtnts, rxnschemefile, schemename, True)
        if not retIntermeds:
            return res[0]
        else:
            return res[0], res[2][2]

    def Enumerate_Dummy_Scaffold (self, rxnschemefile, schemename, bbsmiles, rxtntnum):
        scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, FullInfo=True)
        if 'scaffold_dummy_structures' not in scheme:
            return'FAIL'
        else:
            inrxtnts = scheme['scaffold_dummy_structures']
            if bbsmiles is not None and rxtntnum is not None:
                inrxtnts [rxtntnum] = bbsmiles
            p, prod_ct, [scheme, rxtants] =  self.RunRxnScheme(inrxtnts, rxnschemefile, schemename, False)
            return p

class EnumerationUI():
    enum = Enumerate()
    rootpath = ''
    initpath = '../CIXTools.init.json'

    def __init__ (self):
        f = open(self.initpath)
        initjson = json.load(f)
        f.close ()
        if 'schemedef' not in st.session_state:
            st.session_state ['scheme_def'] = ''
        self.rootpath = initjson ['lastrootpath']


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
        print ('BBFILES', len(bbfiles), inpath)
        for k in bbfiles.keys ():
            dfs[k] = None
            if 'df' + str(k) not in st.session_state or do_reload:
                if bbfiles[k] != '':
                    dfs [k] = pd.read_csv(bbfiles[k])
                    st.session_state['df' + str(k)] = dfs[k]
            else:
                dfs[k] = st.session_state['df' + str(k)]
        dflist = sorted(dfs.items())
        dflist = [x[1] for x in dflist]
        return dflist

    def SaveToInit(self):
        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        data["lastrootpath"] = st.session_state.rxscheme.replace ('RxnSchemes.json', '')
        if 'schemedef' in st.session_state:
            data["lastschemedef"] = st.session_state.schemedef
        if 'schemepath' in st.session_state:
            data["lastschemepath"] = st.session_state.schemepath
        with open(self.initpath, "w") as jsonFile:
            json.dump(data, jsonFile)

    def SetScheme (self):
        print ('Set Scheme', st.session_state.scheme)
        with open(st.session_state.rxscheme, "r") as jsonFile:
            data = json.load(jsonFile)
            st.session_state.schemedef = json.dumps (data [st.session_state.scheme], indent=4)

    def SaveScheme (self):
        print ('Save Scheme', st.session_state.scheme)
        with open(st.session_state.rxscheme, "r") as jsonFile:
            data = json.load(jsonFile)
        schemejson = json.loads (st.session_state.schemedef)
        data[st.session_state.scheme] = schemejson
        print (data[st.session_state.scheme])
        with open(st.session_state.rxscheme, "w") as jsonFile:
            json.dump(data, jsonFile, indent=4)

    def body(self):
        smilescol='SMILES'
        if 'specstr' not in st.session_state:
            st.session_state['specstr'] = 'Empty'
        rxnschemefile = st.text_input(value=self.rootpath + 'RxnSchemes.json', label='rxscheme', on_change=self.SaveToInit, key='rxscheme')
        if os.path.exists(rxnschemefile) == False:
            st.text (rxnschemefile + ' does not exist. Please adjust path')
            return

        lspath = st.text_input(value=self.rootpath, label='scheme path', on_change=self.SaveToInit, key = 'schemepath')

        ls = lspath.replace(self.rootpath, '')

        f = open(rxnschemefile, 'r')
        schemejson = json.load(f)
        schemelist = []

        if 'schemename' not in st.session_state:
            st.session_state['schemename'] = ls
            bbpath =lspath +  '/BBLists'
        else:
            bbpath = lspath + st.session_state['schemename'] + '/BBLists'

        for s in schemejson:
            schemelist.append(s)
        dfs = self.UpdateBBDfs(bbpath, False)

        with st.expander (label='Scheme Definition Tools'):
            cont1 = st.container ()
        cont2 = st.container()
        cont3= st.container ()


        Enumerate = False
        with cont1:
            ccont = st.container()
            ccontA = st.container ()
            with ccontA:
                newname = st.text_input(label = 'New Scheme Name' )
                if st.button (label = 'Add New Scheme'):
                    st.session_state.schemename= newname
                    ls = newname
                    st.session_state.scheme = newname
                    st.session_state.schemedef = '{}'
                    self.SaveScheme()
                    schemelist.append (newname)

        with cont3:
            colx1, colx2 = st.columns(2)
            with colx1:
                if ls not in schemelist:
                    lsidx = 0
                else:
                    lsidx = schemelist.index(ls)

                schemename = st.selectbox(label='Scheme', options=schemelist, key='scheme', index=lsidx)
                spec = st.text_input(key='spec', label='spec')

                if schemename != st.session_state['schemename'] or spec != st.session_state['specstr']:
                    addspec = ''
                    if spec != '' and spec is not None:
                        addspec = '/' + spec
                    dfs  = self.UpdateBBDfs(self.rootpath + schemename + addspec +'/BBLists', True)
                    self.SetScheme()

                st.session_state['schemename'] = schemename
                st.session_state['specstr'] = spec
                if 'schemedef' not in st.session_state:
                    self.SetScheme()

                for n in range (0, len(dfs)) :
                    if 'bb' + str(n) + 'idx' not in st.session_state:
                        st.session_state['bb' + str(n) + 'idx'] = 0

                rxtnts = [None] * len(dfs)
                for n in range (0, len(dfs)) :
                    df = dfs[n]
                    if df is not None:
                        rxtnts[n] = st.text_input(key='bb' +str(n) + 'txt', label='bb' + str(n) + ' (override)')

                        if rxtnts[n] == '':
                             rxtnts[n] = st.selectbox(label='bb' + str (n), options=dfs[n][smilescol]
                                 , key='bb' + str (n),
                                   index=st.session_state['bb' + str (n) + 'idx'])

        with cont1:
            ccol1, ccol2 = st.columns(2)
            with ccol1:
                if st.button (label='revert'):
                    self.SetScheme()
            with ccol2:
                if st.button (label='save scheme'):
                    self.SaveScheme()
            with ccont:
                st.text_area(height=200, label='Scheme Definition', key='schemedef')

        with cont2:
            col1, col2, col3 = st.columns(3)
            with col1:
                if st.button('random'):
                    for dx in range (0, len(dfs)):
                        if dfs[dx] is not None:
                            st.session_state['bb' + str(dx) + 'idx'] = dfs [dx].index[random.randint(0, len(dfs[dx]))]
                    Enumerate = True

            with col2:
                if st.button('enumerate'):
                    Enumerate = True


            with col3:

                expval = st.button('Export Random Selection')
                countval = st.text_input(label='Count', key='rndcount')
                if expval:
                    try:
                        ct = int (countval)
                    except:
                        ct = 5000

                    addspec = ''
                    if spec != '' and spec is not None:
                        addspec = '/' + spec
                    self.enum.EnumFromBBFiles(schemename, '','', lspath, schemename + addspec, ct, rxnschemefile)

        with cont3:
            if Enumerate == True:
                with colx2:
                    try:
                        res , intermeds= self.enum.TestReactionScheme(schemename, rxtnts, st.session_state.schemedef, True)
                        st.pyplot(MolDisplay.ShowMols(rxtnts))
                        st.pyplot(MolDisplay.ShowMols(intermeds, cols=2, subImgSize=(400,400)))
                        if res is None or res == 'FAIL':
                            fig = plt.figure()
                            st.pyplot(fig)
                        else:
                            st.pyplot(MolDisplay.ShowMol(res))
                    except:
                        st.text ('Error: bad scheme definition')

    def RunUI(self):
        self.head()
        self.body()


if __name__=="__main__":
    if st._is_running_with_streamlit:
        enum = EnumerationUI ()
        enum.RunUI ()