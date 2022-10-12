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
import sys
import gc
import os, psutil
import math
import csv
import time
import copy


class Enumerate:
    rxndict = {}
    chksz = 100000
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
            if r == "product":
                smi = r1
            else:
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


        def rec_bbpull( bdfs, level, cycct, bbslist, ct, reslist, fullct):

            for i, b in bdfs[level].iterrows():
                bb = b['BB_ID']
                bbs = b['SMILES']
                if level == 0:
                    bbslist = []

                bblevellist = copy.deepcopy (bbslist)
                bblevellist.append(bb)
                bblevellist.append(bbs)

                if level < cycct - 1:
                    ct, reslist = rec_bbpull(bdfs, level + 1, cycct, bblevellist, ct, reslist, fullct)
                else:
                    reslist.append(bblevellist)
                    ct += 1
                    if ct % 10000 == 0:
                        print('RECURS CT', ct, '/', fullct, end='\r')

            return ct,  reslist

        if SMILEScolnames is None:
            SMILEScolnames = []
        if BBIDcolnames is None:
            BBIDcolnames = []
        infilelist.sort ()
        cycct = len (infilelist)

        if type(infilelist) is list:
            cycdict = {}
            ct = 1

            cycdict = self.load_BBlists(infilelist, BBIDcolnames, SMILEScolnames)

            bdfs = [None] * cycct
            if bblistfile is not None:
                picklistdf = pd.read_csv(bblistfile)
                for ix in range (0, cycct):
                    bdfs [ix] = picklistdf[picklistdf['Cycle'] == 'BB' + str (ix + 1)]
                    bdfs = bdfs.merge(cycdict['BB' + str (ix +1)], on=['BB_ID'], how='left')
            else:
                for ix in range(0, cycct):
                    bdfs[ix] = cycdict['BB' + str(ix+1)]

            fullct = 1
            for ix in range(0, cycct):
                fullct *= len(bdfs[ix])

            print('molecules to enumerate:', rndct)
            if rndct > fullct:
                rndct = -1
            if int (rndct) == -1:
                reslist= []
                ct ,  reslist = rec_bbpull (bdfs, 0, cycct, [], 0, reslist, fullct)
            else:
                reslist = [[]] * rndct
                ct = 0

                while ct < rndct:
                    bblist = []
                    for ix in range (0, cycct):
                        ri = random.randint (0, len (cycdict['BB' + str (ix + 1)]))
                        # bbdf = cycdict['BB' + str (ix + 1)].sample()
                        b = cycdict['BB' + str (ix + 1)].iloc[ri]
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

            if len (reslist) >  self.chksz:
                with open(outpath + ".EnumList.csv", "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(hdrs, )
                    writer.writerows(reslist)

                enum_in = outpath + ".EnumList.csv"
            else:
                enum_in  = pd.DataFrame (reslist, columns= hdrs)

        else:
            resdf = infilelist
            enum_in = resdf
            hdrs = None

        gc.collect ()
        self.DoParallel_Enumeration (enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct)

    def DoParallel_Enumeration (self, enum_in, hdrs, libname, rxschemefile,outpath, cycct, rndct=-1):
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
        def processchunk (resdf, df):
            pbar = ProgressBar()
            pbar.register()
            ddf = dd.from_pandas(resdf, npartitions=NUM_WORKERS)
            schemeinfo = self.ReadRxnScheme(rxschemefile, libname, False)
            res = ddf.apply(taskfcn, axis=1, result_type='expand',
                            args=(libname, rxschemefile, rndct == 1, schemeinfo, cycct),
                            meta=(0, str)).compute()
            pbar.unregister()

            moddf = resdf.merge(res, left_index=True, right_index=True)
            moddf = moddf.rename(columns={0: 'full_smiles'})

            if outpath is not None:
                moddf[moddf['full_smiles'] != 'FAIL'].to_csv(outpath + '.' + outsuff + '.enum.csv', index=False)
                moddf[moddf['full_smiles'] == 'FAIL'].to_csv(outpath + '.' + outsuff + '.fail.csv', index=False)
                moddf.to_csv(outpath + '.' + outsuff + '.all.csv', mode='a', index=False, header=False)
                gc.collect()
                return None
            else:
                if df is None:
                    df = moddf
                else:
                    df.append(moddf, ignore_index=True)
                gc.collect()
                return df

        print('starting enumeration')
        NUM_WORKERS = 16

        df  = None
        outsuff = 'full'
        if rndct != -1:
            outsuff = str(rndct)
        f = open(outpath + '.' + outsuff + '.all.csv', 'w')
        f.close()

        if type (enum_in) is str:
            reader = pd.read_csv(enum_in, chunksize=self.chksz)
            df = None
            for chunk in reader:
                resdf = pd.DataFrame (chunk, columns = hdrs)
                df = processchunk(resdf, df)
        else:
            df = None
            df = processchunk(enum_in, df)

        if outpath is not None:
            return outpath + '.' + outsuff + '.all.csv'
        else:
            return df

    def EnumFromBBFiles(self, libname, bbspec, outspec, inpath, foldername, num_strux, rxschemefile, picklistfile=None, SMILEScolnames = [], BBcolnames = []):
        print ('inpath: ', inpath,  libname,  num_strux)
        bbspecprefix = ''
        outspecprefix = ''
        if bbspec != '' and bbspec is not None:
            outspecprefix = '/' + outspec
        bbpath = inpath + libname + outspecprefix + '/BBLists'
        flist = pathlib.Path(bbpath).glob('*.' + bbspec + '*.csv')
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
        outpath = inpath + foldername + '/Samples/' + libname
        if  outspec != '' and outspec is not None:
            outpath += '.' + outspec
        outfile = self.enumerate_library_strux(libname, rxschemefile, infilelist, outpath, num_strux, picklistfile, SMILEScolnames=SMILEScolnames, BBIDcolnames=BBcolnames)
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
            cdf = df[df[ln] == 'Y'].drop_duplicates(subset=['BB_ID']).reset_index(drop=True)
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
            if SMILEScolnames is not None and len (SMILEScolnames) > 0:
                for smc in SMILEScolnames:
                    if smc in bbdict[cyc].columns:
                        changecoldict[smc] = 'SMILES'
            if BBIDcolnames is not None and len (BBIDcolnames) > 0:
                for bbidc in BBIDcolnames:
                    if bbidc in bbdict[cyc].columns:
                        changecoldict[bbidc] = 'BB_ID'
            if len(changecoldict) > 0:
                bbdict[cyc] = bbdict[cyc].rename(columns=changecoldict)
            bbdict[cyc] = bbdict[cyc].drop_duplicates(subset='BB_ID', keep="first").reset_index(drop=True)
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
        res = enum.RunRxnScheme(rxtnts, rxnschemefile, schemename, False)

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
    initpath = '../CIxTools.init.json'
    smiles_colnames = None
    bbid_colnames = None

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
        for k in st.session_state.keys ():
            if k.startswith( 'bb') and k.endswith('idx'):
                st.session_state[k] = 0
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
            if 'specstr' in st.session_state:
                specstr = '/' + st.session_state['specstr']
            else:
                specstr = ''
            bbpath = lspath + st.session_state['schemename'] + specstr +  '/BBLists'

        for s in schemejson:
            schemelist.append(s)
        schemelist.sort ()

        dfs = self.UpdateBBDfs(bbpath, False)
        print ('LEN DFS', len (dfs))

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
                        rxtnts[n] = st.text_input(key='bb' +str(n) + 'txt', label='BB' + str(n+1) + ' (override)')
                        if rxtnts[n] == '':
                             rxtnts[n] = st.selectbox(label='BB' + str (n + 1), options=dfs[n][smilescol]
                                 , key='bb' + str (n ),
                                   index=int(st.session_state['bb' + str(n) + 'idx']))
                    print ('INDEX', st.session_state['bb' + str(n) + 'idx'])
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
                countval = st.text_input(label='Count', key='rndcount', value = 5000 )
                if expval:
                    try:
                        ct = int (countval)
                    except:
                        ct = 5000

                    addspec = ''
                    if spec != '' and spec is not None:
                        addspec = '/' + spec
                    self.enum.EnumFromBBFiles(schemename, spec, spec, lspath, schemename + addspec, ct, rxnschemefile, SMILEScolnames=self.smiles_colnames, BBcolnames=self.bbid_colnames)

        with cont3:
            if Enumerate == True:
                with colx2:
                    if rxtnts is None or len(rxtnts) == 0:
                        st.text('Reactants not correctly loaded')
                    else:
                        try:
                            res , intermeds= self.enum.TestReactionScheme(schemename, rxtnts, st.session_state.schemedef, True)
                            if res is None or res == 'FAIL':
                                st.text ('Reaction Failure')
                            else:
                                st.pyplot(MolDisplay.ShowMol(res))
                            with st.expander(label='Reaction Info'):
                                st.pyplot(MolDisplay.ShowMols(rxtnts))
                                st.pyplot(MolDisplay.ShowMols(intermeds, cols=2, subImgSize=(400,400)))

                        except Exception as e:
                            st.text ('Error: bad scheme definition')
                            exc_type, exc_obj, exc_tb = sys.exc_info()
                            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                            st.text (str(exc_type) + ' ' +  fname + ' ' + str( exc_tb.tb_lineno))
                            st.text (e)

    def RunUI(self):
        self.head()
        self.body()


if __name__=="__main__":
    if st._is_running_with_streamlit:
        enum = EnumerationUI ()
        enum.RunUI ()