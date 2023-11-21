from rdkit import Chem
from rdkit.Chem import rdChemReactions
import json
import toml
import pandas as pd
from numpy import random
from multiprocessing.pool import ThreadPool as Pool
import threading
import MolDisplay
import ChemUtilities
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from rdkit.Chem import SaltRemover
import pathlib
from tqdm import tqdm
import re
import gc
import os, psutil
import csv
import time
import copy
import argparse
import yaml
import numexpr

CPU_COUNT = os.cpu_count()
NUM_WORKERS = CPU_COUNT * 2
chksz = 50000
numexpr.set_num_threads(numexpr.detect_number_of_cores())
if "NUMEXPR_MAX_THREADS" not in os.environ:
    os.environ["NUMEXPR_MAX_THREADS"] = '16'
else:
    NUM_WORKERS = int(os.environ["NUMEXPR_MAX_THREADS"] )


class Enumerate:
    smilescol = '__SMILES'
    idcol = '__BB_ID'
    rxndict = {}

    named_reactions = None
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
        prod_ct = 0
        for r in SM_rxn:
            if r == "product":
                smi = r1
            else:
                rxn = rdChemReactions.ReactionFromSmarts(r)
                products = rxn.RunReactants(reacts)
                if len(products) == 0:
                    product = None
                else:
                    ct = 0
                    if len(products) > 1:

                        proddict= {}

                        for cx in range (0, len (products)):
                            if Chem.MolToSmiles(products[cx][0]) not in proddict:
                                proddict [Chem.MolToSmiles(products[cx][0])] = 1
                                prod_ct += 1
                            else:
                                proddict[Chem.MolToSmiles(products[cx][0])] += 1

                    product = products [0][0]

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

        return smi, prod_ct, products

    def Show_Reactants (self, reactants):
        for r in reactants:
            MolDisplay.ShowMol (r)

    def ReadRxnScheme (self, fname_json, schemename, verbose = True, FullInfo = False):
        try: #case fname_json is json
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
                    if rxn == 'C' or rxn == '[C;H4]':
                        print ('HERE2')
                        if reactants[0] == 'C' or reactants [1] == 'C':
                            reaction = step["Rxns"][rxn]
                            print ('HERE', reaction)
                            break
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
            print ('ENDRXN:', reaction)
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
            return 'NOSCHEMEFAIL',0, None
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
                if self.named_reactions is not None:
                    rlist = []
                    if type (reaction) == str:
                        reaction = [reaction]
                    for r in reaction:
                        if r in self.named_reactions:
                            if type (self.named_reactions [r] ) == str:
                                rlist.append (self.named_reactions [r])
                            else:
                                rlist.extend (self.named_reactions [r])
                        else:
                            rlist.append (r)
                    reaction = rlist
                try:
                    p, outct, products = self.React_Molecules(reactants [0], reactants [1],  reaction,  showmols)
                except:
                    p = None
                    outct = 0
                if outct > 0:
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
                rvals.append (row [self.smilescol])
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
                rvals.append (rcx[self.smilescol])
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

        scheme, reactants= self.ReadRxnScheme(schemepath, scheme_name)
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
            pool_size = NUM_WORKERS
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

    def enumerate_library_strux(self, libname, rxschemefile, infilelist, outpath, rndct=-1, bblistfile=None, SMILEScolnames = [], BBIDcolnames = [], removeduplicateproducts = False, outtype = 'filepath', write_fails_enums = True, retIntermeds=False):
        def rec_bbpull( bdfs, level, cycct, bbslist, ct, reslist, fullct, hdrs, currct = 0, appendmode = False, retIntermeds = False):
            if reslist is None:
                reslist = [[]] * min(chksz, fullct)
                print ('RESLIST', len(reslist))

            for i, b in bdfs[level].iterrows():
                bb = b[self.idcol]
                bbs = b[self.smilescol]
                if level == 0:
                    bbslist = []

                bblevellist = copy.deepcopy (bbslist)
                bblevellist.append(bb)
                bblevellist.append(bbs)

                if level < cycct - 1:
                    ct, reslist, currct, appendmode = rec_bbpull(bdfs, level + 1, cycct, bblevellist, ct, reslist, fullct, hdrs, currct = currct, appendmode=appendmode, retIntermeds=retIntermeds)
                else:

                    reslist[currct] = bblevellist

                    ct += 1
                    currct += 1
                    if ct % chksz == 0:
                        print('RECURS CT', ct, '/', fullct, ' currct:', currct, end='\r')

                    if currct == chksz or ct == fullct:
                        enum_in = pd.DataFrame (reslist, columns=hdrs)
                        self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct,
                                                    removeduplicateproducts, appendmode=appendmode, retIntermeds=retIntermeds)
                        reslist = [[]] * min(chksz, fullct - ct)
                        currct = 0
                        appendmode = True
                        gc.collect ()


            return ct, reslist,  currct, appendmode

        if SMILEScolnames is None:
            SMILEScolnames = []
        if BBIDcolnames is None:
            BBIDcolnames = []

        infilelist.sort ()
        cycct = len (infilelist)

        if type(infilelist) is list:
            cycdict = self.load_BBlists(infilelist, BBIDcolnames, SMILEScolnames)

            bdfs = [None] * cycct
            if bblistfile is not None:
                picklistdf = pd.read_csv(bblistfile)
                for ix in range (0, cycct):
                    bdfs [ix] = picklistdf[picklistdf['Cycle'] == 'BB' + str (ix + 1)]
                    bdfs = bdfs.merge(cycdict['BB' + str (ix +1)], on=[self.idcol], how='left')
            else:
                for ix in range(0, cycct):
                    bdfs[ix] = cycdict['BB' + str(ix+1)]

            fullct = 1
            for ix in range(0, cycct):
                fullct *= len(bdfs[ix])

            print('molecules to enumerate:', rndct)
            if rndct > fullct:
                rndct = -1

            hdrs = []
            for ix in range(0, cycct):
                hdrs.append('bb' + str(ix + 1))
                hdrs.append('bb' + str(ix + 1) + '_smiles')

            appendmode = False

            if rndct == -1 :
                hdrs = []
                for ix in range(0, cycct):
                    hdrs.append('bb' + str(ix + 1))
                    hdrs.append('bb' + str(ix + 1) + '_smiles')
                with open(outpath + ".EnumList.csv", "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(hdrs)
                rec_bbpull (bdfs, 0, cycct, [], 0, None, fullct, hdrs, appendmode = appendmode, retIntermeds=retIntermeds)
            else:
                reslist = [[]] * min(rndct, chksz)
                ct = 0
                currct = 0


                while ct < rndct:
                    bblist = []
                    for ix in range (0, cycct):
                        ri = random.randint (0, len (cycdict['BB' + str (ix + 1)]))
                        b = cycdict['BB' + str (ix + 1)].iloc[ri]
                        bb = b[self.idcol]
                        bbs = b[self.smilescol]
                        bblist.append (bb)
                        bblist.append (bbs)
                    if bblist not in reslist:
                        reslist[currct] = bblist
                        ct += 1
                        currct += 1
                        if ct % 1000 == 0:
                            print(ct, '/', fullct, end='\r')
                        if currct ==  chksz  or ct == rndct:
                            enum_in = pd.DataFrame(reslist, columns=hdrs)
                            outpath = self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct, removeduplicateproducts, appendmode=appendmode, retIntermeds=retIntermeds)
                            reslist = [[]] * min (chksz, rndct - ct)
                            currct = 0
                            appendmode = True

        else:
            resdf = infilelist
            enum_in = resdf
            hdrs = None
            outpath = self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct,
                                        removeduplicateproducts, appendmode = False, write_fails_enums=write_fails_enums, retIntermeds=retIntermeds)

        return outpath

    def DoParallel_Enumeration (self, enum_in, hdrs, libname, rxschemefile,outpath, cycct, rndct=-1, removeduplicateproducts = False, appendmode = False, write_fails_enums=True, retIntermeds = False):

        def taskfcn(row, libname, rxschemefile, showstrux, schemeinfo, cycct, retIntermeds = False):
            rxtnts = []
            for ix in range(0, cycct):
                rxtnts.append(row['bb' + str(ix + 1) + '_smiles'])
            try:
                res, prodct, schemeinfo = self.RunRxnScheme(rxtnts, rxschemefile, libname, showmols=showstrux, schemeinfo=schemeinfo)
                if prodct > 1:
                    print ('FAIL--MULTIPLE PRODUCTS', res)
                    res =  ['FAIL--MULTIPLE PRODUCTS']
                else:
                    res = [res]
            except Exception as e:
                res =  ['FAIL']

            if retIntermeds:
                for k in range (0, len(schemeinfo[0])):
                    res.append(schemeinfo[0][k])
                    if  k < len(schemeinfo[2]):
                        res.append (schemeinfo[2][k])
                    else:
                        res.append (None)
                return res
            return res

        def processchunk (resdf, df, outpath, retIntermeds, schemeinfo):
            pbar = ProgressBar()
            pbar.register()
            ddf = dd.from_pandas(resdf, npartitions=CPU_COUNT * 10)
            if len (ddf) <  5000:
                wrkrs = 1
            else:
                wrkrs  = NUM_WORKERS
            metaser = [(0, object)]
            colnames = {0: 'full_smiles'}
            if retIntermeds:
                rxnct = len(schemeinfo [0])
                for r in range(0,rxnct):
                    metaser.append ((2*r + 1, str))
                    metaser.append((2 * r + 2, str))
                    colnames [2*r + 1] = 'step' + str(r) + '_rxn'
                    colnames [2*r + 2] = 'step' + str(r) + '_intermed'
            res = ddf.apply(taskfcn, axis=1, result_type='expand',
                            args=(libname, rxschemefile, rndct == 1, schemeinfo, cycct, retIntermeds),
                            meta=metaser).compute(scheduler='processes',  num_workers=wrkrs)
            pbar.unregister()
            moddf = resdf.merge(res, left_index=True, right_index=True)
            moddf = moddf.rename(columns=colnames)

            if outpath is not None:
                enumdf = moddf[~moddf['full_smiles'].isin( ['FAIL','FAIL--MULTIPLE PRODUCTS'])]
                enumdf.to_csv(outpath + '.' + outsuff + '.enum.csv', mode='a', index=False, header=False)
                faildf = moddf[moddf['full_smiles'].isin(['FAIL','FAIL--MULTIPLE PRODUCTS'])]
                if write_fails_enums:
                    faildf.to_csv(outpath + '.' + outsuff + '.fail.csv', mode='a', index=False, header=False)
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
        print('\nstarting enumeration')

        df  = None
        outsuff = 'full'
        if rndct != -1:
            outsuff = str(rndct)
        hdrstr = ','.join (hdrs)
        schemeinfo = self.ReadRxnScheme(rxschemefile, libname, False)
        if (schemeinfo is None):
            print ('Failed, scheme not detected')
            raise(13)
        if outpath is not None:
            if not write_fails_enums:
                flist = [outpath + '.' + outsuff + '.all.csv',None, None]
            else:
                flist = [outpath + '.' + outsuff + '.all.csv', outpath + '.' + outsuff + '.fail.csv', outpath + '.' + outsuff + '.enum.csv']
            if appendmode == False:
                for fname in flist:
                    f = open(fname, 'w')
                    f.write(hdrstr +',full_smiles')
                    rxnct = len(schemeinfo[0])
                    if retIntermeds:
                        for r in range(0, rxnct):
                            f.write( ',step' + str(r) + '_rxn')
                            f.write(',step' + str(r) + '_intermed')
                    f.write ('\n')
                    f.close()

        if type (enum_in) is str:
            reader = pd.read_csv(enum_in, chunksize=chksz)
            df = None
            cct = 0
            for chunk in reader:
                resdf = pd.DataFrame (chunk, columns = hdrs)
                print('CHUNK', cct)
                df = processchunk(resdf, df, outpath , schemeinfo=schemeinfo, retIntermeds=retIntermeds)
                cct += 1
        else:
            df = None
            df = processchunk(enum_in, df, outpath, schemeinfo=schemeinfo,retIntermeds=retIntermeds)

        if outpath is not None:
            path = outpath + '.' + outsuff + '.all.csv'
            if removeduplicateproducts:
                path = self.Deduplicate (outpath, outsuff)
            return path
        else:
            if removeduplicateproducts:
                df = df.drop_duplicates(keep='first', subset = ['full_smiles'])
            return df

    def Deduplicate (self, outpath, outsuff):
        df = pd.read_csv(outpath + '.' + outsuff + '.all.csv')
        df = df.drop_duplicates(keep='first', subset=['full_smiles'])
        df.to_csv(outpath + '.' + outsuff + '.all.dedup.csv')
        return outpath + '.' + outsuff + '.all.dedup.csv'

    def Get_BBFiles (self, bb_specialcode, lib_subfolder, enumsfolder, libname):
        libspecprefix = ''
        if lib_subfolder != '' and lib_subfolder is not None:
            libspecprefix = '/' + lib_subfolder
        bbpath = enumsfolder + libname + libspecprefix + '/BBLists'
        if bb_specialcode is not None and bb_specialcode != '':
            srchstr = libname + '.' + bb_specialcode +  '.BB?.csv'
            print (bbpath)
            flist = pathlib.Path(bbpath).glob(srchstr)
        else:
            srchstr = libname + '.BB?.csv'
            flist = pathlib.Path(bbpath).glob(srchstr)
        infilelist = []

        for f in flist:
            c = str(f)
            result = re.search(r'\.BB[0-9]+\.', c)
            if result is not None:
                infilelist.append(c)
        infilelist.sort()
        if len(infilelist) == 0:
            print('FAIL: No BB files found with format ' + srchstr + ' found in ' + bbpath, infilelist)
            print ('Inputs:BBcode:' , bb_specialcode, ' lib_subfldr:', lib_subfolder, ' enum folder:', enumsfolder, ' lib:', libname)
            return ''
        else:
            print (infilelist)
            print (srchstr)
        return infilelist

    def EnumFromBBFiles(self, libname, bb_specialcode, out_specialcode, enumsfolder, lib_subfolder, num_strux, rxschemefile, picklistfile=None,
                        SMILEScolnames = [], BBcolnames = [], rem_dups = False, returndf = False, write_fails_enums = True, retIntermeds = False):

        infilelist = self.Get_BBFiles (bb_specialcode, lib_subfolder, enumsfolder, libname)
        if type (infilelist) == str:
            return

        if rxschemefile is None:
            rxschemefile = enumsfolder + 'RxnSchemes.json'

        samplespath = enumsfolder  + libname + '/' + lib_subfolder + '/Samples/'
        print ('SAMPLESPATH', samplespath)
        if not os.path.exists(samplespath):
            os.makedirs(samplespath)
        outpath = samplespath + libname
        if  out_specialcode != '' and out_specialcode is not None:
            out_specialcode += '.' + out_specialcode

        if returndf is True:
            outpath = None
        outfile = self.enumerate_library_strux(libname, rxschemefile, infilelist, outpath, num_strux, picklistfile,
                                               SMILEScolnames=SMILEScolnames, BBIDcolnames=BBcolnames, removeduplicateproducts=rem_dups, write_fails_enums=write_fails_enums, retIntermeds = retIntermeds)
        return outfile

    def FilterBBs(self, bbdict, filterfile):
        filtersdf = pd.read_csv(filterfile)
        patterndict = {}
        removeidxs = {}
        for ix, r in filtersdf.iterrows():
            pattern = r['smarts']
            patterndict[pattern] = Chem.MolFromSmarts(pattern)

        for dx in bbdict.keys():
            removeidxs[dx] = []
            for idx, row in bbdict[dx].iterrows():
                try:
                    m = Chem.MolFromSmiles(row[self.smilescol])
                except Exception as e:
                    print (filterfile)
                    print (row [self.idcol], row[self.smilescol])
                    raise (e)
                for v in patterndict.values():
                    if m.HasSubstructMatch(v) == True:
                        removeidxs[dx].append(idx)
                        continue
            bbdict[dx].drop(removeidxs[dx], axis=0, inplace=True)
            bbdict[dx] = bbdict[dx].reset_index ()
        return bbdict

    def Get_MoleculesFromSMARTSFilters (self, infile, outfile, inclSMARTS, exclSMARTS ):
        df = pd.read_csv (infile)
        keep_idxs = []
        filters_dict = {}
        filtername = 'f1'
        filters_dict [filtername] ={}
        if inclSMARTS is None:
            inclSMARTS = []
        filters_dict[filtername]['include'] = inclSMARTS
        if exclSMARTS is None:
            exclSMARTS = []
        filters_dict[filtername]['exclude'] = exclSMARTS
        for idx, row in df.iterrows ():
            input_mol = Chem.MolFromSMILES (row [self.smilescol])
            if (self.PassFilters(input_mol, filters_dict, filtername)):
                keep_idxs.append(idx)
        df = df.iloc [keep_idxs].reset_index ()
        print (df)
        df.to_csv (outfile)

    def PassFilters(self, input_mol, filters_dict, filtname):
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

    def generate_BBlists(self, inpath, outpath, libname, rxnschemefile):
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
                input_mol = Chem.MolFromSmiles(row[self.smilescol])
                input_mol = saltstrip.StripMol(input_mol)
                df.iloc[idx][self.smilescol] = Chem.MolToSmiles(input_mol)
            except:
                input_mol = None
            if input_mol is not None:
                for ln in listnames:
                    if self.PassFilters(input_mol, filters_dict, ln):
                        df.iloc[idx][ln] = 'Y'
                    else:
                        df.iloc[idx][ln] = 'N'
        print('loop done')
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        for ln in listnames:
            cdf = df[df[ln] == 'Y'].drop_duplicates(subset=[self.idcol]).reset_index(drop=True)
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
                    if '.' not in l:
                        bbdict [cyc] = None
                    else:
                        bbdict[cyc] = pd.read_csv(l)

        for cyc in bbdict.keys():
            if bbdict [cyc] is not None:
                changecoldict = {}
                if SMILEScolnames is not None and len (SMILEScolnames) > 0:
                    for smc in SMILEScolnames:
                        if smc in bbdict[cyc].columns:
                            changecoldict[smc] = self.smilescol
                            break
                if BBIDcolnames is not None and len (BBIDcolnames) > 0:
                    for bbidc in BBIDcolnames:
                        if bbidc in bbdict[cyc].columns:
                            changecoldict[bbidc] = self.idcol
                            break
                if len(changecoldict) > 0:
                    bbdict[cyc] = bbdict[cyc].rename(columns=changecoldict)
                bbdict[cyc] = bbdict[cyc].drop_duplicates(subset=self.idcol, keep="first").reset_index(drop=True)
        return bbdict

    def pull_BBs(self, inpath, idcol, smilescol):
        globlist = pathlib.Path(inpath).glob('*.csv')
        full_df = None
        for f in tqdm(globlist):
            if full_df is None:
                full_df = pd.read_csv(f)
                full_df = full_df[[idcol, smilescol]]
                full_df = full_df.rename(columns={idcol: self.idcol, smilescol: self.smilescol})
            else:
                df = pd.read_csv(f)
                df = df.rename(columns={idcol: self.idcol, smilescol: self.smilescol})
                try:
                    df = df[[self.idcol, self.smilescol]].dropna()
                except:
                    print('exiting:', f)
                    exit()
                full_df = full_df.append(df, ignore_index=True)
        return full_df

    def TestReactionScheme(self,schemename, rxtnts, rxnschemefile, retIntermeds = False):
        res = self.RunRxnScheme(rxtnts, rxnschemefile, schemename, False)
        rxnslist = []
        for k in res [2][0]:
            print ('X:', k)
            rxnslist.append([str(k['Reactants']) + '<p>' + str(k['Rxns'])])
        if res[1] > 1:
            if not retIntermeds:
                return 'FAIL--MULTIPLE PRODUCTS'
            else:
                return 'FAIL--MULTIPLE PRODUCTS', None, None
        if not retIntermeds:
            return res[0]
        else:
            if res [2] is not None:
                return res[0], res[2][2], rxnslist
            else:
                return res[0], None, None

    def Enumerate_Dummy_Scaffold (self, rxnschemefile, schemename, bbsmiles, rxtntnum):
        scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, FullInfo=True)
        if 'scaffold_dummy_structures' not in scheme:
            return'FAIL'
        else:
            inrxtnts = scheme['scaffold_dummy_structures']
            if bbsmiles is not None and rxtntnum is not None:
                inrxtnts [rxtntnum] = bbsmiles
            p, prod_ct, [scheme, rxtants, intermeds] =  self.RunRxnScheme(inrxtnts, rxnschemefile, schemename, False)
            return p

    def Get_CycleList (self, rxnschemefile, schemename):
        scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, verbose=True, FullInfo=True)
        return scheme ['BB_Cycles']

    def Get_LibCycle_BBCriteria (self, rxnschemefile, schemename, cycle, fieldname='filters'):
        scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename,  FullInfo=True)
        if fieldname in scheme:
            if cycle in scheme [fieldname]:
                return scheme [fieldname][cycle]
            else:
                if 'names' in scheme [fieldname]:
                    for n in scheme [fieldname]['names']:
                        if scheme [fieldname]['names'][n] == cycle:
                            return scheme [fieldname]['BB_filters'][n]
                    return None
        else:
            return None

class EnumerationCLI :
    @staticmethod
    def Run_CLI (SMILEScolnames = None, BBcolnames = None):
        paramdefaults = [ ('rxnschemefile', './RxnSchemes.json'), ('schemepath','.'), ('scheme',''), ('schemespec',''), ('numstrux', 5000), ('removedups', False)]
        parser = argparse.ArgumentParser(description='Enumeration Options')
        parser.add_argument('-p', '--paramfile', nargs='?', default=None, type=str,
                            help='optional .yaml file for commandline paramaters')
        parser.add_argument('-r', '--rxnschemefile', nargs='?', default=None, type=str,
                            help='Rxnschemes.json file path')
        parser.add_argument('-sp', '--schemepath', nargs='?', default=None, type=str, help='Enumerations folder path')
        parser.add_argument('-s', '--scheme', nargs='?', default=None, type=str, help='Scheme Name')
        parser.add_argument('-sx', '--schemespec', nargs='?', default=None, type=str, help='sub-scheme Specifier')
        parser.add_argument('-bx', '--bbspecialcode', nargs='?', default=None, type=str, help='BB files special code')
        parser.add_argument('-n', '--numstrux', nargs='?', default=None, type=int,
                            help='number of structures to enumerate (-1 for all)')
        parser.add_argument('-rd', '--removedups', nargs='?', default=None, type=str,
                            help='Remove duplicate structures (True/False)')
        parser.add_argument('-wfe', '--write_fails_enums', nargs='?', default=None, type=str,
                            help='Write Fails and Enumerated Molecules in separate files (True/False)')
        parser.add_argument('-rxntest', '--rxntest', nargs=2, default=None, type=str,
                            help='2 reactants SMILES')
        parser.add_argument('-rxn', '--rxn', nargs=1, default=None, type=str,
                            help='Reaction SMARTS')
        args = vars(parser.parse_args())

        if args['paramfile'] is not None:
            with open(args['paramfile'], 'r') as stream:
                try:
                    prms = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
            for k in prms.keys ():
                if k not in args or args [k] is None or args[k] == '':
                    args [k] = prms[k]
            for x in paramdefaults:
                if x[0] not in args or args [x[0]] is None or args[x[0]] == '':
                    args [x[0]] = x[1]
        specstr = args['schemespec']
        addspec = ''
        if specstr != '' and specstr is not None:
            addspec = '/' + specstr

        if 'rxntest' in args:
            enum = Enumerate()
            res = enum.React_Molecules(args ['rxntest'] [0],args ['rxntest'] [1] , args['rxn'], True)
            print (res[0])
            exit ()

        if 'removedups' in args:
            if args['removedups'] not in [True, False]:
                if args['removedups'] == 'True':
                    rd = True
                elif args['removedups'] == 'False':
                    rd = False
                else:
                    'Fail: rd must be True or False'
                    exit()
            else:
                rd = args["removedups"]

        if 'write_fails_enums' in args:
            if args['write_fails_enums'] not in [True, False]:
                if args['write_fails_enums'] == 'True':
                    write_fails_enums = True
                elif args['write_fails_enums'] == 'False':
                    write_fails_enums = False
                else:
                    print ('Fail: write_fails_enums must be True or False')
                    exit()
            else:
                write_fails_enums = args["write_fails_enums"]

        enum = Enumerate()
        print ('Starting Enumeration')
        outf = enum.EnumFromBBFiles(args['scheme'], args['bbspecialcode'], args['schemespec'], args['schemepath'],
                             addspec, args['numstrux'],
                             args['rxnschemefile'], SMILEScolnames=SMILEScolnames, BBcolnames=BBcolnames, rem_dups=rd, write_fails_enums=write_fails_enums)
        print ('output prefix', outf)
        print('End Enumeration')

if __name__=="__main__":
    EnumerationCLI.Run_CLI()




