import pandas as pd
from rdkit.Chem import AllChem as Chem
import rdkit.Chem.Scaffolds.MurckoScaffold as MScaff
import rdkit.DataStructs as DataStructs
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import numpy as np
import ChemUtilities
import random
from tqdm import tqdm
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import Enumerate
from itertools import islice
from multiprocessing.pool import ThreadPool as Pool
import operator
import math
import streamlit as st
import json
import copy
import os
import numexpr
import time
import NamingStandardizer as ns
import sys
from st_aggrid import AgGrid
from st_aggrid.grid_options_builder import  GridOptionsBuilder
from st_aggrid import GridUpdateMode, DataReturnMode
import MolDisplay
from streamlit import runtime

CPU_COUNT = os.cpu_count()
NUM_WORKERS = CPU_COUNT * 2
if "NUMEXPR_MAX_THREADS" not in os.environ:
    os.environ["NUMEXPR_MAX_THREADS"] = '16'
numexpr.set_num_threads(numexpr.detect_number_of_cores())

class Diversity:
    enum = Enumerate.Enumerate()
    def DiversitySelection (self, infile, out_file, n_analyzed, ss_file, deprotect_specfile,
                            exclude_file, n_cmpds, keeporder, priority_col, priority_val, propfilters):
        print (infile)
        infile = open (infile)
        df = pd.read_csv (infile)
        infile.close ()
        #randomize
        if not keeporder:
            df = df.sample(frac=1)

        if 'SMILES' in df.columns:
            smilescol = df.columns.get_loc ('SMILES')
        else:
            smilescol = df.columns.get_loc('smiles')
        if 'VENDOR_ID' in df.columns:
            idcol = df.columns.get_loc('VENDOR_ID')
        else:
            idcol = df.columns.get_loc('vendor_id')
        df = df.assign(ORIGSMILES='')
        origsmilescol = df.columns.get_loc('ORIGSMILES')

        scaffdict, incmpd_list = self.filters_scaffolds (df, n_analyzed, idcol,smilescol,origsmilescol,
                                                         deprotect_specfile, ss_file, exclude_file, propfilters)

        if n_cmpds > len(list(scaffdict.keys())):
            n_cmpds = len(list(scaffdict.keys())) -2

        scaff_clusters, medoids, n_cmpds = self.Kmeans (list(scaffdict.keys()), n_cmpds)

        f = open(out_file, 'w')
        f.write(','.join(df.columns))  # id,SMILES
        f.write(',cluster,scaffold,orig_smiles\n')

        out_list = self.choose_compounds (scaffdict, scaff_clusters, medoids, f, smilescol, origsmilescol, priority_col, priority_val)

        self.Generate_UMAP(incmpd_list, out_list, out_file.replace('.csv', '.png'), out_file.replace('.csv', '.umap.csv'))
        f.close()
        return
    def filters_scaffolds (self, df, n_analyzed, idcol, smilescol,origsmilescol,
                           deprotect_specfile, ss_file,  exclude_file, propfilters):
        scaffdict = {}
        incmpd_list = []
        excllist = []

        out_deprot_file = open('Archive/deprotlog.csv', 'w')
        out_deprot_file.write('SMILES,deprot_smiles,deprot_res_ct\n')
        out_filter_file = open('Archive/filterlog.csv', 'w')
        out_filter_file.write('SMILES,SS_found\n')

        if exclude_file is not None:
            exclude_file = open(exclude_file)
            excldf = pd.read_csv(exclude_file)
            exclude_file.close()
            excldf.columns = excldf.columns.str.lower()
            excllist = excldf['vendor_id'].tolist()

        del_rows = []
        tct = 0
        filterct = 0
        uct = 0
        errct = 0
        ssfilters = ChemUtilities.Read_FilterFile(ss_file)

        for index, row in df.iterrows() :
            if origsmilescol == -1:
                origsmilescol = len (row)
            if not (n_analyzed == -1 or tct < n_analyzed):
                break

            if row [idcol] in excllist:
                print ('excluded')
                continue

            m = Chem.MolFromSmiles(row[smilescol])
            if deprotect_specfile is not None:
                try:
                    res, deprotected = ChemUtilities.Deprotect (m, deprotect_specfile)
                except:
                    deprotected = False
                    res = [m]
            else:
                deprotected = False

            rowcopy = None
            rowlist = []

            if deprotected:
                out_deprot_file.write(row[smilescol] + ',')
                if len (res) ==1:
                    row [origsmilescol] = row[smilescol]
                    row [smilescol] = Chem.MolToSmiles(res [0])
                    rowlist.append ([row, res[0]])
                    out_deprot_file.write(row[smilescol] + ',1\n')
                else:
                    orig = row[smilescol]
                    row[smilescol] = Chem.MolToSmiles(res[0])
                    rowlist.append([row, res[0]])
                    for ixr in range (1,len(res)):
                        rowcopy = row.copy()
                        rowcopy[smilescol]  = Chem.MolToSmiles(res [ixr])
                        rowcopy [origsmilescol] = orig
                        rowlist.append ([rowcopy, res[ixr]])
                    out_deprot_file.write(rowcopy[smilescol] + ',' + str (ixr + 1) +  '\n')
            else:
                rowlist.append ([row, m])


            if rowcopy is None:
                rows = [row]
            else:
                rows = [row, rowcopy]

            for rowinfo in rowlist:
                row = rowinfo [0]
                m = rowinfo [1]
                try:
                    fail_filters = ChemUtilities.Filters(m,propfilters )
                    if fail_filters:
                        del_rows.append (index)
                    else:
                        has_filtered_ss= ChemUtilities.Substructure_Filters(m, out_filter_file, ssfilters, False, True)
                        if has_filtered_ss:
                            del_rows.append (index)

                    if not fail_filters and not has_filtered_ss:
                        scaffdict, addnum = self.CheckListScaffoldMatch (row, scaffdict, m)
                        incmpd_list.append (row [smilescol])
                    else:
                        addnum = -1
                except Exception as e:
                    print ('divsel')
                    print ('exception:' + str(e) )
                    errct += 1
                    addnum = 0
                if addnum == -1:
                    filterct += 1
                else:
                    uct += 1
            if (tct%1000 == 0):
                print (filterct, uct, tct, errct)
            tct += 1

        out_filter_file.close()
        out_deprot_file.close()
        return scaffdict, incmpd_list
    def choose_compounds (self, scaffdict, scaff_clusters, medoids, out_file, smilescol, origsmilescol,
                          priority_col, priority_val):

        out_fps = []
        out_list = []
        sim_cutoff = .5

        cluster_num = 0

        for idx in medoids:
            key = list(scaffdict.keys())[idx]
            curr_cluster = scaff_clusters[key]
            cscaffs = [k for k,v in scaff_clusters.items() if v ==curr_cluster]
            cmpds = scaffdict[key]
            cdf = pd.DataFrame(cmpds)
        #use medoid unless priority col is specified, prefer medoid priority, then priority, then non-priority medoid compounds
            if (priority_col is not None):
                if (priority_val not in ['numerical_high', 'numerical_low', 'numerical_prob']):
                    if (priority_val in cdf[priority_col].value_counts()):
                        cdf = cdf.loc[cdf[priority_col] == priority_val]
                    else:
                        priority_list = []
                        for s in cscaffs:
                            cmpds = scaffdict[s]
                            tempdf = pd.DataFrame(cmpds)
                            tempdf = tempdf.loc[tempdf[priority_col] == priority_val]
                            for idx, c in tempdf.iterrows():
                                if c[priority_col] == priority_val:
                                    priority_list.append(c)
                        if len(priority_list) != 0:
                            cdf = pd.DataFrame(priority_list)
                else:
                    if priority_val in ['numerical_high', 'numerical_low', 'numerical_prob']:
                        allmax = None
                        allmin = None
                        for s in cscaffs:
                            cmpds = scaffdict[s]
                            tempdf = pd.DataFrame(cmpds)
                            maxval = tempdf[priority_col].max()
                            minval = tempdf[priority_col].min()
                            if allmax is None or  maxval > allmax:
                                allmax = maxval
                            if allmin is None or maxval > allmax:
                                allmin = minval
                        if priority_val in ['numerical_high',  'numerical_prob']:
                            searchval = allmax
                        if priority_val == 'numerical_low':
                            searchval = allmin

                        if (searchval in cdf[priority_col].value_counts()):
                            cdf = cdf.loc[cdf[priority_col] == searchval]
                        else:
                            priority_list = []
                            for s in cscaffs:
                                 cmpds = scaffdict[s]
                                 tempdf = pd.DataFrame(cmpds)
                                 tempdf = tempdf.loc [tempdf[priority_col] == searchval]
                                 for idx, c in tempdf.iterrows ():
                                      if c [priority_col] == searchval:
                                            priority_list.append(c)
                            if len (priority_list) != 0:
                                cdf = pd.DataFrame (priority_list)


        #get structure w/ min # similars
            fps = []
            simcts = []

            if priority_val == 'numerical_prob':
                prob = searchval
                rnum = random.random()
                if rnum > prob:
                    continue

            for ix, c in cdf.iterrows ():
                mol = Chem.MolFromSmiles(c[smilescol])
                fp = Chem.GetHashedMorganFingerprint(mol, radius=2, nBits=1024, useFeatures=False)
                fps.append (fp)
                res = Chem.DataStructs.BulkTanimotoSimilarity(fp, out_fps)
                simct = len([element for element in res if element > sim_cutoff])
                simcts.append (simct)

            minidx = simcts.index (min (simcts))

            out_fps.append(fps [minidx])
            c = cdf.iloc[minidx]
            if c[smilescol] in out_list:
                print("Duplicate Error")
            out_list.append(c[smilescol])

            rowstring = ''
            for e in c:
                if ',' in str (e):
                    e = '\"' + e + '\"'
                rowstring = rowstring  + str (e) + ','
            out_file.write ( rowstring)
            out_file.write ( str(cluster_num) + ',' + key  +',' +  str(c[origsmilescol]) + '\n')
            cluster_num += 1

        return out_list
    def Kmeans (self, cmpdlist, num_clusters):
        fplist, uniquect = self.Get_FPArray(cmpdlist)
        if num_clusters > uniquect:
            num_clusters = uniquect

        kmeans = KMeans(n_clusters=num_clusters, init='k-means++', max_iter=1000, n_init=3, verbose=False)
        y_pred = kmeans.fit_predict(fplist)
        res_dict = {}
        for idx in range(0, len(cmpdlist)):
            res_dict[cmpdlist[idx]]= y_pred[idx]
        closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, fplist)

        return res_dict, closest, num_clusters
    def Analyze_SimStats (self, cmpd_list):
        fp_list = []
        for smi in cmpd_list:
            mol = Chem.MolFromSmiles(smi)
            fp = Chem.GetHashedMorganFingerprint(mol, radius=2, nBits=1024, useFeatures=False)
            fp_list.append (fp)
        for idx in range (0,len(fp_list)):
            res = Chem.DataStructs.BulkTanimotoSimilarity(fp_list [idx], fp_list)
            print (len([element for element in res if element > .5]))
    def CheckListScaffoldMatch (self,rowinfo, scaffdict, m):
        s = MScaff.MakeScaffoldGeneric(m)
        s = MScaff.GetScaffoldForMol(s)
        scaffsmi = Chem.MolToSmiles(s)
        canon = Chem.CanonSmiles(scaffsmi)
        if canon in scaffdict:
            scaffdict[canon].append (rowinfo)

            ret = -1
        else:
            scaffdict [canon] = []
            scaffdict[canon].append (rowinfo)
            ret = 1

        return scaffdict, ret
    # def Filters (self, m, filter_dict):
    #     if filter_dict is None:
    #         return False;
    #     if ('MW_high' in filter_dict) or ('MW_low' in filter_dict):
    #         mw = Descr.MolWt(m)
    #     if 'MW_high' in filter_dict:
    #         if mw > filter_dict['MW_high']:
    #             #print('MW high', mw)
    #             return True
    #     if 'MW_low' in filter_dict:
    #         if mw < filter_dict['MW_low']:
    #             #print('MW low', mw)
    #             return True
    #     if ('tPSA' in filter_dict):
    #         tpsa = Descr.TPSA(m)
    #         if tpsa > filter_dict['tPSA']:
    #             # print('tpsa', tpsa)
    #             return True
    #     if ('logP_high' in filter_dict) or ('logP_low' in filter_dict):
    #         logP = Descr.MolLogP(m)
    #     if ('logP_high' in  filter_dict):
    #         logphigh = filter_dict['logP_high']
    #         if logP > logphigh: #aLogP
    #             # print ('logphigh', logP)
    #             return True
    #     if ('logP_low' in filter_dict):
    #         logplow = filter_dict['logP_low']
    #         if logP < logplow: #aLogP
    #             # print('loglow', logP)
    #             return True
    #
    #     return False
    def Get_FPArray(self, cmpdlist):
        def get_numpy_fp(smi):
            mol = Chem.MolFromSmiles(smi)
            fp = Chem.GetHashedMorganFingerprint(mol, radius=2, nBits=1024, useFeatures=False)
            arr = np.zeros((4, 1024), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            return arr, fp

        fplist = []
        fpvecs = []
        dupct = 0
        for smi in cmpdlist:
            arr, fp = get_numpy_fp(smi)
            if fp in fpvecs:
                dupct +=1
            fpvecs.append(fp)
            fplist.append(arr)
        uniquect = len(cmpdlist)- dupct
        return fplist, uniquect
    def Get_FPArrayAsVstack(self, cmpdlist, hidetqdm=False, verbose = True):
        fplist = []
        deletedfps = []
        ct =0
        cmpdlist = [None if type (x) ==float and math.isnan(x) else x for x in cmpdlist]
        for smi in tqdm(cmpdlist, total=len(cmpdlist), disable=hidetqdm):

            if not smi is None and smi != 'FAIL' and smi != 'ENUMERATION_FAIL':
                if type(smi) != str:
                    mol =None
                    print('FAIL: smiles not in correct format: ', smi)
                    return None,None
                mol = Chem.MolFromSmiles(smi)
                try:
                    fp = Chem.GetMorganFingerprintAsBitVect(mol, radius = 2, nBits = 1024, useFeatures = False)
                except:
                    if verbose:
                        print ('Error in structure')
                    fp = None
                    deletedfps.append (ct)
            else:
                fp = None
                deletedfps.append(ct)
            if fp is not None:
                fplist.append(fp)
            ct += 1
        return np.vstack(fplist), deletedfps
    def Get_FPArray_dask(self, inlist, rettype = 'vstack', NUM_WORKERS=16, showprog = True, hidetqdm=False, verbose = True):
        def fptaskfcn(row):
            if row == 'FAIL' or row == 'ENUMERATION_FAIL':
                return None
            try:
                m = Chem.MolFromSmiles(row)
                fp = Chem.GetMorganFingerprintAsBitVect(m, radius = 2, nBits = 1024, useFeatures = False)
            except:
                fp = None
            return fp

        if (type(inlist[0]) == str):
            slist = inlist
        else:
            slist = []
            for c in inlist:
                slist.append(c[0])
        inseries = pd.Series(slist)

        ddf = dd.from_pandas(inseries, npartitions=CPU_COUNT * 10)
        if showprog:
            pbar = ProgressBar()
            pbar.register()
        res = ddf.apply(fptaskfcn, args=(), meta=(0, object)).compute(scheduler='processes',  num_workers=NUM_WORKERS)
        if showprog:
            pbar.unregister()

        print ('find remove list')
        deletedfps = []
        fplist = []
        for ix in range (0, len (res)):
            if res [ix] is None:
               deletedfps.append (ix)
            else:
                fplist.append(res[ix])
        if rettype == 'vstack':
            print('array')
            fplist = np.array (fplist)
        print ('fps complete')
        return fplist, deletedfps
    def DiSE_Selection (self, infile, outfile, rankcol, simcutoff, scaffold_sim_cutoff, n_chosen, minrankscore, maxviolations, hardsimcutoff, maxbbusage):
        df = pd.read_csv(infile)
        df = df.sort_values (rankcol, ascending=False)
        df ['chosen'] = 'Unknown'
        chosen = 0

        chosen_strux = []
        min_score = 1
        numexamined = 0
        bbusage = {}
        bbusageelim = 0
        for ix, r in df.iterrows ():
            if (r [rankcol] < minrankscore ):
                break
            numexamined += 1
            smi = r ['smiles']
            df.loc [ix,'chosen'] ='Yes'

            m = Chem.MolFromSmiles(smi)
            mfp = Chem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024, useFeatures=False)
            sc = MScaff.MakeScaffoldGeneric(m)
            sc = MScaff.GetScaffoldForMol(sc)
            scfp = Chem.GetMorganFingerprintAsBitVect(sc, radius=2, nBits=1024, useFeatures=False)
            violations = 0
            if r['BB1'] in  bbusage and bbusage [r['BB1']] >= maxbbusage:
                df.loc[ix, 'chosen'] = 'No'
                bbusageelim += 1
            else:
                for comp_mfp in chosen_strux:
                    tan = Chem.DataStructs.FingerprintSimilarity(mfp, comp_mfp [0] )
                    if tan >= simcutoff:
                        violations += 1
                        if (violations >= maxviolations or tan >= hardsimcutoff):
                            df.loc[ix, 'chosen'] = 'No'
                        continue
                    scaff_tan = Chem.DataStructs.FingerprintSimilarity(scfp, comp_mfp[1])
                    if scaff_tan >= scaffold_sim_cutoff:
                        violations += 1
                        if (violations >= maxviolations ):
                            df.loc[ix, 'chosen'] = 'No'
                        continue
            if df.loc[ix, 'chosen'] != 'No':
                df.loc[ix, 'chosen'] = 'Yes'
                chosen_strux.append([mfp, scfp])
                chosen += 1
                if r[rankcol] < min_score:
                    min_score = r[rankcol]
                if chosen > n_chosen:
                    break
                bb1 = r['BB1']
                bb2 = r['BB2']
                if bb1 not in bbusage:
                    bbusage [bb1] = 1
                else:
                    bbusage [bb1] += 1
                if bb2 not in bbusage:
                    bbusage [bb2] = 1
                else:
                    bbusage [bb2] += 1
                print (chosen, ' out of ', (numexamined), min_score, bbusageelim)

        print('final:' , min_score)
        exdf = df[df['chosen'] == 'Yes']
        exdf.to_csv (outfile, index = False)
            # chosen += 1
            # if chosen >= n_chosen:
            #     break
    def generatefps_dask(self, inlist, showprog=True, bitvec=False, SMILEScolname = 'SMILES'):
        def fptaskfcn(row, bitvec, showprog):
            try:
                if type (row) == str:
                    m = Chem.MolFromSmiles(row)
                else:
                    m = Chem.MolFromSmiles(row[SMILEScolname])
                if not bitvec:
                    fp = Chem.GetHashedMorganFingerprint(m, radius=2, nBits=1024, useFeatures=False)
                else:
                    fp =Chem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024, useFeatures=False)
            except:
                if showprog:
                    print('bad row:', row)
                return None
            return fp
        dfmode = False
        if (type(inlist) == pd.DataFrame):
            inseries = inlist
            if 'SMILES' not in inseries.columns:
                print ('FAIL: No SMILES column')
                return None
            dfmode = True
        else:
            if (type(inlist[0]) == str):
                slist = inlist
            else:
                slist = []
                for c in inlist:
                    slist.append(c[0])
            inseries = pd.Series(slist)

        ddf = dd.from_pandas(inseries, npartitions=CPU_COUNT * 10)
        if showprog:
            pbar = ProgressBar()
            pbar.register()

        if dfmode:
            res = ddf.apply(fptaskfcn, axis = 1, args=([bitvec, showprog]), meta=(0, object)).compute(scheduler='processes',  num_workers=NUM_WORKERS)
        else:
            res = ddf.apply(fptaskfcn, args=([bitvec, showprog]), meta=(0, object)).compute(scheduler='processes',  num_workers=NUM_WORKERS)

        if showprog:
            pbar.unregister()
        fplist = []
        for f in res:
            fplist.append(f)
        return fplist
    def compfps_dask(self, fplist, compfplist, showprog=True):
        def comptaskfcn(row, compfplist):
            try:
                tanlist = [round (x,2) for x in Chem.DataStructs.BulkTanimotoSimilarity(row, compfplist)]
            except:
                return [None] * len(compfplist)
            return tanlist
        dfmode = False
        if type(fplist) == pd.DataFrame:
            fp_series = fplist
            dfmode = True
        else:
            fp_series = pd.Series(fplist)
        ddf = dd.from_pandas(fp_series, npartitions=CPU_COUNT * 10)

        if showprog:
            pbar = ProgressBar()
            pbar.register()
        if dfmode:
            res = ddf.apply(comptaskfcn, axis = 1 , args=([compfplist]), meta=(0, object)).compute(scheduler='processes',  num_workers=NUM_WORKERS)
        else:
            res = ddf.apply(comptaskfcn, args=([compfplist]), meta=(0, object)).compute(scheduler='processes',  num_workers=NUM_WORKERS)
        if showprog:
            pbar.unregister()

        return list(res)




class Similarity:
    chksz = 10000000
    enum = Enumerate.Enumerate ()
    def PrepBBs(self, bbdict, libname, rxnschemefile, sim_column='dummy_smiles'):
        scheme, rxtnts = self.enum.ReadRxnScheme(rxnschemefile, libname, FullInfo=True)
        dummyrxtnts = scheme['scaffold_dummy_structures']
        cycs = []
        for cyc in bbdict.keys():
            cycs.append(cyc)
        for cx in range(0, len(cycs)):
            if bbdict [cycs[cx]] is not None:
                bbdict[cycs[cx]]['dummy_smiles'] = None
                bbdict[cycs[cx]]['skel_smiles'] = None
                bbdict[cycs[cx]]['scaff_smiles'] = None
                rxtnts = copy.deepcopy(dummyrxtnts)
                for idx, row in tqdm(bbdict[cycs[cx]].iterrows(), total=len(bbdict[cycs[cx]])):
                    rxtnts[cx] = row['SMILES']
                    res = self.enum.RunRxnScheme(rxtnts, rxnschemefile, libname, False)

                    if res[0] != 'FAIL':
                        sdx = Chem.CanonSmiles(res[0])
                        bbdict[cycs[cx]].loc[idx, 'dummy_smiles'] = sdx
                        m = Chem.MolFromSmiles(sdx)
                        try:
                            scaff = MScaff.GetScaffoldForMol(m)
                            skel = MScaff.MakeScaffoldGeneric(m)
                            skel = MScaff.GetScaffoldForMol(skel)
                            scx = Chem.MolToSmiles(scaff)
                            skx = Chem.MolToSmiles(skel)
                            bbdict[cycs[cx]].loc[idx, 'skel_smiles'] = skx
                            bbdict[cycs[cx]].loc[idx, 'scaff_smiles'] = scx
                        except:
                            bbdict[cycs[cx]].loc[idx, 'skel_smiles'] = 'FAIL'
                            bbdict[cycs[cx]].loc[idx, 'scaff_smiles'] = 'FAIL'
                    else:
                        bbdict[cycs[cx]].loc[idx, 'dummy_smiles'] = 'FAIL'
                        bbdict[cycs[cx]].loc[idx, 'skel_smiles'] = 'FAIL'
                        bbdict[cycs[cx]].loc[idx, 'scaff_smiles'] = 'FAIL'
                bbdict[cycs[cx]] = bbdict[cycs[cx]][bbdict[cycs[cx]][sim_column] != 'FAIL']
        return bbdict

    def Cluster_BBs(self, path, infiles, libname, bbspec, libspec, BBIDcolnames, SMILEScolnames, rxnschemefile, nclusters, mw_cutoff, simcolumn = 'dummy_smiles'):
        if infiles is not None:
            infilelist = infiles
        else:
            infilelist = self.enum.Get_BBFiles(bbspec, libspec, path, libname)
        cycdict = self.enum.load_BBlists(infilelist, BBIDcolnames=BBIDcolnames, SMILEScolnames=SMILEScolnames)

        for k in cycdict.keys ():
            if cycdict[k] is not None:
                print('START LEN', len(cycdict[k]))
                cycdict[k] = cycdict[k][cycdict[k]['MWT']<=mw_cutoff]
                print('END LEN', len(cycdict[k]))
                print (len(cycdict[k]))
        bbdict = self.PrepBBs (cycdict, libname, rxnschemefile, simcolumn )
        div = Diversity ()

        for k in bbdict.keys ():
            if bbdict[k] is not None:
                fps = div.generatefps_dask (bbdict [k], bitvec=True, SMILEScolname = simcolumn)
                for x in range (0, len (bbdict[k])):
                    if fps [x] is None:
                        print  (bbdict[k].iloc [x], 'None')


                k_means = KMeans(n_clusters=nclusters)
                k_means.fit(fps)
                bbdict [k]['clustnum'] = k_means.labels_
                bbdict [k].to_csv (path + '/' + libname + '/BBLists/' +libname + '.'  + k + '.clust.csv')
        return
    def Find_MostSimilarBB(self, bbsmiles, comp_bbdf, rxnschemefile, schemename, rxtntnum, retct = 1, excludelist = None, addtlcols=None):
        dummyscaff  = self.enum.Enumerate_Dummy_Scaffold (rxnschemefile, schemename, bbsmiles, rxtntnum)
        reslist = []
        for idx, row in comp_bbdf.iterrows ():
            compsmiles = row ['SMILES']
            if compsmiles == 'FAIL':
                reslist.append([row['BB_ID'], 0, 'FAIL'])
            else:
                comp_dummyscaff = self.enum.Enumerate_Dummy_Scaffold (rxnschemefile, schemename, compsmiles, rxtntnum)
                tnm = self.TanimotoComparison(dummyscaff, comp_dummyscaff)
                if addtlcols is not None:
                    newrow = [row['BB_ID'], tnm, row['SMILES']]
                    for c in addtlcols:
                        newrow.append (row [c])
                    reslist.append(newrow)
                else:
                    reslist.append ([row ['BB_ID'],  tnm, row ['SMILES']])
        reslist = sorted(reslist, key=operator.itemgetter(1), reverse=True)

        if excludelist is None:
            return reslist[0:retct]
        ct = 0
        next = 0
        outlist = []
        while ct < retct:
            if reslist[next] [0] not in excludelist:
                outlist.append(reslist[next] )
                ct += 1

            next += 1
        return outlist

    def TanimotoComparison (self, smiles1, smiles2 ):
        #currently hardcoded to Morgan2
        try:
            if smiles1 == 'FAIL' or smiles2 == 'FAIL':
                return -1
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            fp1 = Chem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
            fp2 = Chem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            tnm_score = round(DataStructs.TanimotoSimilarity(fp1, fp2), 2)
        except:
            tnm_score = -1

        return tnm_score
    def SimSearchLargeFile (self, smileslist, catfile, outpath, sep):
        def processchunk (smileslist, df, outpath, smilesfps):
            searchdf, smilesfps = self.Run_SimilarityComp(smileslist, df, outpath, wmode, smilesfps)
            return searchdf, smilesfps


        reader = pd.read_csv(catfile, chunksize=self.chksz, sep=sep)
        cct = 0
        wmode = 'w'
        smilesfps = None
        for chunk in reader:
            resdf = copy.deepcopy (pd.DataFrame(chunk))
            resdf = resdf.reset_index()
            if 'smiles' in resdf.columns:
                resdf = resdf.rename (columns = {'smiles':'SMILES'})

            print('CHUNK', cct)

            df, smilesfps = processchunk(smileslist, resdf, outpath, smilesfps)

            wmode = 'a'
            cct += 1
        return
    def Run_SimilarityComp(self, smileslist, searchstrux, outpath, wmode, smilesfps):
        div = Diversity()

        if type(smileslist) == str:
            smileslist = [smileslist]

        if type(searchstrux) == str:
            searchdf = pd.read_csv(searchstrux)
        else:
            searchdf = searchstrux
        if 'SMILES' not in searchdf.columns:
            return None, None
        if smilesfps is None:
            smilesfps = div.generatefps_dask(smileslist, showprog=True)

        start = time.perf_counter()
        searchfps = div.generatefps_dask(searchdf, showprog=True)
        next = time.perf_counter()
        print (next-start)
        start = next
        res = div.compfps_dask( smilesfps, searchfps )
        coldict = {}
        for s in range (0,len(smileslist)):
            coldict [s] ='query' + str(s)
        resdf = pd.DataFrame(res).transpose ()
        resdf= resdf.rename (columns = coldict)
        searchdf = pd.concat([searchdf, resdf], axis=1)
        next = time.perf_counter()
        print(next - start)
        start = next
        if wmode == 'a':
            writehdrs = False
        else:
            writehdrs = True
        search_ddf = dd.from_pandas(searchdf, npartitions=CPU_COUNT)
        search_ddf.to_csv(outpath, mode=wmode, index=False, header=writehdrs)
        next = time.perf_counter()
        print(next - start)
        start = next
        return searchdf, smilesfps
    def CompareMolecule(self, line, lct, cutoff, probefps, probesmiles, outfile):
        splitstr = line.split('\t')
        id = splitstr[1]
        smiles = splitstr[0]
        mol1 = Chem.MolFromSmiles(smiles)
        fp2 = Chem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)

        probenum = 0
        probelist = []
        tvscores = []
        tnmscores = []
        for fp1 in probefps:
            tv_score = round(DataStructs.TverskySimilarity(fp1, fp2, 0, 1), 2)
            if (tv_score > cutoff):
                probelist.append(probenum)
                tvscores.append(tv_score)
                tnm_score = round(DataStructs.TanimotoSimilarity(fp1, fp2), 2)
                tnmscores.append(tnm_score)
            probenum += 1

        if len(probelist) > 0:
            if len(probelist) == 1:
                maxidx = 0
            else:
                max_value = max(tvscores)
                maxidx = tvscores.index(max_value)
            splitsmiles = smiles.split(' ')
            smiles = splitsmiles[0]
            if len(splitsmiles) > 1:
                extrainfo = splitsmiles[1].replace(',', '&#44')
            else:
                extrainfo = ''
            outline = ",".join([smiles, extrainfo, str(lct), str(probelist).replace(',', ';'),
                                probesmiles[probelist[maxidx]], str(probelist[maxidx]), str(tvscores[maxidx]),
                                str(tnmscores[maxidx]),
                                str(len(probelist))]) + '\n'
            outfile.write(outline)
            return True
        else:
            return False
    def Run_SimList(self, probefile, catfile, outfile):
        N = 10000
        ct = 0

        hdrread = False
        probefps = []
        probesmiles = []
        with open() as probes:
            for mlist in iter(lambda: list(islice(probes, N)), []):
                for r in mlist:
                    if hdrread == True:
                        splitstr = r.split(',')
                        id = splitstr[0]
                        smiles = splitstr[32]
                        mol1 = Chem.MolFromSmiles(smiles)
                        info = {}
                        fp1 = Chem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=1024)
                        probefps.append(fp1)
                        probesmiles.append(smiles)
                    else:
                        hdrread = True
        probes.close()

        hdrread = False
        pool_size = 16
        lct = 0
        simct = 0
        cutoff = .77
        outfile = open(outfile, 'w')

        with open(catfile) as collection:
            for mlist in iter(lambda: list(islice(collection, N)), []):
                pool = Pool(pool_size)
                block = ''
                for line in mlist:
                    if hdrread == True:
                        hassim = pool.apply_async(self.CompareMolecule,args=(line, lct, cutoff, probefps, probesmiles, outfile)).get()
                        if hassim == True:
                            simct += 1
                    else:
                        hdrread = True
                        block += 'Enamine SMILES,Addtl Info,line,probe list,max probe smiles,probe list #,tversky score, tanimoto score\tNum sims\n'
                    lct = lct + 1
                outfile.write(block)
                pool.close()
                pool.join()
                ct = ct + 1

        collection.close()
        outfile.close()
    def SubSearch(self, hdrlist, smilescol, line, lct, ss1, splitchar):
        fields = line.split(splitchar)
        smi = fields[smilescol]
        try:
            mol1 = Chem.MolFromSmiles(smi)
        except:
            return False
        try:
            in_struct = mol1.GetSubstructMatches(ss1)
        except:
            return False
        if len(in_struct) > 0:
            return True
        else:
            return False
    def Search_Substruct(self, catfile, substruct, outfilename, splitchar, Keep=True, Smilescolname="smiles"):
        N = 10000
        pool_size = 8
        subct = 0
        ss1 = Chem.MolFromSmarts(substruct)
        hdrread = False
        lct = 0

        outfile = open(outfilename, 'w')
        ct = 0
        with open(catfile) as collection:
            for mlist in iter(lambda: list(islice(collection, N)), []):
                pool = Pool(pool_size)
                block = ''
                for line in mlist:
                    if line == '\n':
                        continue
                    if hdrread == True:
                        hassub = pool.apply_async(self.SubSearch, args=(hdrlist, smilescol, line, lct, ss1, splitchar)).get()
                        if (Keep == True and hassub == True) or (Keep == False and hassub == False):
                            subct += 1
                            block += line.strip().replace(' ', ',') + '\n'
                    else:
                        hdrlist = line.split(splitchar)
                        matching = [s for s in hdrlist if Smilescolname in s.lower()]
                        smilescol = hdrlist.index(matching[0])
                        hdrread = True
                        block += line.strip().replace(' ', ',') + '\n'
                    lct = lct + 1
                outfile.write(block)
                pool.close()
                pool.join()
                ct = ct + 1

        collection.close()
        outfile.close()
    def FindClosestList_SimMatch (self, testlist, complists):
        res = []
        for c in complists:
            sumtanscores = 0
            tanct = 0
            for t in testlist:
                settanlist = []
                for cs in c:
                    tanscore = self.TanimotoComparison(t,cs)
                    if tanscore == -1:
                        print ('error:', cs, t)
                    else:
                        settanlist.append (tanscore)
                maxtan = max(settanlist)
                tanct += 1
                sumtanscores += maxtan
            res.append (sumtanscores/tanct)
        return res

class Chem_Similarity_DiversityUI:
    initpath = '../CIxTools.init.json'
    paramslist = ['simdiv_inputfile', 'simdiv_inputsmiles', 'bbsim_inputfile', 'bbsim_inputsmiles', 'bbsim_rxnschemefile','bbsim_schemename',
                  'bbsim_rxtntnum', 'bbsim_retct', 'bbsim_filtcol', 'bbsim_filtvals', 'bbsim_addtlcols', 'bbsim_savetofldr']
    sim = Similarity ()

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

