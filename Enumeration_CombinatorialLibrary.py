
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import Enumerate
enum = Enumerate.Enumerate()
NUM_WORKERS = 16
import random

def Enumerate_FullSpecifx_File (infile, rxschemefile, schemename):
    df = pd.read_csv (infile)
    df ['full_smiles' ] = None
    for ix, row in df.iterrows ():
        in_reactants = [row ['bb1_smiles'], row['bb2_smiles'], row ['bb3_smiles']]
        res = enum.RunRxnScheme(in_reactants, rxschemefile, schemename, False, None)
        df.loc [ix, 'full_smiles'] = res [0]
    df.to_csv (infile)

def EnumTuples_dask (rxnschemefile, libname, tuples, bbdict, cycs, tupleinfodict, showprog = True):
    def taskfcn ( row, rxnschemefile, bbdict, cycs, tupleinfodict, schemeinfo):
        if str(row) in tupleinfodict :
            return tupleinfodict[str(row)], row
        rxtnts = []
        for cx in range(0, len(cycs)):
            rxtnts.append(bbdict[cycs[cx]].iloc[row[cx]]['SMILES'])

        res, prodct, schemeinfo = enum.RunRxnScheme(rxtnts, rxnschemefile, libname, False, schemeinfo)
        tupleinfodict[str(row)] = res
        return res, row


    tuple_series = pd.Series(tuples)
    ddf = dd.from_pandas(tuple_series, npartitions=NUM_WORKERS)
    if showprog:
        pbar = ProgressBar()
        pbar.register()
    schemeinfo = enum.ReadRxnScheme(rxnschemefile, libname, False)
    res = ddf.apply(taskfcn, args=(rxnschemefile, bbdict, cycs, tupleinfodict, schemeinfo), meta=(0, object)).compute()
    if showprog:
        pbar.unregister()
    return res

def recurs_tuples (indict, intuples, cyclist, cyc ):
    tuples = []
    for b in indict [cyc]:
        if intuples is not None:
            for t in intuples:
                temp = t.copy ()
                temp.append (b)
                tuples.append (temp)
        else:
            tuples.append ([b])
    if cyclist.index(cyc) == len(cyclist)-1:
        return tuples
    else:
        return recurs_tuples(indict, tuples, cyclist, cyclist [cyclist.index(cyc) + 1] )

def Generate_BB_Tuples (rxnschemefile,libname,  statcyc, initset, bbixes, cycs, test_BBN, tupleinfodict, bbdict):
    testsamp = {}
    # choose random sample of structures to test  with curr static BB
    for c in cycs:
        if c != statcyc:
            testsamp[c] = []
            idxs = random.choices(range(0, len(initset[c])), k=test_BBN)
            for ix in idxs:
                testsamp[c].append(initset[c][ix])

    currlist_dict = {}

    # # generate tuples from the list
    for bbix in bbixes:
        # add the BB id numbers
        testsamp [statcyc] = [bbix]
        currset = {}
        for c in cycs:
            currset[c] = []
            for i in range(0, len(testsamp[c])):
                currset[c].append(testsamp[c][i])
        testsamp[statcyc] = [bbix]
        currtuples = []
        for i in range(0, test_BBN):
            slist = []
            for c in cycs:
                if c == statcyc:
                    slist.append(testsamp[c][0])
                else:
                    slist.append(testsamp[c][i])
            currtuples.append(slist)
        currlist = EnumTuples_dask(rxnschemefile, libname, currtuples, bbdict, cycs, tupleinfodict, False)
        currlist_dict [bbix] = currlist
    return currlist_dict