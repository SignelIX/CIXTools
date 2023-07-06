import multiprocessing
from rdkit.Chem import AllChem as Chem
import pandas as pd
import tqdm
# from threading import Thread
import multiprocessing
import time
from operator import itemgetter
from multiprocessing import Semaphore

def fpchunkfcn(args):
    chunk, bitvec , ct, queue, SMILEScolname = args [0], args[1], args [2], args [3], args [4]
    fplist = []
    for idx, row in chunk.iterrows() :
        try:
            if type(row) == str:
                m = Chem.MolFromSmiles(row)
            else:
                m = Chem.MolFromSmiles(row[SMILEScolname])
            if not bitvec:
                fp = Chem.GetHashedMorganFingerprint(m, radius=2, nBits=1024, useFeatures=False)
            else:
                fp = Chem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024, useFeatures=False)
            fplist.append([row [SMILEScolname], fp])
        except Exception as e:
            print (str(e))
            fplist.append([ row [SMILEScolname],None])
    queue.put ( fplist)

def generatefps_multiproc( filename, SMILEScolname='SMILES'):
    chunksize = 10 ** 5
    ct = 0
    tdf = pd.read_csv(filename)
    print (len (tdf))
    proclist = []
    x = time.perf_counter()
    queue = multiprocessing.Queue()
    for chunk in pd.read_csv(filename, chunksize=chunksize):
        #fpchunkfcn([chunk, False, ct])
        p =  multiprocessing.Process (target=fpchunkfcn, args = [[chunk, False, ct, queue, SMILEScolname]])
        proclist.append (p)
        p.start ()
        ct += 1


    rets = []
    for p in proclist:
        ret = queue.get()
        rets.extend (ret)
        print (len(rets))
    for p in proclist:
        p.join ()

    y = time.perf_counter()
    print(y - x)
    print ('complete')
    rets = [i for i in rets if i[1] is not None]
    return rets

def comptaskfcn(args):
    fplist, compfplist, queue, sema= args[0], args[1], args[2], args [3]
    try:
        tanlist = [round(x, 2) for x in Chem.DataStructs.BulkTanimotoSimilarity(fplist, compfplist)]
    except Exception as e:
        tanlist = [None] * len(compfplist)
        print (str(e))

    queue.put ( tanlist)
    sema.release ()

def compfps_multiproc( fplist, compfplist):
    proclist = []
    n = 10 ** 5
    x = time.perf_counter()
    queue = multiprocessing.Queue()
    ct = 0
    numprocs = 10
    sema = Semaphore (numprocs)
    for i in range(0, len(compfplist), n):
        print ('starting', end= '\r')
        p = multiprocessing.Process(target=comptaskfcn, args=[[ fplist, compfplist [i:i + n], queue, sema]])
        proclist.append(p)
        p.start()
        print(ct, end='\r')
        ct += 1
    print('\n')
    rets = []
    ct = 0
    for p in proclist:
        ret = queue.get()
        rets.extend(ret)
        print (ct, end = '\r')
        ct += 1
    ct = 0
    for p in proclist:
        p.join()
        print(ct, end='\r')
        ct += 1

    y = time.perf_counter()
    print(y - x)
    print('complete compfps')
    return list(rets)


def Run_MultiSearch (infile,  smilescol, compfile, compsmilescol):
    df = pd.read_csv(
       infile)
    rets = generatefps_multiproc(compfile, SMILEScolname=compsmilescol)

    for idx, row in df.iterrows():
        print (row[smilescol])
        mol = Chem.MolFromSmiles(row[smilescol])
        fp = Chem.GetHashedMorganFingerprint(mol, radius=2, nBits=1024, useFeatures=False)
        compfplist = list(map(itemgetter(1), rets))
        print (len(compfplist))
        tanimotos = compfps_multiproc(fp, compfplist)
        print(len(tanimotos))
        tanimotos = sorted(tanimotos, key=lambda x: (not (x is None), x), reverse=True)
        print(tanimotos[0:100])
