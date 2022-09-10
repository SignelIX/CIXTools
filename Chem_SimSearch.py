from rdkit import Chem
from itertools import islice
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from multiprocessing.pool import ThreadPool as Pool

def CompareMolecule(line, lct, cutoff, probefps, probesmiles, outfile):
    splitstr = line.split('\t')
    id = splitstr[1]
    smiles = splitstr[0]
    mol1 = Chem.MolFromSmiles(smiles)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)

    probenum = 0
    probelist = []
    tvscores = []
    tnmscores = []
    for fp1 in probefps:
        tv_score = round (DataStructs.TverskySimilarity(fp1, fp2, 0, 1),2)
        if (tv_score > cutoff):
            probelist.append (probenum)
            tvscores.append (tv_score)
            tnm_score = round(DataStructs.TanimotoSimilarity(fp1, fp2), 2)
            tnmscores.append (tnm_score)
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
            extrainfo = splitsmiles[1].replace (',', '&#44')
        else:
            extrainfo = ''
        outline = ",".join([ smiles , extrainfo , str(lct), str(probelist).replace(',', ';'),
                             probesmiles[probelist[maxidx]] , str(probelist[maxidx]), str(tvscores[maxidx]),str(tnmscores[maxidx]),
                            str(len(probelist))]) + '\n'
        outfile.write(outline)
        return True
    else:
        return False

def Run_SimList (probefile, catfile, outfile):
    N= 10000
    ct = 0
    writect = 1

    hdrread = False
    probefps = []
    probesmiles = []
    with open () as probes:
        for mlist in iter(lambda: list(islice(probes, N)), []):
            print (len(mlist))
            for r in mlist:
                if hdrread == True:
                    splitstr = r.split(',')
                    id = splitstr[0]
                    smiles = splitstr[32]
                    mol1 = Chem.MolFromSmiles(smiles)
                    info = {}
                    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=1024)
                    probefps.append(fp1)
                    probesmiles.append (smiles)
                else:
                    hdrread = True
    probes.close ()

    hdrread = False
    pool_size =  8
    lct = 0
    simct = 0
    cutoff = .77
    outfile = open (outfile, 'w')


    with open (catfile) as collection:
        for mlist in iter(lambda: list(islice(collection, N)), []):
            pool = Pool(pool_size)
            block = ''
            for line in mlist:
                if hdrread == True:
                    hassim = pool.apply_async(CompareMolecule, args= (line, lct, cutoff, probefps, probesmiles, outfile)).get()
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
            print (ct * N)
            print (simct)
        print ('joined')
    collection.close ()
    outfile.close ()

def SubSearch ( hdrlist, smilescol, line, lct, ss1, splitchar):
    fields = line.split (splitchar)
    smi = fields [smilescol]
    try:
        mol1 = Chem.MolFromSmiles(smi)
    except:
        return False
    try:
        in_struct = mol1.GetSubstructMatches (ss1)
    except:
        return False
    if len(in_struct) > 0:
        return True
    else:
        return False

def Search_Substruct (catfile, substruct, outfilename, splitchar,  Keep = True, Smilescolname = "smiles"):
    N = 10000
    pool_size =  8
    subct = 0
    ss1 = Chem.MolFromSmarts(substruct)
    hdrread = False
    lct = 0

    outfile = open (outfilename, 'w')
    ct = 0
    with open (catfile) as collection:
        for mlist in iter(lambda: list(islice(collection, N)), []):
            pool = Pool(pool_size)
            block = ''
            for line in mlist:
                if line == '\n':
                    continue
                if hdrread == True:
                    hassub = pool.apply_async(SubSearch, args= (hdrlist, smilescol,  line, lct, ss1, splitchar)).get()
                    if (Keep == True and hassub == True) or (Keep == False and hassub == False):
                        subct += 1
                        block += line.strip ().replace(' ',',') + '\n'
                else:
                    hdrlist = line.split (splitchar)
                    matching = [s for s in hdrlist if Smilescolname in s.lower ()]
                    smilescol = hdrlist.index (matching[0])
                    hdrread = True
                    block += line.strip().replace(' ',',') + '\n'
                lct = lct + 1
            outfile.write(block)
            pool.close()
            pool.join()
            ct = ct + 1
            print (ct * N)
            print (subct)
        print ('joined')
    collection.close ()
    outfile.close ()
