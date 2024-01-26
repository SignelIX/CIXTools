from rdkit.Chem import SaltRemover
from rdkit import Chem
from rdkit.Chem import Descriptors as Descr
import matplotlib.pyplot as plt
from itertools import islice
from multiprocessing.pool import ThreadPool as Pool
from tqdm import tqdm
import pandas as pd
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from Archive import Enumerate
import Chem_CalcProps
from rdkit import RDLogger
import time
import threading

def SaltStripFile (infile, smilescol, outfile):
    df  = pd.read_csv (infile)
    df = SaltStripMolecules(df, smilescol)
    df.to_csv(outfile, index=False)

def SaltStripMolecules (molecules: pd.DataFrame, smilescol='SMILES', neutralize = True):
    tqdm.pandas ()
    print ('salt stripping')
    molecules[smilescol] = molecules[smilescol].progress_apply(lambda smi:SaltStrip(smi, neutralize))
    print ('completed salt stripping')
    return molecules

def SaltStrip (molec, neutralize = True):
    try:
        saltstrip = SaltRemover.SaltRemover()
        if type(molec) != Chem.Mol:
            m = Chem.MolFromSmiles(molec)
        else:
            m = molec
        m=  saltstrip.StripMol(m)
        smi = Chem.MolToSmiles(m, kekuleSmiles = False)
        smi = Neutralize(smi)
    except:
        return molec
    return smi

def Neutralize (smi):
    try:
        m = Chem.MolFromSmiles(smi)
    except:
        return smi
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    try:
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = m.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()

        smi = Chem.MolToSmiles(m, kekuleSmiles=False)
        return smi
    except:
        print (smi)
        return smi

def Kekulize_Smiles (smi):
    m = Chem.MolFromSmiles(smi)
    Chem.Kekulize(m)
    smi = Chem.MolToSmiles(m, kekuleSmiles=False)
    return smi

def Molecule_FromSmiles (smi):
    m = Chem.MolFromSmiles(smi)
    return m

def Molecule_ToSmiles (m):
    smi = Chem.MolToSmiles(m)
    return smi

def Smiles_FromMolecule (m):
    smi = Chem.MolToSmiles(m)
    return smi

def Get_FragCount (m):
    frags = Chem.GetMolFrags(m)
    return len(frags)

def Deprotect (compound, rxnscheme, deprotect_name=None, iterate = False, retsmiles = False):
    deprotected = False
    Enumerator = Enumerate.Enumerate ()
    first_round = True

    has_scheme = Enumerator.ReadRxnScheme(rxnscheme, 'Deprotection',  verbose=True)
    schemeinfo = [[{"Reactants":["r0","None"],"Rxns":{"default":[]}}], ["r0","None"] ]
    if deprotect_name is None:
        for k,v in has_scheme[0][0]['Rxns'].items ():
            if type(v) is list:
                for l in v:
                    schemeinfo[0][0]['Rxns']["default"].append(l)
            else:
                schemeinfo[0][0]['Rxns']["default"].append(v)
    else:
        if type (deprotect_name) is str:
            deprotect_name = [deprotect_name]
        for d in deprotect_name:
            v = has_scheme[0][0]['Rxns'] [d]
            if type(v) is list:
                for l in v:
                    schemeinfo[0][0]['Rxns']["default"].append(l)

    if not has_scheme[0] is  None:
        last_dp = False
        while first_round or  (iterate and last_dp == True):
            res, prod_ct, resinfo = Enumerator.RunRxnScheme([compound, 'None'], None,None, False, schemeinfo )
            if (res != 'FAIL' and res != 'None'):
                if type(res) == str:
                    compound = Chem.MolFromSmiles(res)
                deprotected = True
                last_dp = True
            else:
                last_dp = False

            first_round = False

        first_round = True
        last_dp = False
        res_list = [compound]

    has_scheme = Enumerator.ReadRxnScheme(rxnscheme, 'Deprotect_KeepBoth', False)
    if has_scheme[0] is not None:
        while first_round or (iterate and last_dp == True):
            res, prod_ct = Enumerator.RunRxnScheme([compound, 'None'], rxnscheme, 'Deprotect_KeepBoth', False)
            if res == 'NOSCHEMEFAIL':
                break
            if (res != 'NOSCHEMEFAIL' and res != 'FAIL' and res != 'None'):
                res_list.append (Chem.MolFromSmiles(res))
                compound = Chem.MolFromSmiles(res)
                deprotected = True
                last_dp = True
            else:
                last_dp = False
            first_round = False


    res = res_list
    if retsmiles:
        if deprotected:
            return Chem.MolToSmiles(res[0])
        else:
            return compound
    return res, deprotected

def DeprotectFile ( df, n_analyzed, idcol, smilescol,origsmilescol,deprotect_specfile, outfile):
    tct = 0
    deprot_df =  df.drop(df.index)

    if n_analyzed != -1:
        df = df.head(n_analyzed)

    for index, row in tqdm(df.iterrows(),total=len(df)) :
        if origsmilescol == -1:
            origsmilescol = len (row)
        if not (n_analyzed == -1 or tct < n_analyzed):
            break

        m = Chem.MolFromSmiles(row[smilescol])
        if deprotect_specfile is not None:
            try:
                res, deprotected = Deprotect (m, deprotect_specfile, True)
            except:
                deprotected = False
                res = [m]
        else:
            deprotected = False

        rowcopy = None
        rowlist = []

        if deprotected:
            if len (res) ==1:
                row [origsmilescol] = row[smilescol]
                row [smilescol] = Chem.MolToSmiles(res [0])
                rowlist.append ([row, res[0]])
            else:
                orig = row[smilescol]
                row[smilescol] = Chem.MolToSmiles(res[0])
                rowlist.append([row, res[0]])
                for ixr in range (1,len(res)):
                    rowcopy = row.copy()
                    rowcopy[smilescol]  = Chem.MolToSmiles(res [ixr])
                    rowcopy [origsmilescol] = orig
                    rowlist.append ([rowcopy, res[ixr]])
        else:
            rowlist.append ([row, m])


        if rowcopy is None:
            rows = [row]
        else:
            rows = [row, rowcopy]

        # if (len (rows)>2 ):
        #     print (len(rows))
        deprot_df = deprot_df.append(rows, ignore_index=True)

    deprot_df.to_csv (outfile, index = False)

def Add_Properties (filename, prop_list):
    df = pd.read_csv(filename)
    if 'MW' in prop_list:
        df = df.assign(MW=None )
    if ('tPSA' in prop_list):
        df = df.assign (tPSA = None)
    if ('logP' in prop_list):
        df = df.assign (logP = None)

    for idx, row in df.iterrows () :
        smi = row['SMILES']
        try:
            m = Chem.MolFromSmiles(smi)
            if 'MW' in prop_list:
                mw = Descr.MolWt (m)
                df.at[idx,'MW'] = mw
            if ('tPSA' in prop_list):
                tpsa = Descr.TPSA(m)
                df.at[idx,'tPSA'] = tpsa
            if ('logP' in  prop_list):
                logP = Descr.MolLogP (m)
                df.at[idx,'logP'] = logP
            if idx%1000 == 0:
                print (idx, end = '\r')
        except:
            continue

    df.to_csv (filename.replace ('.csv', '.props.csv'))

def PlotProperties (filename, cols, bins):
    df = pd.read_csv(filename)
    for c in cols:
        df.hist (column = c, bins = bins [c], rwidth=0.9)
        plt.show ()

def Filters ( m, filter_dict, stats):
    cp = Chem_CalcProps.Chem_CalcProps ()
    return cp.Filters(m, filter_dict, stats)

def ApplyFilters ( smi, filter_dict, ssfilters, useChirality, AmbiguousChirality, Keep, retbool=True):
    try:
        m = Chem.MolFromSmiles(Kekulize_Smiles( smi))
    except:
        return Keep, [None, None, None]
    propFiltered,stats = Filters(m , filter_dict, stats =None)
    if ssfilters is not None:
        ssFiltered, df = Substructure_Filters (m, None, ssfilters, useChirality, Keep)
    else:
        ssFiltered = False
    ambig = False

    if (AmbiguousChirality == True):
        ambig = ContainsUnresolvedChirality(m)

    if (Keep):
        if retbool:
            return propFiltered and ssFiltered and ambig
        return propFiltered and ssFiltered and ambig, [propFiltered, ssFiltered, ambig]
    else:
        if retbool:
            return propFiltered or ssFiltered or ambig
        return propFiltered or ssFiltered or ambig, [propFiltered, ssFiltered, ambig]

def Process_Block (mlist, smilescol, splitchar, filter_dict, ssfilters, useChirality, AmbigChirality, Keep, blocknum):
    block = ''
    removedblock = ''
    ct = 0
    for line in mlist:
        if line == '\n':
            continue
        smi = line.split(splitchar)[smilescol]
        isFiltered = ApplyFilters(smi, filter_dict, ssfilters, useChirality, AmbigChirality, Keep)
        if (Keep == True and isFiltered == True) or (Keep == False and isFiltered == False):
            block += line.strip().replace(' ', ',') + '\n'
        else:
            removedblock += line.strip().replace(' ', ',') + '\n'
        ct += 1
    return block, removedblock

class CompleteBlock ():
    tic = None
    lock = None
    outfile = None
    blockct = 0

    def __init__ (self, outf, remf):
        self.tic = time.perf_counter()
        self.lock = threading.Lock()
        self.outfile = outf
        self.removedfile = remf

    def CompleteBlockAsync (self, blocklist):
        if blocklist is not None:
            self.lock.acquire()
            self.blockct = self.blockct + 1
            self.outfile.write(blocklist [0] + '\n')
            self.removedfile.write(blocklist[1] + '\n')
            self.lock.release()
            print ('Completing BLOCK:', self.blockct, '          ', end = '\r' )


def Filter_File (catfile, outfilename, splitchar, filter_dict, ss_file, useChirality, AmbigChirality, Keep = False, smilescol = 'SMILES'):
    print ('Filtering,', filter_dict)
    N = 10000
    pool_size =  40
    pool = Pool(pool_size)
    outfile = open (outfilename, 'w')
    remf =open (outfilename.replace ('.csv','') + 'removed.csv', 'w')
    if 'excl_ss_smiles' in filter_dict:
        ssfilters = Read_FilterFile(filter_dict ['excl_ss_smiles'])
    else:
        ssfilters = Read_FilterFile(ss_file)
    CB = CompleteBlock(outfile, remf)

    with open (catfile) as collection:
        line = collection.readline()
        hdrlist = line.split(splitchar)
        hdrlist = [sx.upper() for sx in hdrlist]
        matching = [s for s in hdrlist if "SMILES" in s]
        smilescol = hdrlist.index(matching[0])
        blocknum = 0
        outfile.write (line + '\n')

        for mlist in iter(lambda: list(islice(collection, N)), []):
            print('Starting BLOCK: ', blocknum, '           ', end='\r')
            blocknum +=1
            pool.apply_async(Process_Block, args=(mlist, smilescol, splitchar,  filter_dict, ssfilters, useChirality, AmbigChirality, Keep, blocknum), callback=CB.CompleteBlockAsync)

        pool.close()
        pool.join()

    collection.close ()
    outfile.close ()
    print('Filtering Complete, ' , outfilename)

def Read_FilterFile (filter_file_or_list):
    if filter_file_or_list is None:
        return None
    if type (filter_file_or_list) is list:
        ss_filters = pd.DataFrame (filter_file_or_list, columns = ['SMILES'])
    else:
        ss_filters = pd.read_csv(filter_file_or_list)
    ss_filters = ss_filters.assign(SS_ct=0, Molecule=None)
    for idx, row in ss_filters.iterrows():
        mol = Chem.MolFromSmarts(row['SMILES'])
        ss_filters.loc[idx, 'Molecule'] = mol
    return ss_filters

def SubstructureCheckMolecules (molecules: pd.DataFrame, ss_smarts):
    passixs = []
    smartsmol = Chem.MolFromSmarts(ss_smarts)
    for ix, row in molecules.iterrows ():
        m=Chem.MolFromSmiles(row['SMILES'])
        if m is not None and  m.HasSubstructMatch(smartsmol, useChirality=False) == True:
            passixs.append (ix)
    return molecules.iloc [passixs].reset_index ()


def SubstructureCheckSMILES(smiles, smarts, useChirality=False):
    m=Chem.MolFromSmiles(smiles)
    smartsmol = Chem.MolFromSmarts(smarts)
    return m.HasSubstructMatch(smartsmol, useChirality=useChirality)

def Substructure_Filters (input_mol, outfile, ss_filters, useChirality, Keep_IfNoFilter):
    if type (input_mol) == str:
        input_mol = Chem.MolFromSmiles(input_mol)
    if ss_filters is None:
        return Keep_IfNoFilter, None
    if 'Molecule' not in ss_filters.columns:
        ss_filters ['Molecule'] = ss_filters.apply (lambda x: Chem.MolFromSmarts(x['SMILES']), axis=1, result_type='expand')
    if 'SS_ct' not in ss_filters.columns:
        ss_filters ['SS_ct'] = 0
    for idx, row in  ss_filters.iterrows():
        if input_mol.HasSubstructMatch(row['Molecule'], useChirality = useChirality)== True:
            ss_filters.loc[idx, 'SS_ct'] = row.SS_ct + 1
            if outfile is not None:
                outfile.write (Chem.MolToSmiles(input_mol) + ',' + row.SMILES + '\n')
            return True, ss_filters
    return False, ss_filters

def ContainsUnresolvedChirality (m):
    isomers = tuple(EnumerateStereoisomers(m))
    if len(isomers) > 1:
        return True
    else:
        return False

def RemoveChirality_FromSMILES (smi):
    mol = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

def RemoveChirality (infile_or_df, outfile, smilescol):
    df = pd.read_csv(infile_or_df)
    df['achiral_smiles'] = ''
    df = df.rename(columns = {df.columns [0]:'orig_rownum'})

    for idx, row in df.iterrows():
        mol = Chem.MolFromSmarts(row[smilescol])
        try:
            Chem.RemoveStereochemistry(mol)
            df.loc[idx, 'achiral_smiles'] = Chem.MolToSmiles(mol)
        except:
            df.loc[idx, 'achiral_smiles'] = row [smilescol]
    df = df.sort_values (by=['orig_rownum'])
    cols = df.columns.tolist ()
    rarrcols = cols [0:1] + cols[-1:] + cols[1:-1]
    df = df[rarrcols]
    df.to_csv(outfile, index=None)

def SDFtoFile (infile, fix, outfile, idcol='idnumber' ):
    print ('starting conversion')
    # if fix == True:
    #     f=open (infile, 'rb')
    #     text = f.read().decode(errors='replace')
    #     f.close()
    #     g = open (infile, 'w')
    #     g.write (text)
    #     g.close ()
    #     print('fix completed')
    currrec = ""

    proplist = []
    ct = 0
    with open(outfile, 'w') as outf:
        outf.write ('smiles,idnumber\n')
        with open(infile, 'r', buffering=100000) as f:
            for line in f:
                if line.startswith('$$$$'):
                   currrec +=line
                   try:
                       currrec = currrec.replace('\n>','\n\n>')
                       sds = Chem.SDMolSupplier ()
                       sds.SetData(currrec)
                       mol = next(sds)
                       smi = Chem.MolToSmiles(mol)

                       prop_dict = mol.GetPropsAsDict()

                       for k in prop_dict.keys ():
                           if k not in proplist:
                               proplist.append (k)
                       outf.write (smi + ',' + prop_dict [idcol] + '\n')
                       ct += 1
                   except Exception as e:
                       print ('Error:'  ,  str(e))
                   currrec = ""
                else:
                   currrec += line
                outf.flush ()

                print ('line: ' + str(ct), end='\r')
        print ('completed')

def Canonicalize (smiles):
    RDLogger.DisableLog('rdApp.*')
    try:
        return Chem.CanonSmiles (smiles)
    except:
        return smiles

if __name__== "__main__":
    print ('HERE')
    sdfinfile = '/Users/eric/Enamine_screening_collection_202312.sdf'
    csvfile = '/Users/eric/Enamine_screening_collection_202312.csv'
    ssfile = '/Users/eric/Enamine_screening_collection_202312.ss.csv'
    outfile = '/Users/eric/Enamine_screening_collection_202312.ss.filt.csv'
    SDFtoFile (sdfinfile, True, csvfile)
    SaltStripFile(csvfile, 'smiles', ssfile)
    filterdict = {'MW_low': 275, 'MW_high': 650, 'overwrite': True}
    Filter_File(ssfile, outfile, ',', filterdict, None, False, False)