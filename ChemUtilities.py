from rdkit.Chem import SaltRemover
from rdkit import Chem
from rdkit.Chem import Descriptors as Descr
import matplotlib.pyplot as plt
from itertools import islice
from multiprocessing.pool import ThreadPool as Pool
from tqdm import tqdm
import pandas as pd
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import Enumerate
import pathlib

def SaltStripMolecules (molecules: pd.DataFrame):
    tqdm.pandas ()
    print ('salt stripping')
    print (molecules.columns)
    molecules['SMILES'] = molecules['SMILES'].progress_apply(lambda smi:SaltStrip(smi))
    print ('completed salt stripping')
    return molecules

def SaltStrip (molec):
    saltstrip = SaltRemover.SaltRemover()
    if type(molec) != Chem.Mol:
        m = Chem.MolFromSmiles(molec)
    else:
        m = molec
    m=  saltstrip.StripMol(m)
    smi = Chem.MolToSmiles(m, kekuleSmiles = False)
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

def Deprotect (compound, rxnscheme, iterate = False):
    deprotected = False
    Enumerator = Enumerate.Enumerate ()
    first_round = True

    has_scheme = Enumerator.ReadRxnScheme(rxnscheme, 'Deprotection', False)
    if has_scheme[0] is not None:
        last_dp = False
        while first_round or  (iterate and last_dp == True):
            res, prod_ct = Enumerator.RunRxnScheme([compound, 'None'],rxnscheme,'Deprotection', False)

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
                print (idx)
        except:
            continue

    df.to_csv (filename.replace ('.csv', '.props.csv'))

def PlotProperties (filename, cols, bins):
    df = pd.read_csv(filename)
    for c in cols:
        df.hist (column = c, bins = bins [c], rwidth=0.9)
        plt.show ()

def Filters ( m, filter_dict, stats):
    if stats is None:
        stats = {'MW':0,'tPSA':0, 'logP': 0 }
    if filter_dict is None:
        return False, stats;
    if ('MW_high' in filter_dict) or ('MW_low' in filter_dict):
        mw = Descr.MolWt(m)
    if 'MW_high' in filter_dict:
        if mw > filter_dict['MW_high']:
            stats ['MW'] += 1
            return True, stats
    if 'MW_low' in filter_dict:
        if mw < filter_dict['MW_low']:
            stats ['MW'] += 1
            return True, stats
    if ('tPSA' in filter_dict):
        tpsa = Descr.TPSA(m)
        if tpsa > filter_dict['tPSA']:
            stats ['tPSA'] += 1
            return True, stats
    if ('logP_high' in filter_dict) or ('logP_low' in filter_dict):
        logP = Descr.MolLogP(m)
    if ('logP_high' in  filter_dict):
        logphigh = filter_dict['logP_high']
        if logP > logphigh: #aLogP
            # print ('logphigh', logP)
            stats ['logP'] += 1
            return True, stats
    if ('logP_low' in filter_dict):
        logplow = filter_dict['logP_low']
        if logP < logplow: #aLogP
            stats ['logP'] += 1
            return True, stats

    return False, stats


def ApplyFilters ( line, smilescol, filter_dict, splitchar, ssfilters, useChirality, AmbiguousChirality, Keep):
    smi = line.split (splitchar) [smilescol]
    try:
        m = Chem.MolFromSmiles(Kekulize_Smiles( smi))
    except:
        return Keep
    propFiltered = Filters(m , filter_dict)
    ssFiltered = Substructure_Filters (m, None, ssfilters, useChirality, Keep)
    ambig = False
    if (AmbiguousChirality == True):
        ambig = ContainsUnresolvedChirality(m)

    if (Keep):
        return propFiltered and ssFiltered and ambig
    else:
        return propFiltered or ssFiltered or ambig


def Filter_File (catfile, outfilename, splitchar, filter_dict, ss_file, useChirality, AmbigChirality, Keep = False):
    print ('Filtering')
    N = 10000
    pool_size =  8
    subct = 0
    hdrread = False
    lct = 0

    outfile = open (outfilename, 'w')
    ct = 0
    ssfilters = Read_FilterFile(ss_file)

    with open (catfile) as collection:
        for mlist in iter(lambda: list(islice(collection, N)), []):
            pool = Pool(pool_size)
            block = ''
            for line in mlist:
                if line == '\n':
                    continue
                if hdrread == True:
                    isFiltered = pool.apply_async(ApplyFilters, args= (line, smilescol, filter_dict, splitchar, ssfilters,useChirality, AmbigChirality, Keep)).get()
                    if (Keep == True and isFiltered == True) or (Keep == False and isFiltered == False):
                        subct += 1
                        block += line.strip ().replace(' ',',') + '\n'
                else:
                    hdrlist = line.split (splitchar)
                    hdrlist = [sx.upper() for sx in hdrlist]
                    matching = [s for s in hdrlist if "SMILES" in s]
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

def Read_FilterFile (filter_file):
    if filter_file is None:
        return None
    ss_filters = pd.read_csv(filter_file)
    ss_filters = ss_filters.assign(SS_ct=0, Molecule=None)
    for idx, row in ss_filters.iterrows():
        mol = Chem.MolFromSmarts(row['SMILES'])
        ss_filters.loc[idx, 'Molecule'] = mol
    return ss_filters

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

def RemoveChirality (infile, outfile, smilescol):
    df = pd.read_csv(infile)
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

def SDFtoDF (infile, fix):
    if fix ==True:
        f=open (infile, 'rb')
        text = f.read().decode(errors='replace')
        f.close()
        g = open (infile, 'w')
        g.write (text)
        g.close ()

    sdf = Chem.SDMolSupplier (infile)
    df = pd.DataFrame ()

    for sdrec in sdf:
        props = sdrec.GetPropsAsDict ()
        rowidx = len(df.index)

        for p in props:
            if p not in df:
                df [p]= None
        df.loc[len(df)] = [None] * len(df.columns)
        for p in props:
            df.loc[rowidx][p] = props[p]

    return df

def ConvertSDFilesFromDir (inpath, outfile):
    globlist = pathlib.Path(inpath).glob('*.sdf')
    alldf = pd.DataFrame ()
    for f in globlist:
        print (f)
        df = SDFtoDF(str(f), True)
        alldf = pd.concat([alldf, df], axis=0, ignore_index=True)

    alldf.to_csv (outfile)

