from rdkit.Chem import AllChem as Chem
from rdkit.Chem import SDMolSupplier
from tqdm import tqdm

class Format_Conversion:
    def sdf_to_smiles (self, infile, outfile, n, addsmiles = True):
        print('loading')
        suppl = SDMolSupplier(infile)

        print('processing')
        outf = open (outfile, 'w')
        ct = 0
        smi_list  = []
        proplist = []
        print ('N:', n)

        if addsmiles:
            proplist.append ('SMILES')
        for mol in tqdm(suppl):
            for p in list(mol.GetPropNames ()):
                if p not in proplist:
                    proplist.append (p)
            dict = mol.GetPropsAsDict()
            if addsmiles:
                smi = Chem.MolToSmiles(mol)
                dict ['SMILES'] = smi
            smi_list.append( dict)
            ct = ct + 1
            if n != -1 and ct == n:
                break

        print ('smiles:' , smi_list)

        outf.write ('id')

        for p in proplist:
            outf.write(',' + p)
        outf.write ('\n')
        ct = 0


        for s in smi_list:
            outf.write ( str(ct) )
            for p in proplist:
                if p in s:
                    val = str (s[p])
                    if ',' in val:
                        val = '\"' + val + '\"'
                    outf.write (',' + val)
                else:
                    outf.write(',')
            outf.write ('\n')
            ct = ct + 1
        outf.close ()

    def smiles_to_sdf (self, SMILES, outfile):
        mol = Chem.MolFromSmiles(SMILES)
        Chem.MolToMolFile(mol,outfile)



