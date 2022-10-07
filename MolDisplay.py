from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt


def ShowMol( smiles):
    m = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(m, size=(400, 400))
    img.save('test.png')
    plt.imshow(img)
    plt.show()
    return plt

def ShowMols( smiles_list, cols = 4):
    print (smiles_list)
    list = []
    for smi in smiles_list:
        if smi is not None and smi != 'None' and smi != 'NOSMILES':
            if not (isinstance(smi, Chem.rdchem.Mol)):
                list.append(Chem.MolFromSmiles(smi))
            else:
                list.append (smi)
    print (list)
    img = Draw.MolsToGridImage(list, molsPerRow=cols, subImgSize=(200, 200) )
    img.save('test_mols.png')
    plt.imshow(img)
    plt.show()
    return plt
