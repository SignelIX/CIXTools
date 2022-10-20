from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
import numpy as np
from PIL import Image


def ShowMol( smiles):
    if smiles is not None:
        try:
            m = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(m, size=(400, 400))
        except Exception as eq:
            img = np.zeros([100, 100, 3], dtype=np.uint8)
            img.fill(255)  # or img[:] = 255
            img = Image.fromarray(img)
    else:
        img = np.zeros([100, 100, 3], dtype=np.uint8)
        img.fill(255)  # or img[:] = 255
    img.save('test.png')
    plt.imshow(img)
    plt.show()
    return plt

def ShowMols( smiles_list, cols = 4, subImgSize = (200,200), outtype = 'plt'):
    list = []
    if smiles_list is not None:
        for smi in smiles_list:
            if smi is not None and smi != 'None' and smi != 'NOSMILES':
                if not (isinstance(smi, Chem.rdchem.Mol)):
                    try:
                        m = Chem.MolFromSmiles(smi)
                    except:
                        m = Chem.MolFromSmiles ('')
                    list.append(m)
                else:
                    list.append (smi)
            img = Draw.MolsToGridImage(list, molsPerRow=cols, subImgSize=subImgSize )

            if outtype == 'bytesio':
                buf = io.BytesIO()
                img.save(buf, format='png')
                buf.seek(0)
                return buf.read()
            else:
                img.save('test_mols.png')
                plt.imshow(img)
                plt.show()
                return plt
    else:
        return None
