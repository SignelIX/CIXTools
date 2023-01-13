from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
import numpy as np
from PIL import Image
import base64
import streamlit as st
from st_aggrid import AgGrid, JsCode
from st_aggrid.grid_options_builder import  GridOptionsBuilder
from st_aggrid import GridUpdateMode, DataReturnMode
import sys
from streamlit import runtime
import pandas as pd


def ShowMol( smiles, outtype=None):
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
    if outtype == 'b64_datauri':
        buf = io.BytesIO()
        img.save(buf, format='png')
        buf.seek(0)
        data = buf.read()
        sb = ''
        sb += "data:image/png;base64,"
        sb += base64.b64encode(data).decode('ascii')
        return sb
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
        elif outtype == 'b64_datauri':
            buf = io.BytesIO()
            img.save(buf, format='png')
            buf.seek(0)
            data = buf.read()
            sb = ''
            sb += "data:image/png;base64,"
            sb += base64.b64encode(data).decode('ascii')
            return sb
        else:
            img.save('test_mols.png')
            plt.imshow(img)
            return plt
    else:
        return None
    
def ShowMols_StreamLit_Grid (df, smilescols = ['SMILES'], rowheight = 100):
    dispdf = df
    for s in smilescols:
        dispdf [s + '_Image'] = None
    image_render = JsCode("""function (params) {
                         this.params = params;
                         this.img_render = document.createElement('div');
                         this.img_render.innerHTML = `
                             <img src=${this.params.value}
                                 width=""" + str(rowheight) +"""
                                 height=""" + str(rowheight) +"""
                             >
                         `; 
                         return this.img_render;
                         }""")
    gb = GridOptionsBuilder()
    print(smilescols)
    for c in df.columns:
        if c not in smilescols:
            gb.configure_column(c, headerName=c, width=50)
        else:
            print (c)
            gb.configure_column(c + "_Image", headerName=s, cellRenderer=image_render)
    for ix, row in dispdf.iterrows():
        for s in smilescols:
            print (row[s])
            if row[s] is not None:
                 dispdf.at [ix,s + '_Image'] =ShowMol(row[s],outtype='b64_datauri')
            else:
                 dispdf.at[ix, s + '_Image'] = None

    gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, editable=True,
                                enableRangeSelection=True, )

    gridOptions = gb.build()
    gridOptions['rowHeight'] = rowheight

    AgGrid(dispdf, allow_unsafe_jscode=True, enable_enterprise_modules=True, gridOptions=gridOptions, height=500)
    return

class Chem_ShowMols_UI:
    def body(self):
        df = pd.DataFrame ([['CCCC', 'c1ccccc1', 5], ['CCNC1CCCC1', 'C1=NCC=N1', 10]], columns=['SMILES','SMILES2',  'Test'])
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
               MolDisplay</h1>""", unsafe_allow_html=True)

        ShowMols_StreamLit_Grid (df, ['SMILES','SMILES2'])


    def RunUI(self):
        self.body()

if __name__ == "__main__":
    if runtime.exists():
        csm = Chem_ShowMols_UI()
        csm.RunUI()