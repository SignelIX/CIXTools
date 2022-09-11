import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
#import Enumerate

class CIXTools_MainUI:
    def __init__(self):
        st.session_state['MAIN'] = True
    def head(self):
        st.markdown("""
              <h1 style='text-align: center; margin-bottom: -35px;'>
              CIXTools
              </h1>
          """, unsafe_allow_html=True
                    )
    def body (self):
        st.session_state['MAIN'] = False
        import streamlit as st
        from rdkit import Chem
        from rdkit.Chem import Draw
        from PIL import Image

        compound_smiles = 'c1cc(C(=O)O)c(OC(=O)C)cc1'
        m = Chem.MolFromSmiles(compound_smiles)
        im = Draw.MolToImage(m)

        st.image(im)
        # if 'ENUM' in st.session_state and st.session_state['ENUM'] == True:
        #     enum = Enumerate.EnumerationUI()
        #     enum.RunUI('SMILES', '')
        # elif 'CHEMSPACE' in st.session_state and st.session_state['CHEMSPACE'] == True:
        #     st.write ('not implemented')
        # else:
        #     st.session_state['MAIN'] = True
        #     testenum = st.button(label='Test Enumerator', key='ENUM')
        #     if testenum:
        #         st.session_state['Page'] = 'ENUM'
        #     chemspace = st.button(label='Display Chemical Space', key='CHEMSPACE')
        #     if chemspace:
        #         st.session_state['Page'] = 'CHEMSPACE'

    def RunUI (self):
        self.head ()
        self.body ()


with st.echo(code_location='below'):
    st.text ('Hello')

# if __name__=="__main__":
#     if st._is_running_with_streamlit:
#         cim = CIXTools_MainUI()
#         cim.RunUI()
