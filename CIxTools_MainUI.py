import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import Enumerate
import Chem_SpaceVisualization
import streamlit.components.v1 as components

class CIXTools_MainUI:
    page = 'HOME'
    def __init__(self):
        st.session_state['MAIN'] = True
    def head(self):
        st.markdown("""
              <h1 style='text-align: center; margin-bottom: -35px;'>
              CIXTools
              </h1>
              
          """, unsafe_allow_html=True
                    )
    def Get_Page (self):
        self.page = ''
        pg = None
        btnlist = [
            ('ENUM', 'Enumeration', Enumerate.EnumerationUI),
            ('CHEMSPACE', 'Display Chemical Space', Chem_SpaceVisualization.Chem_SpaceVisualization)
        ]
        for b in btnlist:
            btn  = st.sidebar.button(b[1])
            if btn == True:
                pg  = b[2]
                self.page = b[0]
                for bx in btnlist:
                    if bx[0] != b[0]:
                        st.session_state[bx[0]] = False
                    else:
                        st.session_state[bx[0]] = True

        if self.page == '':
            for b in btnlist:
                if (b[0] in st.session_state and st.session_state[b[0]] == True):
                    print (b[0], st.session_state[b[0]])
                    self.page = b[0]
                    pg  = b[2]
                    st.session_state[b[0]] = True
                else:
                    st.session_state[b[0]] = False

        if self.page != '' and pg is not None:
            pg ().RunUI()
        elif self.page == '':
            st.write('Welcome to CIxTools, select a button from the side bar to begin')
        else:
            st.write('not implemented')


    def body (self):
        st.session_state['MAIN'] = False
        st.sidebar.markdown("## Controls")
        self.Get_Page ()


    def RunUI (self):
        self.head ()
        self.body ()



if __name__=="__main__":
    if st._is_running_with_streamlit:
        cim = CIXTools_MainUI()
        cim.RunUI()
