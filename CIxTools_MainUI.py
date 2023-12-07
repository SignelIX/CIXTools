import streamlit as st
from Archive import Enumerate_UI
import Chem_SpaceVisualization_UI
import Chem_CalcProps_UI
import Chem_Similarity_Diversity_UI
from streamlit import runtime


class CIXTools_MainUI:
    page = 'HOME'
    def __init__(self):
        st.set_page_config(layout="wide")
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
            ('ENUM', 'Enumeration', Enumerate_UI.EnumerationUI, None),
            ('CHEMSPACE', 'Display Chemical Space', Chem_SpaceVisualization_UI.Chem_SpaceVisualizationUI, None),
            ('ADD_PROPS', 'Add Chem Properties', Chem_CalcProps_UI.Chem_CalcPropsUI, None),
            ('SIMILARITY', 'Similarity Analysis', Chem_Similarity_Diversity_UI.Chem_Similarity_DiversityUI, 'Sim Search'),
            ('BBSIMILARITY', 'BB Similarity', Chem_Similarity_Diversity_UI.Chem_Similarity_DiversityUI, 'BB Similarity')
        ]
        for b in btnlist:
            btn  = st.sidebar.button(b[1])
            if btn == True:
                pg  = b[2]
                param = b[3]
                self.page = b[0]
                for bx in btnlist:
                    if bx[0] != b[0]:
                        st.session_state[bx[0]] = False
                    else:
                        st.session_state[bx[0]] = True

        if self.page == '':
            for b in btnlist:
                if (b[0] in st.session_state and st.session_state[b[0]] == True):
                    self.page = b[0]
                    pg  = b[2]
                    param = b[3]
                    st.session_state[b[0]] = True
                else:
                    st.session_state[b[0]] = False

        if self.page != '' and pg is not None:
            if param is not None:
                pg().RunUI(param)
            else:
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
    if runtime.exists():
        cim = CIXTools_MainUI()
        cim.RunUI()
