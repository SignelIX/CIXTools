import Chem_Similarity_Diversity
import pandas as pd
import streamlit as st
import json
import NamingStandardizer

class Chem_SpaceVisualization:
    initpath = '../CIxTools.init.json'
    div = Chem_Similarity_Diversity.Diversity()
    paramslist = ['chemspace_outpath', 'chemspace_outprefix', 'chemspace_infilelist', 'chemspace_pklfile']

    def __init__ (self):
        f = open(self.initpath)
        initjson = json.load(f)
        f.close()
        for p in self.paramslist:
            if p not in st.session_state or st.session_state[p] == None or st.session_state [p] == '' :
                if p in initjson:
                    st.session_state[p] = initjson [p]
                else:
                    st.session_state[p] = ''

    def ChemicalSpace_UMap(self, outpath, outprefix, files_list, clr_list, in_pkl_file ):
        fig_file = outpath + outprefix + '.png'
        out_file = outpath + outprefix + '.csv'

        libdf = None
        legend_labels = {}
        for ix in range(0, len(files_list)):
            legend_labels = {}
            l = files_list[ix]
            if ix  >= len(clr_list):
                clr = 'gray'
            else:
                clr = clr_list[ix]
            if libdf is None:
                libdf = pd.read_csv(l)
                libdf['Color'] = clr
                libdf = libdf.sample(frac=5000 / len(libdf))
            else:
                df = pd.read_csv(l)
                df['Color'] = clr
                libdf = libdf.append(df.sample(frac=5000 / len(df)))
            fsplit = files_list[ix].split('/')
            fname = fsplit [len (fsplit) - 1]
            legend_labels [fname] = clr

        libdf = NamingStandardizer.Standardize_Colnames(libdf, 'SMILES' )

        libdf = libdf.sample(frac=1.0)

        smileslist = list(libdf['SMILES'])
        colorlist = list(libdf['Color'])
        print('generating umap/coordinates')

        self.div.Generate_UMAP(None, smileslist, outfname=out_file, fig_fname=fig_file, modelfile=in_pkl_file,
                          transformcolorlist=colorlist, legendlabels=legend_labels)
        return fig_file

    def Generate_UMap (self, outpath, outprefix, infilelist):
        out_img_file = outpath + outprefix + '.png'
        out_pkl_file = outpath + outprefix + '.pkl'
        out_csv_file = outpath + outprefix + '.csv'
        df = None
        for infile in infilelist:
            if df is None:
                df = pd.read_csv(infile, sep=',')
            else:
                df = df.append ( pd.read_csv(infile, sep=','))
        df = NamingStandardizer.Standardize_Colnames(df, 'SMILES')
        if df is not None:
            df = df.sample(frac=.3)
        self.div.Generate_UMAP(list(df['full_smiles']), None, fig_fname=out_img_file, outfname=out_csv_file,
                          modeloutfile=out_pkl_file, transformcolorlist=None)
        return out_pkl_file, out_img_file, out_csv_file

    def streamlit (self):
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
            Umap Chemspace Visualization</h1>""", unsafe_allow_html=True)

        outpath = st.text_input (label ='chemspace_outpath', key='chemspace_outpath', on_change=self.SaveToInit)
        outprefix = st.text_input (label='chemspace_outprefix', key='chemspace_outprefix',on_change=self.SaveToInit)
        infilelist = st.text_area (label='chemspace_infilelist', key='chemspace_infilelist', on_change=self.SaveToInit).split (',')
        pkfile = st.text_input(label='chemspace_pklfile', key='chemspace_pklfile', on_change=self.SaveToInit )
        if st.button(label='run umap generation', key='RunUMAP'):
            pklfile, imgfile, csvfile  = self.Generate_UMap (outpath, outprefix, infilelist)
            st.text (pklfile)
            st.image (imgfile)
        if st.button(label='run umap vs model', key='RunUMAPvModel'):
            with st.spinner('Building Map...'):
                fig_file = self.ChemicalSpace_UMap(outpath, outprefix, infilelist, ['red', 'blue', 'aqua', 'purple'],
                                                   pkfile)
                st.image(fig_file)
            st.success('Done!')

    def SaveToInit(self):

        with open(self.initpath, "r") as jsonFile:
            data = json.load(jsonFile)
        for p in self.paramslist:
            if p in st.session_state:
                data[p] = st.session_state[p]
        with open(self.initpath, "w") as jsonFile:
            print (json.dumps(data, indent=4))
            jsonFile.write(json.dumps(data, indent=4))

    def RunUI (self):
        self.streamlit ()

if __name__ == "__main__":
    csvis = Chem_SpaceVisualization ()
    csvis.RunUI()