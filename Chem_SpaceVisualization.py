import Chem_Similarity_Diversity
import pandas as pd
import streamlit as st

class Chem_SpaceVisualization:
    div = Chem_Similarity_Diversity.Diversity()

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
            print(len(libdf))
        libdf = libdf.sample(frac=1.0)
        print(len(libdf))

        smileslist = list(libdf['full_smiles'])
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
        if df is not None:
            df = df.sample(frac=.3)
        self.div.Generate_UMAP(list(df['full_smiles']), None, fig_fname=out_img_file, outfname=out_csv_file,
                          modeloutfile=out_pkl_file, transformcolorlist=None)
        return out_pkl_file, out_img_file, out_csv_file

    def streamlit (self):
        st.markdown("""<h1 style='text-align: center; margin-bottom: -35px;'>
            Umap Chemspace Visualization</h1>""", unsafe_allow_html=True)
        outpath = st.text_input (label ='outpath', key='outpath')
        outprefix = st.text_input (label='outprefix', key='outprefix')
        infilelist = st.text_area (label='infilelist', key='infilelist').split (',')
        pkfile = st.text_input(label='pklfile', key='pklfile')
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

    def RunUI (self):
        self.streamlit ()

if __name__ == "__main__":
    csvis = Chem_SpaceVisualization ()
    csvis.RunUI()