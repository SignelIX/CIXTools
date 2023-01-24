import Chem_Similarity_Diversity
import pandas as pd
import NamingStandardizer
import matplotlib.pyplot as plt
import umap
import pickle
import seaborn as sns
import Enumeration_CombinatorialLibrary as cle
import numpy as np
import matplotlib
import pathlib


class Chem_SpaceVisualization:
    clrs = ['red','blue','purple','black', 'orange', 'aqua', 'darkgreen', 'darkblue']
    div  = Chem_Similarity_Diversity.Diversity ()
    def ChemicalSpace_UMap(self, outpath, outprefix, files_list, clr_list, in_pkl_file, numpts, underlying_mapfile ):
        fig_file = outpath + outprefix + '.png'
        out_file = outpath + outprefix + '.csv'
        if underlying_mapfile != None:
            umdf = pd.read_csv (underlying_mapfile)
        libdf = None
        legend_labels = {}
        if clr_list is None:
            clr_list = self.clrs

        for ix in range(0, len(files_list)):
            l = files_list[ix]
            if ix  >= len(clr_list):
                clr = 'gray'
            else:
                clr = clr_list[ix]
            if libdf is None:
                libdf = pd.read_csv(l)
                libdf['Color'] = clr
                if numpts == -1:
                    numpts = len (libdf)
                libdf = libdf.sample(frac=min(1, numpts / len(libdf)))
            else:
                df = pd.read_csv(l)
                df['Color'] = clr
                if numpts == -1:
                    numpts = len (libdf)
                libdf = libdf.append(df.sample(frac=min(1,numpts / len(df))))
            fsplit = files_list[ix].split('/')
            fname = fsplit [len (fsplit) - 1]
            legend_labels [fname] = clr

        if underlying_mapfile is not None:
            legend_labels[pathlib.Path(underlying_mapfile).name ] = 'gray'
        libdf = NamingStandardizer.Standardize_Colnames(libdf, 'SMILES' )

        libdf = libdf.sample(frac=1.0)

        smileslist = list(libdf['SMILES'])
        colorlist = list(libdf['Color'])
        print('generating umap/coordinates')

        pltlyfig = self.Generate_UMAP(None, smileslist, outfname=out_file, fig_fname=fig_file, modelfile=in_pkl_file,
                          transformcolorlist=colorlist, legendlabels=legend_labels, underlyingmaplist= umdf.values.tolist())
        return fig_file, pltlyfig

    def CreateMultiFile_UMap (self, outpath, outprefix, infilelist, frac):
        out_img_file = outpath + outprefix + '.png'
        out_pkl_file = outpath + outprefix + '.pkl'
        out_csv_file = outpath + outprefix + '.csv'
        df = None
        legendlabels = {}
        colorlist = []
        clrct = 0
        for infile in infilelist:
            if df is None:
                df = pd.read_csv(infile, sep=',')
                colorlist.extend([self.clrs[clrct]] * len(df))
            else:
                tempdf = pd.read_csv(infile, sep=',')
                colorlist.extend([self.clrs[clrct]] * len(tempdf))
                df = df.append ( tempdf)

            legendlabels[pathlib.Path (infile).name] = self.clrs[clrct]

            clrct +=1


        df = NamingStandardizer.Standardize_Colnames(df, 'SMILES')
        df ['color'] = colorlist
        if df is not None:
            df = df.sample(frac=frac)

        plotlyfig = self.Generate_UMAP(list(df['SMILES']), None, fig_fname=out_img_file, outfname=out_csv_file,
                            modeloutfile=out_pkl_file, fit_colorlist=df['color'], legendlabels=legendlabels)
        return out_pkl_file, out_img_file, out_csv_file, plotlyfig

    def Generate_UMAP (self, fit_cmpdlist, transform_cmpdlist, fig_fname, outfname, modelfile = None,
                       transformcolorlist = None, fit_colorlist = None, modeloutfile=None, legendlabels = None, underlyingmaplist=None ):

        if outfname is not None:
            print(outfname)
            outf = open(outfname, 'w')
            outf.write ('SMILES,x,y,color\n')

        colorlist = []
        alphalist = []
        sizelist = []

        if fit_cmpdlist is not None:
            if fit_colorlist is not None:
                print (len(fit_colorlist), len (fit_cmpdlist))
                colorlist.extend (fit_colorlist)
            else:
                colorlist = ['gray']*(len(fit_cmpdlist))
            alphalist = [1.0]  *(len(fit_cmpdlist))
            sizelist = [4.0]*(len(fit_cmpdlist))
        if underlyingmaplist is not None:
            colorlist.extend(['gray'] * len(underlyingmaplist))
            alphalist = [1.0] * (len(underlyingmaplist))
            sizelist = [4.0] * (len(underlyingmaplist))

        if transform_cmpdlist is not None:
            if transformcolorlist is not None:
                colorlist.extend(transformcolorlist)
                for cx in transformcolorlist:
                    if cx == 'red':
                        alphalist.append (.2)
                        sizelist.append (4.0)
                    else:
                        alphalist.append(.7)
                        sizelist.append(8.0)
            else:
                colorlist.extend(['red'] * len(transform_cmpdlist))
                alphalist.extend([.2] * len(transform_cmpdlist))
                sizelist.extend([4.0] *len(transform_cmpdlist))

        all_cmpds_list = []
        if fit_cmpdlist is not None:
            all_cmpds_list.extend (fit_cmpdlist)

        if transform_cmpdlist is not None:
            all_cmpds_list.extend(transform_cmpdlist)

        fitdata, deletefps = self.div.Get_FPArrayAsVstack(all_cmpds_list)
        if fitdata is None:
            return
        deletefps = sorted (deletefps, reverse=True)
        if len (deletefps) > 0:
            for d in  deletefps:
                colorlist.pop (d)
                all_cmpds_list.pop (d)
                alphalist.pop (d)
                sizelist.pop(d)
        sns.set(style='white', context='poster', rc={'figure.figsize': (14, 10)})

        print ('modelfile', modelfile)
        if modelfile is None:
            print ('calculating UMAP')
            fit= umap.UMAP(metric = "jaccard",
                          n_neighbors = 25,
                          n_components = 2,
                          low_memory = False,
                          min_dist = 0.001,
                          verbose = True)
            print('calculation complete. Fitting...')
            u = fit.fit_transform(fitdata)
            print('umap fitting complete')
        else:
            f = open (modelfile, 'rb')
            fit = pickle.load (f)
            f.close ()
            u = fit.transform(fitdata)

        if (len(all_cmpds_list) != len(u)):
            print('array mismatch')
            return

        if outfname is not None:
            for i in range (0, len (all_cmpds_list)):
                outf.write(','.join(str(a) for a in [all_cmpds_list[i], u[i, 0], u[i, 1], colorlist [i], '\n']))
            outf.close ()


        fitdata = []
        um_pts = []

        if underlyingmaplist is not None:
            for uix in underlyingmaplist:
                um_pts.append ([uix[0],uix[1]])

            if u is None:
                u = np.array(um_pts)
            else:
                u = np.append( np.array(um_pts), u,0)

        plt.scatter(u[:, 0], u[:, 1], color = colorlist, marker='o', s=sizelist, alpha= alphalist,  edgecolors='k', linewidths=0.2)

        pltlyfig = None

        if transform_cmpdlist is not None:
            len_transform = str(len(transform_cmpdlist))
        else:
            len_transform = 'N/A'
        plt.title('UMAP embedding of Diverse Set: fit points:' +  str(len(fitdata)) + ' diversity points:' +  len_transform )
        patches = []

        if legendlabels is not None:
            for cx in legendlabels.keys():
                print ('LABELS:', cx)
                patches.append (matplotlib.patches.Patch(color=legendlabels[cx], label=cx))
        plt.legend(loc='upper right', handles = patches)
        plt.savefig(fig_fname)
        if fig_fname.endswith ('.png'):
            svgname = fig_fname.replace ('.png', '.svg')
        else:
            svgname = fig_fname + '.svg'

        plt.savefig(svgname, format='svg')
        plt.show ()

        if modeloutfile is not None:
            print('writing model')
            pickle.dump(fit, open(modeloutfile, 'wb'))
            f = open(modeloutfile.replace('.pkl', '.mdl.csv'), 'w')
            for ux in u:
                f.write(str(ux[0]) + ',' + str(ux[1]) + '\n')
            f.close()
        return pltlyfig

    def Generate_UMAP_From_File (self, fullfile, cmpdfile, outfile, smilescol, colorcol, umap_outfile, numpts = -1):
        infile = open(fullfile)
        dffull = pd.read_csv(infile)
        infile.close()

        infile = open(cmpdfile)
        dfset = pd.read_csv(infile)
        infile.close()

        if numpts == -1:
            fulllist = dffull[smilescol].tolist()
            cmpdlist = dfset[smilescol].tolist()
            colorlist = dfset[colorcol].tolist()
        else:
            fulllist = dffull[smilescol].tolist()[0:numpts]
            cmpdlist = dfset[smilescol].tolist()[0:numpts]
            colorlist = dfset[colorcol].tolist()[0:numpts]


        self.Generate_UMAP(fulllist, cmpdlist, outfile, umap_outfile, colorlist  )

    def Generate_UmapEmbedding(self, libname, bbdict, ncomp,
                               rxnschemefile, modeloutfile, umap_maxptct):
        cycs = list(bbdict.keys())

        print('building umap')
        umapsetdict = {}
        umapstruclist = []

        for cyc in bbdict.keys():
            bbdict[cyc] = bbdict[cyc].sample(frac=1.0)
            umapsetdict[cyc] = []
            for i in range(0, min(umap_maxptct, len(bbdict[cyc]))):
                umapsetdict[cyc].append(i)

        tuples = cle.recurs_tuples(umapsetdict, None, cycs, cycs[0])
        tupleinfodict = {}

        umapstruclist = cle.EnumTuples_dask(rxnschemefile, libname, tuples, bbdict, cycs, tupleinfodict)
        struclist = []
        for ux in umapstruclist:
            struclist.append(ux[0])
        fplist = self.generatefps_dask(struclist, True, bitvec=True)
        fitdata = np.vstack(fplist)

        print('generating umap...')
        fit = umap.UMAP(metric="jaccard",
                        n_components=ncomp,
                        n_neighbors=25,
                        low_memory=False,
                        min_dist=0.001,
                        verbose=True)
        u = fit.fit_transform(fitdata)
        print('writing umap to pkl')

        if modeloutfile is not None:
            pickle.dump(fit, open(modeloutfile, 'wb'))

        colorlist = ['b'] * len(u)
        plt.scatter(u[:, 0], u[:, 1], s=9, c=colorlist, label='test', edgecolors='k', linewidths=0.2)
        plt.show()
        if modeloutfile is not None:
            print('writing umap data')
            f = open(modeloutfile.replace('.pkl', '.csv'), 'w')
            for ix in range(0, len(u)):
                f.write(umapstruclist[ix][0] + ',' + str(u[ix][0]) + ',' + str(u[ix][1]) + '\n')
            f.close()
