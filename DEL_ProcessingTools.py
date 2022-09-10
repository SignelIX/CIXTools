import glob
import pandas as pd
import pathlib
from tqdm import  tqdm
import os
import dask.dataframe as dd
from timeit import default_timer as timer
from dask.diagnostics import ProgressBar

import math
class DELProcessingTools:
    NUM_WORKERS = 100
    def Load_BBDict (self, inpath):
        globlist = glob.glob(inpath + '/*_smiles.txt')

        BBdict = {}
        for f in globlist:
            df = pd.read_csv(f, sep='\t', header=None)
            df = df.rename(columns={0: "SMILES", 1: "BB_ID"})
            parsef = f.split ('/')
            lib = parsef [len(parsef) - 1] .split('_')[0]
            df[['cycle', 'BBnum']] = df['BB_ID'].str.split('_', 1, expand=True)
            BBdict[lib] = df
        return BBdict

    def AddStructuresToAlphaMaFiles_Dask (self, inpath, outpath, bb_dict):
        tqdm.pandas()

        def ProcessFile (filename, outfilename):
            def taskfcn(row, bbdf, bbcols):
                smilist = []
                for c in bbcols:
                    bbid = row [c]
                    cyc = c.replace('bb','')
                    bbid = str(cyc) + '_' + str(bbid)
                    val = bbdf[bbdf.BB_ID == bbid].iloc[0]['SMILES']
                    smilist.append(val)
                return smilist


            resdf = pd.read_csv (filename, sep='\t')
            dtype_dict = {}
            meta_dict = {}
            ct = 0
            for c in resdf.columns:
                if c.startswith ('bb'):
                    dtype_dict [c] = str
                    meta_dict[ct] = str
                    ct += 1
            resdf = resdf.astype(dtype_dict)

            lib  = filename.name.split('_')[0]
            bbdf = bb_dict [lib]
            bbcols = resdf.columns[resdf.columns.str.contains ('bb')].tolist()

            start = timer()
            ddf = dd.from_pandas(resdf, npartitions=self.NUM_WORKERS)
            res = ddf.apply(taskfcn, axis=1, result_type='expand', args = (bbdf, bbcols) , meta = meta_dict).compute ()
            print (res)
            res.columns =  [b + '_SMILES' for b in bbcols]
            df = resdf.merge(res, left_index=True, right_index=True)
            df.to_csv (outfilename, index=False)
            end = timer()
            print(end - start)

        print (inpath)
        globlist = pathlib.Path(inpath).rglob('*.txt')

        for f in tqdm(globlist):
            if f.name.startswith('AMA'):
                fldr  = f.parent.resolve()
                outfldr= str(fldr).replace(inpath, outpath)
                if not pathlib.Path(outfldr).exists():
                      os.makedirs(outfldr, exist_ok=True)
                outfilename = outfldr + '/' + f.name.replace('.txt', '.strx.csv')
                ProcessFile (f, outfilename)

    def AddStructuresToAlphaMaFiles(self, inpath, outpath, bb_dict):
        tqdm.pandas()

        def ProcessFile(filename, outfilename):
            resdf = pd.read_csv(filename, sep='\t')
            lib = filename.name.split('_')[0]
            bbdf = bb_dict[lib]
            bbcols = resdf.columns[resdf.columns.str.contains('bb')].tolist()
            bbsmilescols = [b + '_SMILES' for b in bbcols]
            resdf[bbsmilescols] = [''] * len(bbcols)

            for idx, r in tqdm(resdf.iterrows(), total=len(resdf)):
                for c in bbcols:
                    bbid = r[c]
                    cyc = c.replace('bb', '')
                    bbid = str(cyc) + '_' + str(bbid)
                    val = bbdf[bbdf.BB_ID == bbid].iloc[0]['SMILES']
                    resdf.loc[idx, c + '_SMILES'] = str(val)

            resdf.to_csv(outfilename, index=False)

        globlist = pathlib.Path(inpath).rglob('*.txt')

        for f in globlist:
            if f.name.startswith('AMA'):
                fldr = f.parent.resolve()
                outfldr = str(fldr).replace(inpath, outpath)
                if not pathlib.Path(outfldr).exists():
                    os.makedirs(outfldr)
                outfilename = outfldr + '/' + f.name.replace('.txt', '.strx.csv')
                ProcessFile(f, outfilename)

        return

    def CalculateAlphaMa_Disynthon_Counts(self, library, libcond, infile, synthstats):
        def Make_CountDict (df, syn):
            ct_dict = {}
            resdf = df.groupby(syn)['Counts'].agg('sum')
            resdict = resdf.to_dict()
            for s in resdict:
                ct_dict[s] = {}
                k_counts = resdict[s]
                ct_dict[s]['ct'] = k_counts
            return ct_dict

        df= pd.read_csv (infile, sep = '\t', header= None)

        colidx = 1
        synthstats [library][libcond] = {}

        cols = {}
        synthlist = []
        bbcols = []
        for c in df.columns[0:len(df.columns)]:
            if c < len(df.columns) - 1:
                cols[c] = 'BB' + str(c  +1)
                synthlist.append ([cols[c]])
                bbcols.append(cols[c])
            else:
                cols [c] = 'Counts'
        df= df.rename (columns = cols)

        for b1 in range (0, len (bbcols)-1):
            for b2 in range (b1+1, len (bbcols)):
                synthlist.append ([bbcols[b1], bbcols[b2]])

        for syn in synthlist:
            synthstats [library][libcond] ["|".join(syn)] = Make_CountDict (df, syn)


        synthstats[library][libcond]['total'] =df['Counts'].agg ('sum')
        return synthstats

    def Get_TagStats (self, lib,  libtagfile, tagdict):
        if tagdict is None:
            tagdict = {}
        df = pd.read_csv(libtagfile, sep='\t', header=None)
        tagdict [lib] = {}
        cols = {}
        cols [0] = 'Cycle'
        cols [1] = 'BBnum'
        for c in df.columns[2:]:
            cols [c] = 'tag' + str(c - 1)
        df = df.rename(columns=cols)
        for idx, row in df.iterrows ():
            cyc = 'BB' + str(row['Cycle'])
            if cyc not in tagdict[lib]:
                tagdict [lib][cyc]={}
                tagdict[lib][cyc]['BBs'] = {}
                tagdict[lib][cyc]['cycle_tagct'] = 0
            if row ['BBnum'] not in tagdict[lib][cyc]['BBs']:
                tagdict[lib][cyc]['BBs'][row['BBnum']]= {}
                tagdict[lib][cyc]['BBs'][row['BBnum']]['tags']= []
                tagdict[lib][cyc]['BBs'][row['BBnum']]['tagct'] = 0
            for  c in df.columns:
                if c.startswith ('tag'):
                    tagdict[lib][cyc]['BBs'][row['BBnum']]['tags'].append(row[c])
                    tagdict[lib][cyc]['BBs'][row['BBnum']]['tagct'] += 1
                    tagdict[lib][cyc]['cycle_tagct'] += 1
        return  tagdict

    def Add_EnrichmentsToSynthStats (self, liblist, tagdict, synthstats):

        ct = 0
        for lib in liblist:
            ct += 1
            print (lib,str(ct) + ' of ' + str(len(synthstats.keys())) )
            for libcond in tqdm(synthstats [lib].keys ()):
                for synth in synthstats [lib][libcond]:
                    if synth == 'total':
                        continue
                    cycles = synth.split ('|')
                    cyctagct = 1
                    for c in cycles:
                        cyctagct *= tagdict [lib][c]['cycle_tagct']
                    for b in synthstats [lib][libcond][synth]:
                        idx = 0
                        enc = 1
                        for c in cycles:
                            if len (cycles) == 1:
                                cycenc = tagdict [lib][c]['BBs'][b]['tagct']
                                idx += 1
                            if len(cycles) > 1:
                                cycenc = tagdict[lib][c]['BBs'][b[idx]]['tagct']
                                idx += 1
                            enc = enc * cycenc

                        metric_dict = self.Calculate_Synth_Metrics (synthstats[lib][libcond][synth][b]['ct'], synthstats[lib][libcond]['total'], enc, cyctagct  )
                        synthstats [lib][libcond][synth][b]['metric_dict'] = metric_dict
        return synthstats

    def Calculate_Synth_Metrics (self, k_counts, n_selectioncount, tagct, N_totaltags):
        def round_sig(x, sig=2):
            return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)
        p_o = k_counts/n_selectioncount
        p_i= tagct/N_totaltags
        Faver_z_n = math.sqrt(p_i / (1 - p_i)) * ((p_o / p_i) - 1)
        z_alpha = 2.576  # 98% CI, can adjust for different confidence level
        n_prime = n_selectioncount + math.pow(z_alpha, 2)
        conf_p_o = (1 / n_prime) * (k_counts + pow(z_alpha / 2, 2))
        Faver_low_conf = conf_p_o - z_alpha * math.sqrt((conf_p_o / n_prime) * (1 - conf_p_o))
        Faver_high_conf = conf_p_o + z_alpha * math.sqrt((conf_p_o / n_prime) * (1 - conf_p_o))
        Faver_z_n_low = math.sqrt (p_i/(1-p_i))*((Faver_low_conf/p_i)-1)
        Faver_z_n_high = math.sqrt(p_i / (1 - p_i)) * ((Faver_high_conf / p_i) - 1)
        Faver_enr_low = Faver_low_conf/p_i
        Faver_enr = p_o/p_i
        Faver_enr_high = Faver_high_conf/p_i
        poi_lambda = (tagct / N_totaltags) * n_selectioncount
        effect_size = (k_counts - poi_lambda) / math.sqrt(poi_lambda)

        metric_dict = {}
        metric_dict ['effect_size'] = round(effect_size,2)
        metric_dict ['CI_enr'] = [round(Faver_enr_low, 2), round(Faver_enr, 2), round(Faver_enr_high, 2)]

        return metric_dict


    def Add_SynthStatsToAlphaMaFiles(self, filename, lib, libcond,  synthstats):
        def taskfcn(row, lib, libcond, synthlist, synthstats, cols):
            outvals = {}
            for disyn in synthlist:
                b1 = row[disyn[0]]
                b2 = row[disyn[1]]
                synth = disyn[0] + '|' + disyn[1]

                if (b1, b2) in synthstats[lib][libcond][synth]:
                    metric_dict = synthstats[lib][libcond][synth][(b1, b2)]['metric_dict']
                else:
                    metric_dict = None

                if metric_dict is not None:
                    outvals[synth + '_Eff.Sz.'] = metric_dict['effect_size']
                    outvals[synth + 'CI_enr_low'] = metric_dict['CI_enr'][0]
                    outvals[synth + '_enr'] = metric_dict['CI_enr'][1]
                    outvals[synth + 'CI_enr_high'] = metric_dict['CI_enr'][2]
                else:
                    outvals[synth + '_Eff.Sz.'] = 0
                    outvals[synth + 'CI_enr_low'] = 0
                    outvals[synth + '_enr'] = 0
                    outvals[synth + 'CI_enr_high'] = 0

            for cond in synthstats[lib].keys():
                for disyn in synthlist:
                    synth2 = disyn[0] + '|' + disyn[1]
                    if (b1, b2) in synthstats[lib][cond][synth2]:
                        outvals[cond + '_' + synth2 + '_Eff.Sz.'] = \
                            synthstats[lib][cond][synth2][(b1, b2)]['metric_dict']['effect_size']
                    else:
                        outvals[cond + '_' + synth2 + '_Eff.Sz.'] = 0
            colvals = []
            for c in cols:
                colvals.append(outvals[c])
            return colvals

        df = pd.read_csv(filename)
        coldict = {}
        bbcols = []
        for c in df.columns:
            if 'bb' in c:
                coldict[c] = c.replace ('bb', 'BB')
                if '_' not in c:
                    bbcols.append (c.replace ('bb', 'BB'))
        synthlist = []
        for b1 in range(0, len(bbcols) - 1):
            for b2 in range(b1 + 1, len(bbcols)):
                synthlist.append([bbcols[b1], bbcols[b2]])
        df = df.rename (columns=coldict)


        usedask = True
        if usedask == True:
            ProgressBar().register()
            meta_dict = {}
            print ('using dask')
            ddf = dd.from_pandas(df, npartitions=self.NUM_WORKERS)
            cols = []
            for s in synthlist:
                synth2 = s[0] + '|' + s[1]
                cols.extend ([synth2 + '_Eff.Sz.', synth2 + 'CI_enr_low', synth2 + '_enr', synth2 + 'CI_enr_high'])
            for cond in synthstats[lib].keys():
                for s in synthlist:
                    synth2 = s[0] + '|' + s[1]
                    cols.append (cond + '_' + synth2 + '_Eff.Sz.')
            print (cols)
            for i in range (0, len(cols)):
                meta_dict [i] = 'float'
            res = ddf.apply(taskfcn, axis=1, result_type='expand', args = ( lib, libcond,  synthlist, synthstats, cols), meta = meta_dict).compute ()
            res.columns = cols
            df = df.merge(res, left_index=True, right_index=True)
            print (df)
        else:
            for disyn in synthlist:
                synth = disyn[0] + '|' + disyn[1]
                df[synth + '_Eff.Sz.'] = None
                df[synth + 'CI_enr_low'] = None
                df[synth + '_enr'] = None
                df[synth + 'CI_enr_high'] = None

            for cond in synthstats[lib].keys():
                for disyn in synthlist:
                    synth = disyn[0] + '|' + disyn[1]
                    df[cond + '_' + synth + '_Eff.Sz.'] = None
            for idx, r in tqdm(df.iterrows (), total=len(df), position=0, leave=True):
                for disyn in synthlist:
                    b1 = r[disyn[0]]
                    b2 = r[disyn[1]]
                    synth = disyn[0] + '|' + disyn[1]
                    if (b1, b2 ) in synthstats[lib][libcond][synth]:
                        metric_dict = synthstats [lib][libcond][synth][(b1, b2)]['metric_dict']
                    else:
                        metric_dict = None

                    if metric_dict is not None:
                        df.loc[idx, synth + '_Eff.Sz.'] = metric_dict ['effect_size']
                        df.loc[idx, synth + 'CI_enr_low'] = metric_dict['CI_enr'][0]
                        df.loc[idx, synth + '_enr'] = metric_dict['CI_enr'][1]
                        df.loc[idx, synth + 'CI_enr_high'] = metric_dict['CI_enr'][2]
                    else:
                        df.loc[idx, synth + '_Eff.Sz.'] = 0
                        df.loc[idx, synth + 'CI_enr_low'] = 0
                        df.loc[idx, synth + '_enr'] = 0
                        df.loc[idx, synth + 'CI_enr_high'] = 0
                for cond in synthstats[lib].keys ():
                    for disyn in synthlist:
                        synth2 = disyn[0] + '|' + disyn[1]
                        if (b1, b2) in synthstats[lib][libcond][synth2]:
                            df.loc[idx, cond + '_' + synth2 + '_Eff.Sz.'] = synthstats[lib][libcond][synth2][(b1, b2)]['metric_dict']['effect_size']
                        else:
                            df.loc[idx, cond + '_' + synth2 + '_Eff.Sz.'] = 0

        outfilename = filename.replace('proc_struc_data', 'metric_struc_data').replace ('.csv', '.metrics.csv')
        outdirname = os.path.dirname(outfilename)
        if not pathlib.Path(outdirname).exists():
            os.makedirs(outdirname)
        df.to_csv (outfilename)

    def Calculate_AlphaMa_SynthStats (self, synthstats, statspath,  tagspath, liblist, condmap):
        tagdict = {}
        statslist = []
        if synthstats is None:
            synthstats=  {}

        globlist = pathlib.Path(statspath).rglob('AMA*.txt')
        for f in globlist:
            if liblist is not None:
                lib = f.name.split('_')[0]
                if lib in liblist:
                    statslist.append(str(f))
            else:
                statslist.append(str(f))

        for f in tqdm(statslist):
            parsef = f.split('/')
            libcond = parsef[len(parsef) - 1].replace('.txt', '')
            library = libcond.split('_')[0]
            cond = libcond.split('_', 1)[1]
            if library not in synthstats:
                synthstats[library] = {}

            libtagfile = tagspath + '/' + library + '_Tags.txt'
            if cond in condmap:
                libcond = condmap[cond]

            synthstats = self.CalculateAlphaMa_Disynthon_Counts(library, libcond, f, synthstats)
            if not library in tagdict:
                tagdict = self.Get_TagStats(library, libtagfile, tagdict)

        synthstats = self.Add_EnrichmentsToSynthStats(liblist, tagdict, synthstats)

        return synthstats

    def Generate_AlphaMa_SynthStats (self, filepath, statspath,  tagspath, liblist, condmap):
        filelist = []
        synthstats = {}
        synthstats = self.Calculate_AlphaMa_SynthStats (synthstats, statspath,  tagspath, liblist, condmap)
        # tagdict = {}
        # statslist = []
        #
        # synthstats = {}
        #
        # globlist = pathlib.Path(statspath).rglob('AMA*.txt')
        # for f in globlist:
        #     if liblist is not None:
        #         lib = f.name.split('_')[0]
        #         if lib in liblist:
        #             statslist.append(str(f))
        #     else:
        #         statslist.append (str(f))
        #
        # for f in tqdm(statslist):
        #     parsef = f.split('/')
        #     libcond = parsef[len(parsef) - 1].replace ('.txt', '')
        #     library = libcond.split('_')[0]
        #     cond = libcond.split('_',1)[1]
        #     if library not in synthstats:
        #         synthstats [library] = {}
        #
        #     libtagfile = tagspath + '/' + library + '_Tags.txt'
        #     if cond in condmap:
        #         libcond = condmap [cond]
        #
        #     synthstats = self.CalculateAlphaMa_Disynthon_Counts(library, libcond, f, synthstats)
        #     if not library in tagdict:
        #         tagdict = self.Get_TagStats(library, libtagfile, tagdict)
        #
        #
        # synthstats = self.Add_EnrichmentsToSynthStats(tagdict, synthstats)

        globlist = pathlib.Path(filepath).rglob('AMA*.csv')
        for f in globlist:
            if liblist is not None:
                lib = f.name.split('_')[0]
                if lib in liblist:
                    filelist.append(str(f))
            else:
                filelist.append(str(f))

        for f in tqdm (filelist):
            parsef = f.split('/')
            filename = parsef [len(parsef) - 1]
            splitf = filename.split ('_')
            libcond = filename [0:filename.find('_compare')]
            lib = splitf[0]

            if liblist is not None:
                if lib in liblist:
                    self.Add_SynthStatsToAlphaMaFiles (f, lib, libcond, synthstats)
            else:
                self.Add_SynthStatsToAlphaMaFiles (f, synthstats, libcond, synthstats)

    def Reformat_AlphaMaTagsToCephalogix (self, library, tagfile, strucfile, outfilepath, rxns):
        df= pd.read_csv (tagfile, header = None, sep='\t')
        sdf =  pd.r (strucfile, header = None, sep='\t')
        sdf = sdf.rename (columns = {0:'SMILES', 1:'Cycle_BBID'})
        cols = {}
        for c in df.columns:
            if c == 0:
                cols [c] = 'Cycle'
            elif c == 1:
                cols [c] = 'BBID'
            else:
                cols [c] = 'Tag' + str (c-2)
        df = df.rename(columns=cols)
        Cycles = df['Cycle'].unique()
        Tags = [c for c in df.columns if 'Tag' in c]

        for cyc in Cycles:
            outdf = pd.DataFrame(None, columns=['id', 'codon', 'smiles', 'reaction_name'])
            cdf = df[df['Cycle']==cyc]
            for idx, row in cdf.iterrows ():
                for t in Tags:
                    cyc_bbid = str(cyc) + '_' + str(row['BBID'])
                    smiles = sdf[sdf['Cycle_BBID'] == cyc_bbid].iloc[0]['SMILES']
                    outdf.loc[len(outdf.index)] = [row['BBID'], row[t],smiles,rxns [cyc]]
            outfile = outfilepath + '/' + library + '/bb' + str(cyc) + '.csv'
            if not pathlib.Path(outfilepath + '/' + library).exists():
                os.makedirs(outfilepath + '/' + library)
            outdf.to_csv (outfile, index = False)

    def Remove_arch_from_fastq(self, infile, outfile, seq1, seq2, dist, vardist=2, keep=False):
        def read_block(f, line_ct):
            llist = []
            for i in range (0,line_ct):
                data = f.readline ()
                if data:
                    llist.append (data)
                else:
                    return None
            return llist
        of = open(outfile,'w')

        with open(infile) as f:
            nextlines = read_block(f, 4)
            rct = 0
            ct = 0
            while nextlines:
                nextlines = read_block(f, 4)
                if nextlines is None:
                    break
                start = nextlines [1].find (seq1)
                end = nextlines [1].find (seq2)
                if(start != -1 and end != -1):
                    if (end-start-vardist <= dist and end-start+vardist >= dist  ):

                        if keep== True:
                            of.write(''.join (nextlines))
                        else:
                            rct += 1
                    else:
                        if keep == False:
                            of.write(''.join(nextlines))
                        else:
                            rct+=1
                else:
                    if keep == False:
                        of.write(''.join (nextlines))
                    else:
                        rct +=1
                ct += 1
                if ct % 1000 == 0:
                    print (rct, ct, end='\r')
        of.close ()
    def Add_SynthStatsToYoDELFiles(self, filepath):
        def taskfcn(row, synthlist, synthstats, cols):
            outvals = {}
            lib = row['library.name']
            for disyn in synthlist:
                b1 = str(row[disyn[0]])
                b2 = str(row[disyn[1]])
                currbbs = b1 + '|' + b2
                synth = disyn[0] + '|' + disyn[1]

            #     if currbbs in synthstats[lib][synth]:
            #         metric_dict = synthstats[lib][synth][currbbs]['metric_dict']
            #     else:
            #         metric_dict = None
            #
            #     if metric_dict is not None:
            #         outvals[synth + '_Eff.Sz.'] = metric_dict['effect_size']
            #         outvals[synth + 'CI_enr_low'] = metric_dict['CI_enr'][0]
            #         outvals[synth + '_enr'] = metric_dict['CI_enr'][1]
            #         outvals[synth + 'CI_enr_high'] = metric_dict['CI_enr'][2]
            #     else:
            #         outvals[synth + '_Eff.Sz.'] = 0
            #         outvals[synth + 'CI_enr_low'] = 0
            #         outvals[synth + '_enr'] = 0
            #         outvals[synth + 'CI_enr_high'] = 0
            #
            # for cond in synthstats[lib].keys():
            #     for disyn in synthlist:
            #         synth2 = disyn[0] + '|' + disyn[1]
            #         if currbbs  in synthstats[lib][cond][synth2]:
            #             outvals[cond + '_' + synth2 + '_Eff.Sz.'] = \
            #                 synthstats[lib][cond][synth2][currbbs]['metric_dict']['effect_size']
            #         else:
            #             outvals[cond + '_' + synth2 + '_Eff.Sz.'] = 0
            colvals = []
            for c in cols:
                colvals.append(outvals[c])
            return colvals

        filename = filepath + 'merged_molecule_enrichments.csv'
        df = pd.read_csv(filename)
        df = df.head ()
        coldict = {}
        bbcols = []

        for c in df.columns:
            if c.startswith ('bb'):
                coldict[c] = c.split ('.',1)[0].replace ('bb', 'BB')
                if coldict[c] not in bbcols:
                    bbcols.append (coldict[c])

        synthlist = []
        for b1 in range(0, len(bbcols) - 1):
            for b2 in range(b1 + 1, len(bbcols)):
                synthlist.append([bbcols[b1], bbcols[b2]])
        df = df.rename (columns=coldict)

        synthfilename = filepath + 'merged_nsynthon_enrichments.csv'
        synthdf = pd.read_csv(synthfilename)

        ProgressBar().register()
        meta_dict = {}
        ddf = dd.from_pandas(df.head(), npartitions=self.NUM_WORKERS)
        cols = []

        for s in synthlist:
            synth2 = s[0] + '|' + s[1]
            cols.extend ([synth2 + '_Eff.Sz.', synth2 + '_CI_enr_low', synth2 + '_enr', synth2 + '_CI_enr_high'])

        synthstats = {}
        for lib in synthdf['library.name'].unique ():
            if  lib not in synthstats:
                synthstats [lib] = {}
            for s in synthlist:
                synthstats[lib][s[0] + '|' + s[1]] = {}

        # LINE_UP = '\033[1A'
        # LINE_CLEAR = '\x1b[2K'
        #
        # print ('line 1')
        # print('line 2')
        # print(LINE_UP, end=LINE_CLEAR)

        # for idx, row in synthdf.head ().iterrows ():
        #     lib = row ['library.name']
        #     synth2 = row['nsynth']
        #     bbs = synth2.split ('.')
        #     if len(bbs) == 2:
        #         synth2 = synth2.upper ().replace ('.','|')
        #         bbstr = str(row[bbs[0]+'.id']) + '|' + str(row[bbs[1]+'.id'])
        #         synthstats[lib][synth2] [bbstr] = row
        #         print (idx, end = '\r')
        #
        #
        # for i in range (0, len(cols)):
        #     meta_dict [i] = 'float'
        # res = ddf.apply(taskfcn, axis=1, result_type='expand', args = (synthlist, synthstats, cols), meta = meta_dict).compute ()
        # res.columns = cols
        # df = df.merge(res, left_index=True, right_index=True)
        # print (df)
        #
        #
        # outfilename = filename.replace('proc_struc_data', 'metric_struc_data').replace ('.csv', '.metrics.csv')
        # outdirname = os.path.dirname(outfilename)
        # if not pathlib.Path(outdirname).exists():
        #     os.makedirs(outdirname)
        # df.to_csv (outfilename)

if __name__ == "__main__":
    dlpt = DELProcessingTools ()
    dlpt.Add_SynthStatsToYoDELFiles ('/Users/eric/Consulting/Citadel/Partners/MOMA/AlphaMa_Data/2022_06_27_Reprocessed/PMOM/')