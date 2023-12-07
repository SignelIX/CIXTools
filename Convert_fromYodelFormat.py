import toml
import json
import glob
import pathlib
import os
import shutil
def Read_YodelScheme (schemename, filename, outfile):
    print (schemename)
    f = open (filename, 'r')
    yodelscheme = toml.load (f)
    f.close ()
    schemedict = {}
    steps = []
    schemedict ['Named_Reactions'] = {}
    isfirst = True
    for w in yodelscheme ['reactions'].keys ():
        step = yodelscheme ['reactions'][w]
        dellist = []
        for s in step.keys ():
            if s not in  ['reactants', 'reactions']:
                if s in schemedict ['Named_Reactions']:
                    rct = 1
                    rname = s + '_' + str(rct)
                    while rname in schemedict ['Named_Reactions']:
                        rct +=1
                        rname = s + '_' + str(rct)
                    schemedict['Named_Reactions'][rname] = step[s]
                    step ['reactions'] = [rname for sr in step ['reactions'] if sr == s]
                else:
                    schemedict['Named_Reactions'] [s] = step [s]
                dellist.append (s)

        for d in dellist:
            del step [d]
        step ['Reactants'] = fix_reactants (step ['reactants'], isfirst)
        del step ['reactants']
        step['Rxns'] = {'default':step['reactions']}
        del step['reactions']
        steps.append (step)
        isfirst = False

    schemedict [schemename] = {'steps': steps}
    f = open (outfile, 'w')
    f.write (json.dumps(schemedict, indent = 4))
    f.close ()

def fix_reactants (rxtants, isfirst):
    fixstep = []
    for r in rxtants:
        if r in ['bb1', 'bb2', 'bb3']:
            r = 'r' + str (int (r.replace ('bb', '')) - 1)
            fixstep.append (r)
        else:
            if isfirst:
                fixstep.append ('CN')
            else:
                fixstep.append ('p')
    if len (fixstep) < 2:
        fixstep.append ('None')
    return fixstep

def convert_reaction_store (infile, outfile):
    f = open(infile, 'r')
    yodelrxns = json.load(f)
    f.close()
    rxndict = {}
    for entry in yodelrxns:
        rxndict [entry ['name']] = entry['smarts']
    f = open(outfile, 'w')
    f.write (json.dumps (rxndict, indent=4))
    f.close()

def convert_BBfilesNames (inpath):
    dirlist = glob.glob(inpath + '*/')
    for d in dirlist :
        libraryname = str (pathlib.Path(d).name)
        libraryfldr = inpath + libraryname
        os.makedirs(libraryfldr + '/BBLists/', exist_ok=True)
        flist = glob.glob (libraryfldr + '/bb*.csv')
        print (flist)
        for f in flist:
            print (f)
            shutil.copy(str(f), str(f).replace ('bb','BBLists/' + libraryname + '.BB' ))

if __name__ == "__main__":
    # inpath = '/Volumes/GoogleDrive/Shared drives/Citadel Projects/Vendors/Hitgen/opendel3x/libraries/'
    # convert_BBfilesNames(inpath)

    for i in range (3, 31):
        library = 'OpenDEL00'  + str(i).zfill (2)
        schemename = library
        schemefile ='/Volumes/GoogleDrive/Shared drives/Citadel Projects/Vendors/Hitgen/opendel3x/libraries/' + library + '/reaction.scheme'
        outfile = '/Users/eric/Temp/Hitgen_3.x_json/' + library + '.json'
        Read_YodelScheme(schemename, schemefile, outfile)
    # filename = '/Users/eric/OneDrive/CIx/reactions_store.json'
    # ofilename = '/Users/eric/OneDrive/CIx/yodel_rxns.json'
    # convert_reaction_store(filename, ofilename)