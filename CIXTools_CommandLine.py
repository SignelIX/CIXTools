import argparse
import ChemUtilities
#region Run Functions

def SD_to_SMILES (args):
    infile = args ['sd2smiles']
    outfile = args ['outpath']
    idcol = args['id_column']
    ChemUtilities.SDFtoFile (infile, fix = True, outfile= outfile, idcol=idcol)

def SaltStripFile(args):
    infile =  args ['saltstrip']
    outfile = args ['outpath']
    smilescol = args ['smiles_column']
    ChemUtilities.SaltStripFile(infile, smilescol, outfile)

#endregion Run Functions


def Setup_Args (parser):
    #Actions
    parser.add_argument("-sd2smiles", "--sd2smiles", help="convert sdfile to smiles csv")
    parser.add_argument("-saltstrip", "--saltstrip", help="convert sdfile to smiles csv")


    parser.add_argument("-o", "--outpath", help="output path")
    parser.add_argument("-id", "--id_column", help="column containing molecular id, e.g. vendor id")
    parser.add_argument("-smiles", "--smiles_column", help="column containing smiles strings")

    args = vars(parser.parse_args())
    return args

cmddict = {
           'sd2smiles': SD_to_SMILES,
           'saltstrip': SaltStripFile,
           }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = Setup_Args(parser)

    for k in cmddict.keys ():
        if args [k] is not None:
            cmddict [k] (args)
            break
