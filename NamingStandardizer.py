import json
initjsonfile = '../CIxTools.init.json'
def Pull_ColStandardizer (colname):
    f = open(initjsonfile)
    initjson = json.load(f)
    f.close()
    if 'Column_Standardizer' in initjson:
        if colname in initjson ['Column_Standardizer']:
            return initjson ['Column_Standardizer'][colname]


def Standardize_Colnames ( df, colname):
    colnames = Pull_ColStandardizer(colname)
    changecoldict = {}
    if colnames is not None and len(colnames) > 0:
        for smc in colnames:
            if smc in df.columns:
                changecoldict[smc] = colname
    if len (changecoldict) > 0:
        return df.rename (columns=changecoldict)
    return df