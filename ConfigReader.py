import json
import sys

f = open('./CIxTools.init.json')
jdict = json.loads (f.read())
f.close ()
sys.path.append(jdict ['cixtoolspath'])