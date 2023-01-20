<h2>Background</h2>
<p>
Authors: Eric Sigel
<p>
Contributors: Antoine Dumas

CIxTools is an accumulated set of Cheminformatics tools built on rdkit that I have been building for my consulting work
The tools are developed for DEL (DNA Encoded Library) and Virtual library analysis and enumeration
Many functions are accessible from a lightweight frontend built with streamlit, increasingly functionality is accessible
through CLI options.

<h2>Setup</h2>
<h3>Conda env</h3>
conda env create environment.yml
conda activate CixTools
pip install -r requirements.txt

<h3>Folder structure</h3>
*Where file names can be specified, the names are suggestions not requirements. 
Recommended installation, create a folder called CIx.
At top level of CIx create a folder called Enumerations and a folder called CIxTools
Additionally create a file called CIxTools.init.json, which can be an empty dictionary ('{}')
This file is used to store default values.
Optionally, for some CLI functionality, a file called inparams.yml can be created at the CIx level

<h3>For Enumerations:</h3>
Under the Enumerations folder create a .json file called RxnSchemes.json
For each reaction scheme that is added to RxnSchemes.json, create a subfolder with the scheme name
Under the scheme folder create a subfolder called BBLists
In each subfolder there should be one file per BB cycle named <scheme>.BB<X>.csv
These files should have at a minimum columns BB_ID, SMILES
Additional libraries using the same scheme can be added using a specstr (specification string)
in the scheme folder create a subfolder named <specstr>
under this folder create a subfolder called BBLists
BB files in this set up should be named <scheme>.<specstr>.BB<X>.csv
This specstr can be entered into the UI to refer to the specific BB sets
<p></p>

![img1.png](Readme_Images/img_1.png)


<p></p>
<h2>Running the streamlit UI</h2>
<p>
The streamlit UI is launched with the following command:
streamlit run CIxToolsMainUI.py
This will provide a URL for access in std out.
The UI is set up to provide multiple functionalities that can be accessed through buttons on the side pane.
To access Enumeration, click the top button under controls

![img.png](Readme_Images/img.png)

<h2>Enumeration</h2>
Groups that I work with increasingly are interested in enumerating chemical space with building blocks in
large scale combinatorial libraries. Building SMARTS/SMIRKS reactions is far from my favorite activity. 
The ability to visually inspect an enumeration, examine random examples, and check specific failing examples are 
intended to be simplified using the (simple/clunky) UI. Features include product structure as well as structures 
of building blocks and intermediates. A grid that shows each step improves ability to debug SMARTS/SMIRKS failures.
A further table of random structures, shows if there are any failures featuring clickable rows that place the relevant BBs
into the relevant structure override textbox. The UI provides a mechanism for editing the current reaction
scheme and testing/saving to the active RxnSchemes.json file for further use.

Example json scheme:
</p>

    "cix1": {
        "steps": [
            {
                "Reactants": [
                    "CN",
                    "r0"
                ],
                "Rxns": {
                    "default": "[NX3:1].[C:2](=[O:3])[OX2;H1]>>[N:1][C:2](=[O:3])"
                }
            },
            {
                "Reactants": [
                    "p",
                    "r1"
                ],
                "Rxns": {
                    "default": "C(C)(C)(C)OC(=O)[NX3;H:1].[C:2](=[O:3])[OX2;H1]>>[N:1][C:2](=[O:3])"
                }
            },
            {
                "Reactants": [
                    "p",
                    "r2"
                ],
                "Rxns": {
                    "default": "O=C(OCC1c2c(c3c1cccc3)cccc2)[NX3:1].[C:2](=[O:3])[OX2;H1]>>[N:1][C:2](=[O:3])"
                }
            }
        ]
    }

<p></p>
<h3>Basic Elements</h3>
each scheme can be separated by comma
<p>

    "schemename": {
        "steps":[
             {
                "Reactants": [
                    "SMILES or r[N] or p",   # N is a 0 based number ordered starting with 1st reactant etc.  r0,r1, r2
                                               # a SMILES string is used explicitly
                                               # p indicates to use the product from the previous step
                    "r[N] or None"             # r[N] indicates to use the Nth reactant, None indicates that there is only 1 reactant
                ],
                "Rxns": {
                    "default or SMARTS": ["rxn SMIRKS 1", "rxn SMIRKS 2", etc] # reactions are done in order where 
                }
            }
        ]   
    }
<h3>Simple Enumeration</h3>
In the Rxn Scheme File  text box paste path to RxnSchemes.json file (this can be an alternative name to allow e.g. dev versions)
Re-Initialize button will reload the schemes
Scheme path should point to folder (e.g. Enumerations) where scheme folders are located (e.g. containing BBLists)
requires trailing slash
Entering the correct information should cause the Schemes dropdown to appear and populate
![img_4.png](Readme_Images/img_4.png)

<h3>Export</h3>
Export provides a mechanism to export a random sample of a specified size (positive number)
Or the full library (-1) in the count field
Remove Duplicates check box will identify rows with the same structure and take the first row only
Add Intermediates check box will add reaction intermediate information to the output including the reaction and product (where possible) of each step

![img_2.png](Readme_Images/img_2.png)


--Incomplete ES 1/19/23--

<h2>generation of install files</h2>
pip list --format=freeze > requirements.txt

