from IPython.display import SVG
from rdkit.Chem import rdDepictor as rdd
from rdkit.Chem.Draw import rdMolDraw2D as draw2d
from rdkit import Chem

def get_molecule():
    # load a mol from a SMILES string
    mol = Chem.MolFromSmiles('CC(N)C(=O)O')

    # load a large molecule database (size=440055 mols)
    #suppl = Chem.SDMolSupplier('database/chembl_30.sdf')
    return mol

def draw2D(m):
    molSize = (450, 150)                    # draw area
    mc = Chem.Mol(m.ToBinary())
    if not mc.GetNumConformers():
        rdd.Compute2DCoords(mc)             # compute 2D coordinates of atoms
    drawer = draw2d.MolDraw2DSVG(\
        molSize[0],molSize[1])              # initialize drawer with size
    drawer.DrawMolecule(mc)                 # draw the molecule
    drawer.FinishDrawing()
    #svg = drawer.GetDrawingText()           # get the SVG string
    drawer.WriteDrawingText('test_00.png')
    #display(SVG(svg.replace('svg:','')))    # fix and display in Jupyter notebook


mol = get_molecule()
draw2D(mol)
