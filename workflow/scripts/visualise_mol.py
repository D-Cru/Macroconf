import sys
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdMolDraw2D

IPythonConsole.ipython_useSVG = (
    True  # < set this to False if you want PNGs instead of SVGs
)
# dir(IPythonConsole)
IPythonConsole.molSize = (900, 300)  # (450, 150)
IPythonConsole.drawOptions.addStereoAnnotation = True
IPythonConsole.drawOptions.annotationFontScale = 1.5


def Mol_stereo(mol_amber, molSize=(900, 300), fontSize=0.8, LineWidth=1):
    """
    Draw a molecule with RDKit and return the SVG string
    """

    mol_ref = Chem.MolFromMol2File(
        mol_amber,
        removeHs=False,
    )
    mol_ref.RemoveAllConformers()

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mol_ref)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg


if __name__ == "__main__":
    Mol_stereo(sys.argv[1])
