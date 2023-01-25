from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import rdBase

print(rdBase.rdkitVersion)

mol_ref = Chem.MolFromMol2File(snakemake.input.mol, removeHs=False)
mol_ref.RemoveAllConformers()

Draw.MolToFile(mol_ref, snakemake.output.png)
