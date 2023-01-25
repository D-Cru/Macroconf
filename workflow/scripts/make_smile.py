from rdkit import Chem
import os
import sys
import json

sys.path.append(os.getcwd())
sys.path.append(f"{os.getcwd()}/libs")
from PepLibGen.StructGen import StructGen as sg  # noqa E402
from src.utils import dotdict  # noqa E402


with open(snakemake.input.json, "r") as f:
    compound = json.load(f)

compound = dotdict(compound)
# dict_keys(['index', 'seq_length', 'bonds', 'sequence', 'sequence_1',
# 'natural_cyclic_peptide', 'non_natural_cyclic_peptide', 'smile'])

# Do bond manipulation here if necessary
if len(compound.bonds) == 1:
    bonds = compound.bonds[0]

    compound["smile"] = sg.constrained_peptide_smiles(
        compound.sequence_1, bonds
    )[2]
elif len(compound.bonds) == 2:
    compound["smile"] = compound.custom_smile

mol = Chem.MolFromSmiles(compound.smile)
mol = Chem.AddHs(mol)
Chem.Draw.MolToFile(mol, snakemake.output.mol_image)

with open(snakemake.output.json_all, "w") as f:
    json.dump(compound, f)

with open(snakemake.output.smile, "w") as f:
    f.write(compound.smile)
