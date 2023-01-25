from rdkit import Chem

# sys.path.append(os.getcwd())
# sys.path.append(f"{os.getcwd()}/libs")
# from PepLibGen.StructGen import StructGen as sg
# import json
# import numpy as np
# import pandas as pd
# from src.utils import dotdict

mol_ref = Chem.MolFromMol2File(
    snakemake.input.ref_mol,
    removeHs=False,
)
# mol_ref.RemoveAllConformers()

# mol = Chem.MolFromSmiles(compound.smile)
# mol = Chem.AddHs(mol)
# Chem.Draw.MolToFile(mol, snakemake.output.mol_image)

with open(snakemake.output.smiles, "w") as f:
    f.write(Chem.MolToSmiles(mol_ref))
