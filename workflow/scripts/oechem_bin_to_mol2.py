from openeye import oechem
import numpy as np

# Load the reference molecule into memory, to use it later. Works
ifs = oechem.oemolistream()
ofs = oechem.oemolostream()
ifs.open(snakemake.input.in_file)
mollist = []
for mol in ifs.GetOEMols():
    if ofs.open(snakemake.output.out_file):
        oechem.OEWriteMolecule(ofs, mol)
        # Retrieve relative energies of conformers and write to txt file
        if snakemake.output.conf_energy:
            energy = []
            for confs in mol.GetConfs():
                if oechem.OEHasSDData(confs, "Relative Energy"):
                    energy.append(
                        float(oechem.OEGetSDData(confs, "Relative Energy"))
                    )
            energy = np.array(energy)
            np.savetxt(snakemake.output.conf_energy, energy)
    else:
        oechem.OEThrow.Fatal("Unable to create 'output file'")
    # print(mol.GetTitle())
    # mollist.append(mol)
    # for confs in mol.GetConfs():
    # print(confs)
# for mol in ifs.GetOEGraphMols():
#    mollist.append(oechem.OEGraphMol(mol))

# ref_conf = mollist[0]
# ref_conf_rdkit = Chem.MolFromPDBFile(ref_path)
