# Create tleap input files automatically
import json

with open(snakemake.input.parm) as f:
    compound_info = json.load(f)
    # compound_info = compound_info["compound"]
# print(compound_info)
# compound_info = compound_info["compound"]
seq_length = compound_info["seq_length"]

# Adapt forcefield depending on whether peptide is HT cyclic
if compound_info["bonds"][0] == "HT":
    bond = f"bond ramp.1.N ramp.{seq_length}.C"
    forcefield = snakemake.input.ff
else:  # Not HT cyclic peptide... Has 1, N sidechains..
    # This should work for the example provided but not generally.
    # need to adapt
    # the bond syntax in json file again. Have bond[0]=HT if HT
    # have bond[1]= SS
    bond = ""  # f"bond ramp.1.SG ramp.{seq_length}.SG"
    forcefield = "leaprc.protein.ff14SB"

# if len(compound_info["bonds"]) == 2:
#     ss_bond_info = compound_info["bonds"][1]
#     ss_positions = []
#     for i,c in enumerate(ss_bond_info[2:]):
#         if c == "X":
#             continue
#         elif c == "C":
#             ss_positions.append(i)

#     bond = f"""
#     {bond}
#     bond ramp.{ss_positions[0]+1}.SG ramp.{ss_positions[1]+1}.SG
#     """


# Solvent settings:
print(snakemake.config)
if snakemake.wildcards.solvent == "H2O":
    solvent_ff = "source leaprc.water.tip3p"
    box = "TIP3PBOX"
    box_source = ""
    box_size = "12.0"
    charge = """
    # neutralize system. this will create 1 Warning...
    addIons ramp Na+ 0
    addIons ramp Cl- 0
    """
elif snakemake.wildcards.solvent == "DMSO":
    solvent_ff = f"loadAmberParams {snakemake.config['DMSO_params']}"
    box = "d"
    box_size = "12.0"
    box_source = "loadoff libs/md_solvents/dmso/dmsobox.off"
    charge = ""

elif snakemake.wildcards.solvent == "DMSO_b":
    solvent_ff = f"loadAmberParams {snakemake.config['DMSO_params']}"
    box = "d"
    box_size = "20.0"
    box_source = "loadoff libs/md_solvents/dmso/dmsobox.off"
    charge = ""

elif snakemake.wildcards.solvent == "DMSO_GAFF_BCC":
    solvent_ff = "source leaprc.gaff"
    box = "dmsobox"
    box_source = "loadoff libs/md_solvents/dmso_gaff_bcc/dmsobox.lib"
    box_size = "12.0"
    charge = ""

elif snakemake.wildcards.solvent == "DMSO_GAFF_BCC_b":
    solvent_ff = "source leaprc.gaff"
    box = "dmsobox"
    box_source = "loadoff libs/md_solvents/dmso_gaff_bcc/dmsobox.lib"
    box_size = "20.0"
    charge = ""

elif snakemake.wildcards.solvent == "DMSO_GAFF_RESP":
    solvent_ff = "source leaprc.gaff"
    box = "dmsobox"
    box_source = "loadoff libs/md_solvents/dmso_gaff_resp/dmsobox.lib"
    box_size = "12.0"
    charge = ""

elif snakemake.wildcards.solvent == "DMSO_GAFF_RESP_b":
    solvent_ff = "source leaprc.gaff"
    box = "dmsobox"
    box_source = "loadoff libs/md_solvents/dmso_gaff_resp/dmsobox.lib"
    box_size = "20.0"
    charge = ""

elif snakemake.wildcards.solvent == "Chloroform":
    solvent_ff = f"loadAmberParams {snakemake.config['Chloroform_params']}"
    box = "CHCL3BOX"
    box_source = "loadoff chcl3box.off"
    box_size = "12.0"
    charge = ""

elif snakemake.wildcards.solvent == "Chloroform_b":
    solvent_ff = f"loadAmberParams {snakemake.config['Chloroform_params']}"
    box = "CHCL3BOX"
    box_source = "loadoff chcl3box.off"
    box_size = "25.0"
    charge = ""
else:
    print("error. solvent not available.")
out = f"""#tleap.in
logFile {snakemake.log}
source {forcefield}
{solvent_ff}
ramp=loadpdb {snakemake.input.pdb}

{bond}

saveAmberParm ramp {snakemake.params.path}mc_gas.prmtop {snakemake.params.path}mc_gas.inpcrd

savepdb ramp {snakemake.params.path}mc_gas.pdb

saveMol2 ramp {snakemake.params.path}mc_gas.mol2 0

charge ramp

{charge}

check ramp

{box_source}
solvateOct ramp {box} {box_size}

# add solvent
# addIonsRand ramp Na+ 19 Cl- 19

#check ramp

saveamberparm ramp {snakemake.params.path}mc_sol.prmtop {snakemake.params.path}mc_sol.inpcrd

savepdb ramp {snakemake.params.path}mc_sol.pdb

quit
"""

with open(snakemake.output.tleap, "w") as f:
    f.write(out)
