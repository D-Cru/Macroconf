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

# Solvent settings:
# print(snakemake.config)
# solvent_ff = "source leaprc.water.tip3p"
# box = "TIP3PBOX"
# box_source = ""
# box_size = "12.0"
charge = """
# neutralize system. this will create 1 Warning...
addIons SYS Na+ 0
addIons SYS Cl- 0
"""
# elif snakemake.wildcards.solvent == 'DMSO':
#     solvent_ff = f"loadAmberParams {snakemake.config['DMSO_params']}"
#     box = "d"
#     box_size = "12.0"
#     box_source = "loadoff libs/md_solvents/dmso/dmsobox.off"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'DMSO_b':
#     solvent_ff = f"loadAmberParams {snakemake.config['DMSO_params']}"
#     box = "d"
#     box_size = "20.0"
#     box_source = "loadoff libs/md_solvents/dmso/dmsobox.off"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'DMSO_GAFF_BCC':
#     solvent_ff = "source leaprc.gaff"
#     box = "dmsobox"
#     box_source = "loadoff libs/md_solvents/dmso_gaff_bcc/dmsobox.lib"
#     box_size = "12.0"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'DMSO_GAFF_BCC_b':
#     solvent_ff = "source leaprc.gaff"
#     box = "dmsobox"
#     box_source = "loadoff libs/md_solvents/dmso_gaff_bcc/dmsobox.lib"
#     box_size = "20.0"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'DMSO_GAFF_RESP':
#     solvent_ff = "source leaprc.gaff"
#     box = "dmsobox"
#     box_source = "loadoff libs/md_solvents/dmso_gaff_resp/dmsobox.lib"
#     box_size = "12.0"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'DMSO_GAFF_RESP_b':
#     solvent_ff = "source leaprc.gaff"
#     box = "dmsobox"
#     box_source = "loadoff libs/md_solvents/dmso_gaff_resp/dmsobox.lib"
#     box_size = "20.0"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'Chloroform':
#     solvent_ff = f"loadAmberParams {snakemake.config['Chloroform_params']}"
#     box = "CHCL3BOX"
#     box_source = "loadoff chcl3box.off"
#     box_size = "12.0"
#     charge = ""
#
# elif snakemake.wildcards.solvent == 'Chloroform_b':
#     solvent_ff = f"loadAmberParams {snakemake.config['Chloroform_params']}"
#     box = "CHCL3BOX"
#     box_source = "loadoff chcl3box.off"
#     box_size = "25.0"
#     charge = ""
# else:
#     print('error. solvent not available.')


out = f"""#tleap.in
logFile {snakemake.log}
source {forcefield}
source leaprc.water.tip3p
source leaprc.gaff
ramp=loadpdb {snakemake.input.pdb}

SV2 = loadmol2 {snakemake.input.dmso}
loadamberparams {snakemake.input.dmso_mod}

{bond}

saveAmberParm ramp {snakemake.params.path}mc_gas.prmtop {snakemake.params.path}mc_gas.inpcrd

savepdb ramp {snakemake.params.path}mc_gas.pdb

saveMol2 ramp {snakemake.params.path}mc_gas.mol2 0



check ramp

# load system
SYS = loadpdb {snakemake.input.system}

# set box dimensions
set SYS box {{30,30,30}}

solvatebox ramp TIP3PBOX 2

charge SYS

{charge}

saveamberparm SYS {snakemake.params.path}mc_sol.prmtop {snakemake.params.path}mc_sol.inpcrd

savepdb SYS {snakemake.params.path}mc_sol.pdb

quit
"""

with open(snakemake.output.tleap, "w") as f:
    f.write(out)
