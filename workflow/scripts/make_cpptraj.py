# Create cpptraj files

out = f"""parm {snakemake.input.parm}
trajin {snakemake.input.traj} lastframe
strip :WAT
strip :Na+,Cl-
trajout {snakemake.params.pdb} pdb
run

"""

with open(snakemake.output.cpptraj_file, "w") as f:
    f.write(out)
