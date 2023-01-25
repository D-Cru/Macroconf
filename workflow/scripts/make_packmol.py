# Create tleap input files automatically
import json
import pandas as pd
import sys

out = f"""#packmol input file
# All atoms from diferent molecules will be at least 2.0 Angstroms apart
tolerance 2.0

# The type of the files will be pdb
filetype pdb

# The name of the output file
output {snakemake.params.system}

# put the COM of the solute at the center of the box
structure {snakemake.input.solute}
  number 1
  fixed 15. 15. 15. 0. 0. 0.
  centerofmass
end structure

# # add first type of solvent molecules
# structure {snakemake.input.solvent1}
#   number 150
#   inside cube 0. 0. 0. 30
# end structure

# add second type of solvent molecules
structure {snakemake.input.solvent2}
  number 150
  inside cube 0. 0. 0. 30
end structure
"""

with open(snakemake.output.packmol, "w") as f:
    f.write(out)
