# Compute the aMD parameters from the cMD simulation

# in Pierce:2012 need average dihedral and total potential energy of
# cMD simulation

import numpy as np
import mdtraj as md

# Compute number of atoms and residues via mdtraj

t = md.load(snakemake.input.traj, top=snakemake.input.top)
n_atoms = t.n_atoms
t = t.restrict_atoms(t.topology.select("protein"))
n_residues = t.n_residues  # solute residues

time = np.loadtxt(snakemake.input.epot, comments=["#", "@"])[:, 0]
Epot = np.loadtxt(snakemake.input.epot, comments=["#", "@"])[:, 1]
Edih = np.loadtxt(snakemake.input.edih, comments=["#", "@"])[:, 1]

Epot_avg = np.average(Epot)
Edih_avg = np.average(Edih)

e_dihed = Edih_avg + (4 * n_residues)
alpha_dihed = (1 / 5) * (4 * n_residues)
e_tot = Epot_avg + (0.16 * n_atoms)
alpha_tot = 0.16 * n_atoms

namespace = {
    "e_dihed": e_dihed,
    "alpha_dihed": alpha_dihed,
    "e_tot": e_tot,
    "alpha_tot": alpha_tot,
}

sample_params = snakemake.params.sample
namespace.update(sample_params)

if int(namespace["simtime"]) > 1000:
    namespace["ntpr"] = namespace["simtime"]
    namespace["ntwx"] = namespace["simtime"]
    namespace["ntwr"] = namespace["simtime"]
else:
    namespace["ntpr"] = 1000
    namespace["ntwx"] = 1000
    namespace["ntwr"] = 1000

with open(snakemake.input.template) as f:
    params = f.read()

params_formatted = params.format(**namespace)

with open(snakemake.output.param, "w") as f:
    f.write(params_formatted)
