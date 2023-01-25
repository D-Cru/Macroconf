import os

compounds = ["regular"]  # , "cis_only", "cis_trans"]
methods = ["GaMD", "aMD", "cMD"]

files = []
noes = []
dihedrals = []
ids = []
fdir = fdir = os.path.dirname(os.path.abspath(__file__))
for c in compounds:
    for m in methods:
        files.append(
            (
                os.path.join(fdir, "data", c, m, "traj.netcdf"),
                os.path.join(fdir, "data", c, m, "top_sol.prmtop"),
            )
        )
        ids.append(f"{c}_{m}")

# @pytest.fixture(params=files, ids=ids)
# def mdtraj_traj(request):
#     trajectory, topology = request.param
#     mdtraj_traj = md.load(trajectory, top=topology, stride=1)
#     return mdtraj_traj

# @pytest.fixture
# def mdtraj_traj_protein(mdtraj_traj):
#     mdtraj_traj = mdtraj_traj.atom_slice(
#       mdtraj_traj.topology.select('protein'))
#     return mdtraj_traj

# @pytest.fixture
# def mdtraj_frame(mdtraj_traj_protein):
#     mdtraj_traj_protein = mdtraj_traj_protein[0]
#     return mdtraj_traj_protein


# def test_trajs(mdtraj_traj):
#     assert mdtraj_traj.n_frames == 76000

# def test_traj_prot(mdtraj_traj_protein):
#     assert mdtraj_traj_protein.n_frames == 76000
