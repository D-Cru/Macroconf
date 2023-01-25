"""
Test pca functions.

"""

import pytest
import mdtraj as md
import os
import src.pca
import src.dihedrals
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


@pytest.fixture
def example_mdtraj_frame():
    fdir = os.path.dirname(os.path.abspath(__file__))
    mol = os.path.join(fdir, "data", "regular", "22_solvated.prmtop")
    traj = os.path.join(fdir, "data", "regular", "22_traj.netcdf")
    mdtraj_traj = md.load_frame(traj, 0, top=mol)
    mdtraj_traj = mdtraj_traj.atom_slice(
        mdtraj_traj.topology.select("protein")
    )
    return mdtraj_traj


@pytest.fixture
def example_mdtraj_traj():
    fdir = os.path.dirname(os.path.abspath(__file__))
    top_file = os.path.join(fdir, "data", "regular", "22_solvated.prmtop")
    traj = os.path.join(fdir, "data", "regular", "22_traj.netcdf")
    mdtraj_traj = md.load(traj, top=top_file)
    mdtraj_traj = mdtraj_traj.atom_slice(
        mdtraj_traj.topology.select("protein")
    )
    return mdtraj_traj


def test_getDih(example_mdtraj_frame):
    example = src.dihedrals.getReducedDihedrals(example_mdtraj_frame)
    assert np.allclose(
        example,
        np.array(
            [
                0.53243357,
                0.596079,
                -0.22253884,
                -0.95667714,
                0.34100965,
                -0.8464718,
                -0.8029258,
                -0.97492385,
                -0.29115087,
                -0.9400598,
                -0.9995696,
                0.9848278,
                0.94615984,
                -0.36416814,
                0.7476457,
                -0.02933573,
                -0.17353444,
                -0.32369983,
                -0.93133324,
                -0.66409785,
                -0.9625979,
                -0.99595374,
                -0.9602315,
                -0.93949574,
                -0.9227731,
                0.270934,
                0.08986745,
                -0.2792052,
                0.34256056,
                0.38534364,
            ]
        ),
    )


def test_make_PCA(example_mdtraj_traj):
    pca_object, transformed_output = src.pca.make_PCA(
        example_mdtraj_traj, "dihedral"
    )
    assert type(pca_object) is PCA
    assert np.allclose(
        pca_object.fit_transform(
            src.dihedrals.getReducedDihedrals(example_mdtraj_traj)
        ),
        transformed_output,
    )


def test_plot_PCA():
    input = np.array([[1, 0], [1, 1], [0, 0]])
    fig_axis = src.pca.plot_PCA(input, "cartesian", "compound")
    fig, ax = plt.subplots()
    assert type(fig_axis) == type(ax)


def test_plot_PCA_citra():
    input = np.array([[1, 0], [1, 1], [0, 0]])
    fig_axis = src.pca.plot_PCA_citra(input, input, "cartesian", "Compound 22")
    fig, ax = plt.subplots()
    assert type(fig_axis) == type(ax)
