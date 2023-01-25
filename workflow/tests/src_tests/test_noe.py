"""
Test noe functions.

"""

import pytest
import src.noe
import numpy as np


@pytest.mark.parametrize(
    "input, weight, expected_NOE_average",
    [
        (np.array([0]), None, 0),
        (np.array([0, 0]), None, 0),
        (np.array([1.0, 1.0, 1.0]), None, 1.0),
        (np.array([1, 2, 3]), None, 1.197568),
        (np.array([1, 2, 3]), np.array([1, 1, 1]), 1.197568),
    ],
)
def test_NOE_average(input, weight, expected_NOE_average):
    assert src.noe.NOE_average(input, weight) == pytest.approx(
        expected_NOE_average
    )


# @pytest.mark.parametrize("rw_type",[
#     0,
#     1,
#     2,
# ])
# def test_NOE_mdtraj(mdtraj_traj, rw_type):
#     assert src.noe.NOE_mdtraj(mdtraj_traj_protein, reweigh_type=rw_type) == \
# pytest.approx(1.197568)
