"""
Test dihedrals functions.

"""
import pytest
import src.dihedrals
import numpy as np


@pytest.mark.parametrize(
    "phi, psi, expected_motive",
    [
        (np.array([-50]), np.array([-50]), ["Λ"]),
        (np.array([120]), np.array([-140]), ["β"]),
    ],
)
def test_miao_ramachandran(phi, psi, expected_motive):
    assert src.dihedrals.miao_ramachandran(phi, psi) == expected_motive
