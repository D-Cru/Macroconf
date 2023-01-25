"""
Test reweighting functions.

"""

import pytest
import os
import src.pyreweight
import numpy as np


def get_path(filename):
    fdir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(fdir, "data", "pyreweight_test", filename)
    return input_path


def get_file(filename):
    input = np.loadtxt(get_path(filename))
    return input


@pytest.mark.parametrize(
    "input_2d_file, weight_file, result_2d_file, args",
    [
        (
            "Phi_Psi",
            "weights.dat",
            "pmf-2D-Phi_Psi-reweight-CE2.xvg",
            {
                "job": "amdweight_CE",
                "Xdim": [-180, 180],
                "discX": str(6),
                "Emax": 20,
                "cutoff": 10,
            },
        ),
        (
            "Phi_Psi",
            "weights.dat",
            "pmf-2D-Phi_Psi-reweight-MC.xvg",
            {
                "job": "amdweight_MC",
                "Xdim": [-180, 180],
                "discX": str(6),
                "Emax": 20,
                "cutoff": 10,
            },
        ),
    ],
)
def test_pyreweight_2d(input_2d_file, weight_file, result_2d_file, args):
    from src.utils import dotdict
    import numpy as np

    # Original input command
    # python PyReweighting-2D.py -cutoff 10 -input Phi_Psi -Xdim -180 180
    # -discX 6 -Ydim -180 180 -discY 6 -Emax 20 -job amdweight_CE
    # -weight weights.dat

    # now do reweighting...
    # possible arguments are:
    # args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff,
    #  T, Emax, fit, order}

    # -job amdweight_MC -weight weights.dat -Xdim -3 3 -discX 0.1 -Ydim -3 3
    # -discY 0.1 -Emax 20
    args["weight"] = get_path(weight_file)
    arguments = dotdict(args)
    hist, binsX, binsY = src.pyreweight.pyreweight_2d(
        get_file(input_2d_file), arguments
    )

    result = []

    for jx in range(len(hist[:, 0])):
        for jy in range(len(hist[0, :])):
            result.append([binsX[jx], binsY[jy], hist[jx, jy]])
    result = np.array(result)
    assert np.allclose(get_file(result_2d_file), result)


@pytest.mark.parametrize(
    "input_1d_file, weight_file, result_1d_file, args",
    [
        (
            "Psi.dat",
            "weights.dat",
            "pmf-Psi-reweight-CE2.xvg",
            {
                "job": "amdweight_CE",
                "Xdim": [-180, 180],
                "discX": str(6),
                "Emax": 20,
                "cutoff": 10,
            },
        ),
        (
            "Psi.dat",
            "weights.dat",
            "pmf-Psi-reweight-MC.xvg",
            {
                "job": "amdweight_MC",
                "Xdim": [-180, 180],
                "discX": str(6),
                "Emax": 20,
                "cutoff": 10,
            },
        ),
    ],
)
def test_pyreweight_1d(input_1d_file, weight_file, result_1d_file, args):
    from src.utils import dotdict
    import numpy as np

    # Original input command
    # python PyReweighting-1D.py -input Psi.dat -cutoff 10 -Xdim -180 180
    # -disc 6 -Emax 20 -job amdweight_CE -weight weights.dat

    # now do reweighting...
    # possible arguments are:
    # args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff,
    #  T, Emax, fit, order}
    # -job amdweight_MC -weight weights.dat -Xdim -3 3 -discX 0.1 -Ydim -3 3
    # -discY 0.1 -Emax 20
    # arguments={"job": "amdweight_CE", "weight": example_weight_path,
    # "Xdim":[-180, 180], "discX":str(6), "Emax":20, "cutoff":10}
    args["weight"] = get_path(weight_file)
    arguments = dotdict(args)
    hist, binsX, *_ = src.pyreweight.pyreweight_1d(
        get_file(input_1d_file), arguments
    )

    result = []

    for j in range(len(hist[:])):
        result.append([binsX[j], hist[j]])

    result = np.array(result)
    assert np.allclose(get_file(result_1d_file), result)
