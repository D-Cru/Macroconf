"""
Test pca functions.

"""

import pytest
import numpy as np
import pandas as pd
import src.samples


def test_nan_columns():
    config_dict = {
        "igamd": {"methods": ["GaMD"], "default": "3", "others": "nan"},
        "test": {
            "methods": ["aMD", "GaMD"],
            "default": "3xyz",
            "others": "nan",
        },
        "test2": {
            "methods": ["cMD", "aMD", "GaMD"],
            "default": "2.5",
            "others": "nan",
        },
    }

    columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
        "nstlim",
        "test",
        "test2",
    ]
    data = [
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "nan",
            "2",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "nan",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "nan",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "nan",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2.5",
        ],
    ]

    result = [
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "nan",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "nan",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "2.5",
        ],
    ]

    samples = pd.DataFrame(data, columns=columns)
    reference = pd.DataFrame(result, columns=columns)

    with pytest.raises(Exception) as e_info:
        samples = src.samples.adjust_method_specific_parameters(
            samples, config_dict
        )


def test_default_value_replacements():

    config_dict = {
        "igamd": {"methods": ["GaMD"], "default": "3", "others": "nan"},
        "test": {
            "methods": ["aMD", "GaMD"],
            "default": "3xyz",
            "others": "nan",
        },
        "test2": {
            "methods": ["cMD", "aMD", "GaMD"],
            "default": "2.5",
            "others": "nan",
        },
    }

    columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
        "nstlim",
        "test",
        "test2",
    ]
    data = [
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "nan",
            "2",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "nan",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2.5",
        ],
    ]

    result = [
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "3",
            "1000000000",
            "1xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "3.5",
        ],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "2",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "cMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "nan",
            "2",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "2xyz",
            "2.5",
        ],
        [
            "22",
            "aMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
            "3xyz",
            "2.5",
        ],
    ]

    samples = pd.DataFrame(data, columns=columns)
    reference = pd.DataFrame(result, columns=columns)

    samples = src.samples.adjust_method_specific_parameters(
        samples, config_dict
    )
    # with pytest.raises(Exception) as e_info:
    pd.testing.assert_frame_equal(samples, reference)


def test_hashes_drop():
    """Test that the hash is the same when dropping a column with all nan"""

    columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
        "nstlim",
    ]
    hash_columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
    ]
    data = [
        ["22", "cMD", "H2O", "2000", "0.002", "0", "0.2", "nan", "1000000000"],
        ["22", "aMD", "H2O", "2000", "0.002", "0", "0.2", "nan", "1000000000"],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
        ],
    ]

    df = pd.DataFrame(data, columns=columns)
    hashes = src.samples.compute_hash(df, hash_columns).index.values

    # drop column igamd from df
    df = df.drop(columns=["igamd"])
    hash_columns_drop = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
    ]
    hashes_drop = src.samples.compute_hash(df, hash_columns_drop).index.values
    np.testing.assert_array_equal(hashes, hashes_drop)


def test_hashes_drop_extended():
    """Test that the hash is the same when dropping a column with all nan"""

    columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
        "nstlim",
    ]
    hash_columns = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
        "igamd",
    ]
    data = [
        ["22", "cMD", "H2O", "2000", "0.002", "0", "0.2", "nan", "1000000000"],
        ["22", "aMD", "H2O", "2000", "0.002", "0", "0.2", "nan", "1000000000"],
        [
            "22",
            "GaMD",
            "H2O",
            "2000",
            "0.002",
            "0",
            "0.2",
            "nan",
            "1000000000",
        ],
        ["22", "GaMD", "H2O", "2000", "0.002", "0", "0.2", "3", "1000000000"],
    ]

    df = pd.DataFrame(data, columns=columns)
    hashes = src.samples.compute_hash(df, hash_columns).index.values

    # drop column igamd from df
    df = df.drop(columns=["igamd"])
    hash_columns_drop = [
        "compound",
        "method",
        "solvent",
        "simtime",
        "dt",
        "repeats",
        "other_param",
    ]
    hashes_drop = src.samples.compute_hash(df, hash_columns_drop).index.values
    # Check that igamd nan column hashes are the same
    np.testing.assert_array_equal(hashes[0:2], hashes_drop[0:2])
    # check that hash changes when igamd is not nan, but when column is dropped
    assert hashes[3] != hashes_drop[3]
    assert hashes_drop[2] == hashes_drop[3]
    assert hashes[2] == hashes_drop[3]
