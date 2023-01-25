import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_rdkit_confgen():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/rdkit_confgen/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/rdkit_confgen/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/rdkit/ETKDGv3mmff/mcs.pdb data/interim/tests/24/rdkit/ETKDGv3mmff/conf_energies.txt",
            file=sys.stderr,
        )

        # Copy other file requirements to the temporary workdir.
        shutil.copyfile("samples_tests.tsv", f"{workdir}/samples_tests.tsv")

        # Add local modules to the python path.
        os.environ["PYTHONPATH"] = (
            os.path.abspath(".")
            + os.pathsep
            + os.path.abspath("./libs/")
            + os.pathsep
            + os.path.abspath("./libs/rdkit")
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/interim/tests/24/rdkit/ETKDGv3mmff/mcs.pdb",
                "data/interim/tests/24/rdkit/ETKDGv3mmff/conf_energies.txt",
                "-f",
                "-j1",
                "--keep-target-files",
                "--configfile",
                "snakemake-config_tests.yaml",
                "--use-conda",
                "--directory",
                workdir,
                "--conda-not-block-search-path-envvars",
            ]
        )

        to_delete = [
            "samples_old.tsv",
            "samples_tests.tsv",
            "data/interim/tests/24/rdkit/ETKDGv3mmff/rdkit_confgen_log.py.ipynb",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker_present(data_path, expected_path, workdir).check()
