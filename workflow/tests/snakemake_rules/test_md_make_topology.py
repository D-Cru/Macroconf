import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import common

# @pytest.mark.skip(reason="no way of currently testing this, since no Amber available on CI")
def test_md_make_topology():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(
            "tests/snakemake_rules/md_make_topology/data"
        )
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_make_topology/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy other file requirements to the temporary workdir.
        shutil.copyfile("samples_tests.tsv", f"{workdir}/samples_tests.tsv")
        os.makedirs(
            os.path.dirname(
                f"{workdir}/libs/forcefields/leaprc.protein.ff14SB_noterminal"
            ),
            exist_ok=True,
        )
        shutil.copyfile(
            "libs/forcefields/leaprc.protein.ff14SB_noterminal",
            f"{workdir}/libs/forcefields/leaprc.protein.ff14SB_noterminal",
        )

        # Add local modules to the python path.
        os.environ["PYTHONPATH"] = (
            os.path.abspath(".") + os.pathsep + os.path.abspath("./libs/")
        )

        # dbg
        print(
            "data/interim/tests/24/H2O/1_make_topology/mc_sol.prmtop data/interim/tests/24/H2O/1_make_topology/mc_sol.inpcrd data/interim/tests/24/H2O/1_make_topology/mc_sol.pdb data/interim/tests/24/H2O/1_make_topology/mc_gas.mol2",
            file=sys.stderr,
        )

        # Run the test job.
        try:
            sp.check_output(
                [
                    "python",
                    "-m",
                    "snakemake",
                    "data/interim/tests/24/H2O/1_make_topology/mc_sol.prmtop",
                    "data/interim/tests/24/H2O/1_make_topology/mc_sol.inpcrd",
                    "data/interim/tests/24/H2O/1_make_topology/mc_sol.pdb",
                    "data/interim/tests/24/H2O/1_make_topology/mc_gas.mol2",
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
        except:
            with open(
                f"{workdir}/data/interim/tests/24/H2O/1_make_topology/leap.stdout",
                "r",
            ) as fin:
                print(fin.read())

        to_delete = [
            "samples_old.tsv",
            "samples_tests.tsv",
            "libs/forcefields/leaprc.protein.ff14SB_noterminal",
            "data/interim/tests/24/H2O/1_make_topology/leap.log",
            "data/interim/tests/24/H2O/1_make_topology/leap.stdout",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
