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
def test_md_pre_topology():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/md_pre_topology/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_pre_topology/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/H2O/0_pre_topology/4amb_noH.pdb",
            file=sys.stderr,
        )

        # Copy other file requirements to the temporary workdir.
        shutil.copyfile("samples_tests.tsv", f"{workdir}/samples_tests.tsv")

        # Add local modules to the python path.
        os.environ["PYTHONPATH"] = (
            os.path.abspath(".") + os.pathsep + os.path.abspath("./libs/")
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/interim/tests/24/H2O/0_pre_topology/4amb_noH.pdb",
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
            "data/interim/tests/24/H2O/0_pre_topology/4amb_noH_nonprot.pdb",
            "data/interim/tests/24/H2O/0_pre_topology/4amb_noH_renum.txt",
            "data/interim/tests/24/H2O/0_pre_topology/4amb_noH_sslink",
            "data/interim/tests/24/H2O/0_pre_topology/log.log",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
