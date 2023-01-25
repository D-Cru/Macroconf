import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_md_aMD_cMD_readE():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(
            "tests/snakemake_rules/md_aMD_cMD_readE/data"
        )
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_aMD_cMD_readE/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_cpptraj.in data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_Epot.dat data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_Edih.dat",
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
                "data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_cpptraj.in",
                "data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_Epot.dat",
                "data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_Edih.dat",
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
            "data/interim/tests/24/H2O/9_aMD_cMD/2/0/1b5ac9a7d3a7b8e5_cpptraj.log",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
