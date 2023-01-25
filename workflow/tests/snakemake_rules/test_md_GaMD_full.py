import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import common


@pytest.mark.skip(
    reason="no way of currently testing this, since no Amber available on CI"
)
def test_md_GaMD_full():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/md_GaMD_full/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_GaMD_full/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md.out data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md-1.rst data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md-1.info data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_traj.netcdf data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_gamd.log data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_gamd-restart.dat",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md.out data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md-1.rst data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_md-1.info data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_traj.netcdf data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_gamd.log data/interim/tests/24/H2O/11_GaMD_full/2/0/7f5696daa1cbe87d_gamd-restart.dat",
                "-f",
                "-j1",
                "--keep-target-files",
                "--configfile",
                "workflow/snakemake-config_tests.yaml",
                "--use-conda",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
