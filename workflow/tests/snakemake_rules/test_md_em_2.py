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
def test_md_em_2():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/md_em_2/data")
        expected_path = PurePosixPath("tests/snakemake_rules/md_em_2/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/H2O/3_minim_2/em_2.out data/interim/tests/24/H2O/3_minim_2/em_2.rst data/interim/tests/24/H2O/3_minim_2/mdinfo.info data/interim/tests/24/H2O/3_minim_2/traj.ncdf",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/interim/tests/24/H2O/3_minim_2/em_2.out data/interim/tests/24/H2O/3_minim_2/em_2.rst data/interim/tests/24/H2O/3_minim_2/mdinfo.info data/interim/tests/24/H2O/3_minim_2/traj.ncdf",
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
