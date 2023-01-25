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
    reason="no way of currently testing this, since no OMEGA available on CI"
)
def test_omega_confgen():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/omega_confgen/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/omega_confgen/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("data/interim/tests/24/omega/basic/mcs.oeb.gz", file=sys.stderr)

        # Copy other file requirements to the temporary workdir.
        shutil.copyfile("samples_tests.tsv", f"{workdir}/samples_tests.tsv")
        os.makedirs(
            os.path.dirname(f"{workdir}/libs/omega/basic.param"), exist_ok=True
        )
        shutil.copyfile(
            "libs/omega/basic.param", f"{workdir}/libs/omega/basic.param"
        )

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
                "data/interim/tests/24/omega/basic/mcs.oeb.gz",
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
            "data/interim/tests/24/omega/_log.txt",
            "data/interim/tests/24/omega/_parm.txt",
            "data/interim/tests/24/omega/_rpt.csv",
            "data/interim/tests/24/omega/basic/log.log",
            "libs/omega/basic.param",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker_smart(data_path, expected_path, workdir).check()
