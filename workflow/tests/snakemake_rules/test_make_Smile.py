import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_make_Smile():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/make_Smile/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/make_Smile/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/interim/tests/24/data.json data/interim/tests/24/smile.smi data/interim/tests/24/structure.png",
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
                "data/interim/tests/24/data.json",
                "data/interim/tests/24/smile.smi",
                "data/interim/tests/24/structure.png",
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
            # "data/interim/tests/24/make_smile.log",
            "samples_old.tsv",
            "samples_tests.tsv",
            "data/interim/tests/24/make_smile.log",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker_smart(data_path, expected_path, workdir).check()
