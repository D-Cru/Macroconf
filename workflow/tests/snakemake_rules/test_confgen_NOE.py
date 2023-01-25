import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_confgen_NOE():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/confgen_NOE/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/confgen_NOE/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/processed/tests/results/24/conf_gen/omega/basic/best_NOE.png data/processed/tests/results/24/conf_gen/omega/basic/NOE_distribution.png data/processed/tests/results/24/conf_gen/omega/basic/NOE_fulfilled.json data/processed/tests/results/24/conf_gen/omega/basic/bundle_plot.png",
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
                "data/processed/tests/results/24/conf_gen/omega/basic/best_NOE.png",
                "data/processed/tests/results/24/conf_gen/omega/basic/NOE_distribution.png",
                "data/processed/tests/results/24/conf_gen/omega/basic/NOE_fulfilled.json",
                "data/processed/tests/results/24/conf_gen/omega/basic/bundle_plot.png",
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
        # print(test_run.stdout)
        # print(test_run.stderr)
        # shutil.copyfile(f"{workdir}/*", "test_workdir/")
        # Delete the files a, b, c in the workdir.
        to_delete = [
            "data/processed/tests/notebooks/24/conf_gen/omega_basic_NOE.py.ipynb",
            "samples_old.tsv",
            "samples_tests.tsv",
            "data/processed/tests/results/24/conf_gen/omega/basic/bundle_plot.png",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        # common.OutputChecker(data_path, expected_path, workdir).check()
        common.OutputChecker_smart(data_path, expected_path, workdir).check()
