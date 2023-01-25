import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_md_comp_details():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/md_comp_details/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_comp_details/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-pca_dihe.svg data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-cluster_pca.svg data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-pca_dihe_report.svg data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-shape_comparison.svg data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-cluster_hbonds.svg data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-all_cheminfo_comp.png data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-all_cheminfo_shape.png data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-single_comp_plot.svg",
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
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-pca_dihe.svg",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-cluster_pca.svg",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-pca_dihe_report.svg",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-shape_comparison.svg",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-cluster_hbonds.svg",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-all_cheminfo_comp.png",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-all_cheminfo_shape.png",
                "data/processed/tests/results/24/comparison/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74-omega_basic_rdkit_ETKDGv3mmff-single_comp_plot.svg",
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
            "data/processed/tests/notebooks/24/7f5696daa1cbe87d_a033db1075c6b36c_214b8ec88026de74_omega_basic_rdkit_ETKDGv3mmff_compar.ipynb",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker_smart(data_path, expected_path, workdir).check()
