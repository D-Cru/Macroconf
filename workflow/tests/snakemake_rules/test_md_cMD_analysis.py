import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_md_cMD_analysis():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake_rules/md_cMD_analysis/data")
        expected_path = PurePosixPath(
            "tests/snakemake_rules/md_cMD_analysis/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy other file requirements to the temporary workdir.
        shutil.copyfile("samples_tests.tsv", f"{workdir}/samples_tests.tsv")

        # Add local modules to the python path.
        os.environ["PYTHONPATH"] = (
            os.path.abspath(".") + os.pathsep + os.path.abspath("./libs/")
        )

        # dbg
        print(
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_pca_dihed.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_pmf.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_plot.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_pca.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_min_samp.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_time.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_structs.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_stat_plot.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/clusters.pdb data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/clusters_solvated.pdb data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_result.json data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_stats.json data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dihedrals.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dPCA.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dPCA_weights_MC.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NOE_dist.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_multiple.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/rst/done.done data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape.png data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape_weights.dat data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_overview.svg",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_pca_dihed.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_pmf.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_plot.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_pca.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_min_samp.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_time.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_cluster_structs.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_stat_plot.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/clusters.pdb",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/clusters_solvated.pdb",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_result.json",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_noe_stats.json",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dihedrals.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dPCA.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dPCA_weights_MC.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NOE_dist.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_multiple.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/rst/done.done",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape.png",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_NPR_shape_weights.dat",
                "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_overview.svg",
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
            "data/processed/tests/notebooks/24/H2O/cMD/2/0/214b8ec88026de74_cMD_processed.ipynb",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/pym.pml",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/rst/rst_700.rst",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_clusters/test_0001.png",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dPCA_weights_MC_short.dat",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_dihedrals_short.dat",
            "data/processed/tests/results/24/H2O/cMD/2/0/214b8ec88026de74_traj_short.netcdf",
        ]
        for d in to_delete:
            os.remove(f"{workdir}/{d}")

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker_smart(data_path, expected_path, workdir).check()
