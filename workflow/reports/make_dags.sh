#!/bin/sh
# Create various DAGs for the documentation

# Go to the root of the project
BASEDIR=$(dirname $0)
cd ${BASEDIR}
cd .. 

# Set output directory
report_dir="reports/dags"
mkdir -p ${report_dir}

echo "Creating DAGS"
# Create DAG (rulegraph) showing the topology rules:
snakemake --rulegraph data/interim/refactor-test/22/H2O/1_make_topology/mc_sol.pdb | dot -Tpdf > $report_dir/dag_make_topology.pdf

# eq/em rules
snakemake --rulegraph data/processed/refactor-test/results/22/H2O/em_eq/eq_3.png | dot -Tpdf > $report_dir/dag_em_eq.pdf

# cMD rules
snakemake --rulegraph data/processed/refactor-test/results/22/H2O/cMD/100/0/282f14929802f7ac_pca_dihed.png | dot -Tpdf > $report_dir/dag_cMD.pdf

# aMD rules
snakemake --rulegraph data/processed/refactor-test/results/22/H2O/aMD/100/0/b3332ce08307c920_pca_dihed.png | dot -Tpdf > $report_dir/dag_aMD.pdf

# GaMD rules
snakemake --rulegraph data/processed/refactor-test/results/22/H2O/GaMD/100/0/10236f52ef18b2a3_pca_dihed.png | dot -Tpdf > $report_dir/dag_GaMD.pdf

# Omega rules
snakemake --rulegraph data/processed/refactor-test/results/22/conf_gen/omega/basic/best_NOE.png | dot -Tpdf > $report_dir/dag_Omega.pdf

# RDKit rules
snakemake --rulegraph data/processed/refactor-test/results/22/conf_gen/rdkit/ETKDGv3mmff/best_NOE.png | dot -Tpdf > $report_dir/dag_RDKit.pdf

# Full dag
snakemake --rulegraph | dot -Tpdf > $report_dir/dag_full.pdf