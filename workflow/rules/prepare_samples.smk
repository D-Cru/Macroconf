# Read in tabular configuration from samples.tsv and process into samples dict.

import pandas as pd
import hashlib
import os
import numpy as np
import src.samples
import src.utils


# Run all samples from the file samples.tsv
# Get samples and transform
samples = pd.read_table(config["sample_file"], dtype=str)
columns = list(samples.columns.values)

# Expand samples to account for multiple replicates
for col in columns:
    if col != "repeats":
        samples[col] = samples[col].str.split(",")
samples["repeats"] = samples["repeats"].astype(int)
samples["repeats"] = samples["repeats"].apply(range)

# Expand sample rows if multiple parameters are given
for col in columns:
    samples = samples.explode(col)

# Set data type for each column
samples["simtime"] = samples["simtime"].astype(int)
samples["dt"] = samples["dt"].astype(float)


# Compute nstlim from simtime and dt
samples["nstlim"] = samples.apply(
    lambda x: src.samples.compute_nstlim(x["dt"], x["simtime"]), axis=1
)

# Replace 'native' by the reported solvent used in the experiment


def replace_solvent(compound, solvent):
    if solvent == "native":
        # print(f"replace native with {lookup_exp_solvent(compound)}, C: {compound}")
        exp_solvent = src.samples.lookup_exp_solvent(
            compound, config["dataset_file"]
        )
        if exp_solvent == "DMSO":
            exp_solvent = config["DMSO_default"]
        elif exp_solvent == "CDCl3":
            exp_solvent = config["CDCl3_default"]
        elif exp_solvent == "Chloroform":
            exp_solvent = config["CDCl3_default"]
        elif exp_solvent == "H2O":
            exp_solvent = "H2O"
        else:
            raise ValueError(f"Unknown solvent {exp_solvent}")
        return exp_solvent
    else:
        return solvent


samples["solvent"] = samples.apply(
    lambda x: replace_solvent(x["compound"], x["solvent"]), axis=1
)

samples = samples.astype(str)

# Replace nans with default values, check for other nans in columns where they don't belong
samples = src.samples.adjust_method_specific_parameters(
    samples, config["method_defaults"]
)

samples = samples.reset_index(drop=True)
samples = samples.drop_duplicates()

for col in columns:
    samples[col] = samples[col].str.strip()

# Compute hash for each sample set
samples = src.samples.compute_hash(samples, columns)
samples = samples.drop_duplicates()

# print("The following samples will be run:")
# pd.set_option("display.max_rows", None, "display.max_columns", None)
# print(samples)

# Read in history file, concat, drop duplicates and write hist file again.
# Like this, we don't overwrite old runs with different parameters
src.samples.update_hist_file(samples, config["sample_output"])

# produce list of all specified samples.
outputs = expand(
    "data/processed/{{exp_name}}/results/{compound}/{solvent}/{md_method}/{sim_time}/{repeat}/{index}_{{file}}",
    zip,
    compound=samples["compound"],
    md_method=samples["method"],
    sim_time=samples["simtime"],
    solvent=samples["solvent"],
    repeat=samples["repeats"],
    index=samples.index,
)

# List of single MD run analysis results
if bool(config["run_single_md_analysis"]):
    single_md_analysis_files = (
        expand(outputs, file="pca_dihed.png", exp_name=config["exp_name"]),
    )
else:
    single_md_analysis_files = []

# Produce list of files from various analysis steps, if requested in config file
# Read in manual method comparison file
hash_list = np.array(config["hash_list"])
if bool(config["run_comp_analysis"]):
    comp_analysis_files = expand(
        expand(
            "data/processed/{exp_name}/results/{{compound_dir}}/comparison/{{index_0}}_{{index_1}}_{{index_2}}-{{confgen}}_{{mode}}-cluster_pca.svg",
            exp_name=config["exp_name"],
        ),
        zip,
        index_0=hash_list[:, 0],
        index_1=hash_list[:, 1],
        index_2=hash_list[:, 2],
        compound_dir=hash_list[:, 3],
        confgen=hash_list[:, 4],
        mode=hash_list[:, 5],
    )
else:
    comp_analysis_files = []

if bool(config["run_heatmap_analysis"]):
    heatmap_files = expand(
        "data/processed/{exp_name}/results/{compound}/comparison/heatmap.png",
        compound=config["heatmap_compounds"],
        exp_name=config["exp_name"],
    )
else:
    heatmap_files = []

if not bool(config["run_method_comp"]):
    method_comp_files = []

# Generate convergence check files, if requested
if bool(config["run_convergence_check"]):
    convergence_check_files = []
    for hs in config["convergence_check"]:
        sample_details = samples.loc[hs]
    convergence_check_files.append(
        f"data/processed/conv_check/results/{sample_details.compound}/{sample_details.solvent}/{sample_details.method}/{sample_details.simtime}/{sample_details.repeats}/{hs}_comparison-plot.png",
    )
else:
    convergence_check_files = []


# Generate equilibration check files, if requested for every run performed
if bool(config["run_eq_analysis"]):
    equilibration_check_files = []
    # Get equilibration files for all runs
    for hs in samples.index:
        sample_details = samples.loc[hs]
        equilibration_check_files.append(
            f"data/processed/{config['exp_name']}/results/{sample_details.compound}/{sample_details.solvent}/em_eq/eq_3.png"
        )
    # remove duplicates, since equilibrations can be shared between runs
    equilibration_check_files = list(set(equilibration_check_files))
else:
    equilibration_check_files = []


# Generate files for cheminformatics conformer generator runs
conf_gen_files = []
# OMEGA:
if bool(config["run_omega"]):
    compounds = list(set(samples.compound.tolist()))
    conf_gen_files.extend(
        expand(
            "data/interim/{exp_name}/{compound}/{conf_gen}/{conf_type}/mcs.pdb",
            exp_name=config["exp_name"],
            compound=compounds,
            conf_gen="omega",
            conf_type=config["confgen_parameters"]["omega"],
        )
    )
# RDKIT:
if bool(config["run_rdkit"]):
    compounds = list(set(samples.compound.tolist()))
    conf_gen_files.extend(
        expand(
            "data/interim/{exp_name}/{compound}/{conf_gen}/{conf_type}/mcs.pdb",
            exp_name=config["exp_name"],
            compound=compounds,
            conf_gen="rdkit",
            conf_type=config["confgen_parameters"]["rdkit"],
        )
    )

# Generate files for cheminformatics conformer generator NOE analysis
conf_gen_NOE = []
if bool(config["run_cheminfo_NOE_analysis"]):
    conf_gen_NOE.extend(
        expand(
            "data/processed/{exp_name}/results/{compound}/conf_gen/{conf_gen}/{conf_type}/best_NOE.png",
            exp_name=config["exp_name"],
            compound=compounds,
            conf_gen="omega",
            conf_type=config["confgen_parameters"]["omega"],
        )
    )
    conf_gen_NOE.extend(
        expand(
            "data/processed/{exp_name}/results/{compound}/conf_gen/{conf_gen}/{conf_type}/best_NOE.png",
            exp_name=config["exp_name"],
            compound=compounds,
            conf_gen="rdkit",
            conf_type=config["confgen_parameters"]["rdkit"],
        )
    )

# Generate figures for paper

if bool(config["make_paper_figures"]):
    # paper_figures.csv columns: 'label,file_path_name'
    paper_fig_df = pd.read_csv("reports/paper_figures.csv")
    # remove path from file name
    paper_fig_df["file_names"] = [
        s.split("/")[-1] for s in paper_fig_df.file_path_name.tolist()
    ]
    # combine paper_fig_files and paper_fig_labels into one list of concatenated strings
    paper_fig_df["files_labeled"] = [
        f"reports/paper_figures/{x}_{y}"
        for x, y in zip(
            paper_fig_df.label.tolist(), paper_fig_df["file_names"]
        )
    ]
    paper_fig_files = paper_fig_df.files_labeled.tolist()
else:
    paper_fig_files = []
# TODO: add 19,20, (bokeh output)


# Update TOC file with all relevant notebooks for jupyterbook

if bool(config["update_jupyter_book"]):

    # Generate symlinks
    import os
    from pathlib import Path
    import src.utils
    import yaml

    # Clear symlink directory
    src.utils.remove_symlinks(config["jb-links"])

    # Dataset notebook
    dataset_notebook = os.path.join(config["jb-links"], "dataset.ipynb")
    src.utils.generate_symlink(
        "../../../notebooks/visualize_MacroConf.ipynb", dataset_notebook
    )

    # Single notebooks
    single_notebooks = expand(
        expand(
            "../../../data/processed/{exp_name}/notebooks/{{compound}}/{{solvent}}/{{md_method}}/{{sim_time}}/{{repeat}}/{{index}}_{{md_method}}_processed.ipynb",
            exp_name=config["exp_name"],
        ),
        zip,
        compound=samples["compound"],
        exp_name=config["exp_name"],
        md_method=samples["method"],
        sim_time=samples["simtime"],
        solvent=samples["solvent"],
        repeat=samples["repeats"],
        index=samples.index,
    )
    single_notebook_titles = [s.split(os.sep)[7:12] for s in single_notebooks]
    for nb in single_notebooks:
        # Generate symlinks that reside in reports/jb/links/ and point to the notebooks in data/processed/../notebooks
        src.utils.generate_symlink(
            nb, f"{config['jb-links']}{os.path.basename(nb)}"
        )
    # Comparison notebooks

    comp_notebooks = expand(
        expand(
            "../../../data/processed/{exp_name}/notebooks/{{compound_dir}}/{{index_0}}_{{index_1}}_{{index_2}}_{{confgen}}_{{mode}}_compar.ipynb",
            exp_name=config["exp_name"],
        ),
        zip,
        index_0=hash_list[:, 0],
        index_1=hash_list[:, 1],
        index_2=hash_list[:, 2],
        compound_dir=hash_list[:, 3],
        confgen=hash_list[:, 4],
        mode=hash_list[:, 5],
    )

    comp_notebook_titles = [
        f"Compound {c.split(os.sep)[7]}" for c in comp_notebooks
    ]

    for nb in comp_notebooks:
        # Generate symlinks that reside in reports/jb/links/ and point to the notebooks in data/processed/../notebooks
        src.utils.generate_symlink(
            nb, f"{config['jb-links']}{os.path.basename(nb)}"
        )

    # Full comparison notebook
    full_comp_notebooks = config["full_analysis_files"]

    for nb in full_comp_notebooks:
        # Generate symlinks that reside in reports/jb/links/ and point to the notebooks in data/processed/../notebooks
        src.utils.generate_symlink(
            nb, f"{config['jb-links']}{os.path.basename(nb)}"
        )

    # Generate Jupyterbook TOC file
    toc = {
        "format": "jb-book",
        "root": "index",
        "options": {"numbered": True},
        "chapters": [
            {"file": "links/dataset"},
            {
                "file": "single",
                "sections": [
                    {
                        "file": f"links/{Path(s).stem}",
                        "title": f"{'-'.join(t)}",
                    }
                    for s, t in zip(single_notebooks, single_notebook_titles)
                ],
            },
            {
                "file": "manual_comparison",
                "sections": [
                    {"file": f"links/{Path(s).stem}", "title": f"{t}"}
                    for s, t in zip(comp_notebooks, comp_notebook_titles)
                ],
            },
            {"file": "method_comparison"},
            {
                "file": "full_comparison",
                "sections": [
                    {"file": f"links/{Path(s).stem}"}
                    for s in full_comp_notebooks
                ],
            },
        ],
    }
    with open(config["jb-toc"], "w") as f:
        docs = yaml.dump(toc, f, default_flow_style=False, sort_keys=False)
