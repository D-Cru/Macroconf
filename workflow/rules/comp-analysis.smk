# Rules to compare different MD simulation runs
def compare_inputs(wildcards):
    """Produce valid inputs for the rule:md_comp_analysis.

    Input are the wildcards. They contain values for index_0/1/2
    Outputs the corresponding input files (trajectory, topology, etc) to perform
    the comparison analysis.
    """
    output = {}
    method_dict = {
        "cMD": "8_cMD",
        "aMD": "10_aMD_prod",
        "GaMD": "11_GaMD_full",
    }
    output[
        "noe"
    ] = "data/interim/{wildcards.exp_name}/{wildcards.compound_dir}/NOE.json".format(
        wildcards=wildcards
    )
    output[
        "parm"
    ] = "data/interim/{wildcards.exp_name}/{wildcards.compound_dir}/data.json".format(
        wildcards=wildcards
    )
    indices = [wildcards.index_0, wildcards.index_1, wildcards.index_2]
    for i, idx in zip(range(3), indices):
        w = samples.loc[idx]
        path = f"data/interim/{wildcards.exp_name}/{w.compound}/{w.solvent}/{method_dict[w.method]}/{w.simtime}/{w.repeats}/{idx}"
        if w.method != "cMD":
            output[f"weights_{i}"] = f"{path}_weights.dat"
        output[f"out_{i}"] = f"{path}_md.out"
        output[f"traj_{i}"] = f"{path}_traj.netcdf"
        output[f"traj_ncdf_{i}"] = f"{path}_traj.ncdf"
        output[
            f"clusters_{i}"
        ] = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{idx}_clusters/clusters.pdb"
        output[
            f"top_{i}"
        ] = f"data/interim/{wildcards.exp_name}/{w.compound}/{w.solvent}/1_make_topology/mc_sol.prmtop"
        if wildcards.confgen != "0":
            output[
                "cheminfoconfs"
            ] = f"data/interim/{wildcards.exp_name}/{w.compound}/{wildcards.confgen}/{wildcards.mode}/mcs.pdb"  # mcs_aligned.pdb"

        path_proc = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{idx}"
        output[f"red_dihe_{i}"] = f"{path_proc}_dihedrals.dat"
        output[f"dPCA_{i}"] = f"{path_proc}_dPCA.dat"
        output[f"dPCA_weights_MC_{i}"] = f"{path_proc}_dPCA_weights_MC.dat"
        output[f"NOE_dist_{i}"] = f"{path_proc}_NOE_dist.dat"
        output[f"multiple_{i}"] = f"{path_proc}_multiple.dat"
        output[f"shape_{i}"] = f"{path_proc}_NPR_shape.dat"
        output[f"shape_weights_MC_{i}"] = f"{path_proc}_NPR_shape_weights.dat"
    return output


def comp_details(wildcards):
    """Produce valid inputs for the rule:md_comp_details.

    Input are the wildcards. They contain values for index_0/1/2
    Outputs the corresponding input files (trajectory, topology, etc) to perform
    the comparison analysis.
    """
    output = {}
    method_dict = {
        "cMD": "8_cMD",
        "aMD": "10_aMD_prod",
        "GaMD": "11_GaMD_full",
    }
    output[
        "noe"
    ] = "data/interim/{wildcards.exp_name}/{wildcards.compound_dir}/NOE.json".format(
        wildcards=wildcards
    )
    output[
        "parm"
    ] = "data/interim/{wildcards.exp_name}/{wildcards.compound_dir}/data.json".format(
        wildcards=wildcards
    )
    indices = [wildcards.index_0, wildcards.index_1, wildcards.index_2]
    for i, hash in enumerate(indices):
        w = samples.loc[hash]
        path = f"data/interim/{wildcards.exp_name}/{w.compound}/{w.solvent}/{method_dict[w.method]}/{w.simtime}/{w.repeats}/{hash}"
        if w.method != "cMD":
            output[f"weights_{i}"] = f"{path}_weights.dat"
        output[f"out_{i}"] = f"{path}_md.out"
        # output[f"traj_{i}"] = f"{path}_traj.netcdf"
        # output[f"traj_ncdf_{i}"] = f"{path}_traj.ncdf"
        output[
            f"clusters_{i}"
        ] = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{hash}_clusters/clusters.pdb"
        output[
            f"top_{i}"
        ] = f"data/interim/{wildcards.exp_name}/{w.compound}/{w.solvent}/1_make_topology/mc_sol.prmtop"

        # Get conformer generators if not set to 0
        if wildcards.confgens != "0_0":
            conf_gens = wildcards.confgens.split("_")

            for idx, (j, k) in enumerate(
                zip(conf_gens[0::2], conf_gens[1::2])
            ):
                if not ((j == k) and j == "0"):
                    output[f"confgen{idx}"] = (
                        f"data/processed/{wildcards.exp_name}/results/methods/{j}-{k}-NOE_fulfilled.json",
                    )
                    output[
                        f"cheminfoconfs{idx}"
                    ] = f"data/interim/{wildcards.exp_name}/{w.compound}/{j}/{k}/mcs_aligned.pdb"

        path_proc = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{hash}"
        output[f"red_dihe_{i}"] = f"{path_proc}_dihedrals.dat"
        output[f"dPCA_{i}"] = f"{path_proc}_dPCA.dat"
        output[f"dPCA_weights_MC_{i}"] = f"{path_proc}_dPCA_weights_MC.dat"
        output[f"NOE_dist_{i}"] = f"{path_proc}_NOE_dist.dat"
        output[f"multiple_{i}"] = f"{path_proc}_multiple.dat"
        output[f"shape_{i}"] = f"{path_proc}_NPR_shape.dat"
        output[f"shape_weights_MC_{i}"] = f"{path_proc}_NPR_shape_weights.dat"
    return output


rule md_comp_details:
    input:
        unpack(comp_details),
    output:
        pca_dihe=report(
            "data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-pca_dihe.svg"
        ),
        cluster_pca=report(
            "data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-cluster_pca.svg"
        ),
        report_pca_comparison="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-pca_dihe_report.svg",
        shape_comparsion="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-shape_comparison.svg",
        cluster_hbonds="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-cluster_hbonds.svg",
        cluster_hbonds_nolabel="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-cluster_hbonds_nolabel.svg",
        all_cheminfo_comp_pca="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-all_cheminfo_comp.png",
        all_cheminfo_comp_shape="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-all_cheminfo_shape.png",
        single_comp_plot="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-{confgens}-single_comp_plot.svg",
    params:
        sample_0=lambda w: samples.loc["{}".format(w.index_0)].to_dict(),
        sample_1=lambda w: samples.loc["{}".format(w.index_1)].to_dict(),
        sample_2=lambda w: samples.loc["{}".format(w.index_2)].to_dict(),
    log:
        notebook="data/processed/{exp_name}/notebooks/{compound_dir}/{index_0}_{index_1}_{index_2}_{confgens}_compar.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/template-analyse-amber-compar.py.ipynb"


def stats_input(wildcards):
    output = expand(
        expand(
            "data/processed/{{exp_name}}/results/{{compound}}/{solvent}/{md_method}/{sim_time}/{repeat}/{index}_noe_stats.json",
            zip,
            md_method=samples[samples["compound"] == f"{wildcards.compound}"][
                "method"
            ],
            sim_time=samples[samples["compound"] == f"{wildcards.compound}"][
                "simtime"
            ],
            solvent=samples[samples["compound"] == f"{wildcards.compound}"][
                "solvent"
            ],
            repeat=samples[samples["compound"] == f"{wildcards.compound}"][
                "repeats"
            ],
            index=samples[
                samples["compound"] == f"{wildcards.compound}"
            ].index,
        ),
        exp_name=wildcards.exp_name,
        compound=wildcards.compound,
    )
    return output


rule md_comp_stats:
    input:
        stats=stats_input,
    output:
        heatmap=report(
            "data/processed/{exp_name}/results/{compound}/comparison/heatmap.png"
        ),
    log:
        notebook="data/processed/{exp_name}/notebooks/{compound}/comp_stats_log.ipynb",
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/comp_stats.py.ipynb"


def comp_methods_input(wildcards):
    output = {}
    if wildcards.igamd == "nan" or wildcards.igamd == "3":
        igamd_parameters = ["nan", "3"]
    else:
        igamd_parameters = wildcards.igamd

    if wildcards.solvent == "native":
        # Native solvent. Find all simulation id's matching the parameters
        # irrespceitve of solvent
        simulation_ids = samples.query(
            "method == @wildcards.method and simtime == @wildcards.simtime and igamd == @igamd_parameters"
        ).index.values.tolist()
        sim_ids = []
        for sim_hash in simulation_ids:
            # Lookup parameters
            sim_params = samples.loc[sim_hash]

            # lookup native solvent
            native_solvent = replace_solvent(sim_params["compound"], "native")
            # if solvent matches native, add to simulation_ids
            if native_solvent == sim_params["solvent"]:
                sim_ids.append(sim_hash)
        simulation_ids = sim_ids
    else:
        # get all samples that match the method description...
        simulation_ids = samples.query(
            "method == @wildcards.method and solvent == @wildcards.solvent and simtime == @wildcards.simtime and igamd == @igamd_parameters"
        ).index.values.tolist()

    # get all outputs
    for sim_hash in simulation_ids:
        w = samples.loc[sim_hash]
        output[
            f"run_{sim_hash}"
        ] = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_noe_stats.json"
        output[
            f"conv_{sim_hash}"
        ] = f"data/processed/{wildcards.exp_name}/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_conv_data.json"

    return output


# Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
rule md_comp_methods:
    input:
        unpack(comp_methods_input),
    output:
        plot=report(
            "data/processed/{exp_name}/results/methods/{method}-{solvent}-{simtime}-{igamd}-NOE.png"
        ),
        data="data/processed/{exp_name}/results/methods/{method}-{solvent}-{simtime}-{igamd}-NOE_fulfilled.json",
        conv_plot=report(
            "data/processed/{exp_name}/results/methods/{method}-{solvent}-{simtime}-{igamd}-conv_plot.svg"
        ),
    log:
        notebook="data/processed/{exp_name}/notebooks/methods/{method}-{solvent}-{simtime}-{igamd}-NOE_method_comp.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/comp_compounds.py.ipynb"


def confgen_comp_methods_input(wildcards):
    output = {}
    compounds = list(set(samples.compound.values.tolist()))
    for c in compounds:
        output[
            f"{c}"
        ] = f"data/processed/{wildcards.exp_name}/results/{c}/conf_gen/{wildcards.confgen}/{wildcards.mode}/NOE_fulfilled.json"
    return output


# Aggregate all conformer generator methods
rule confgen_comp_methods:
    input:
        unpack(confgen_comp_methods_input),
    output:
        data="data/processed/{exp_name}/results/methods/{confgen,\w+}-{mode}-NOE_fulfilled.json",
    log:
        notebook="data/processed/{exp_name}/notebooks/methods/{confgen}-{mode}-NOE_fulfilled.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/comp_compounds_confgen.py.ipynb"


def comp_all_methods_input(wildcards):
    output = {}
    methods = wildcards.methods.split("-")
    native = methods[1::2]
    methods = methods[::2]
    conf_gens = wildcards.conf_gens.split("-")

    # create list of methods to compare
    method_comp = [samples.loc[a] for a in methods]
    for idx, m in enumerate(method_comp):
        if native[idx] == "native":
            solvent = "native"
        else:
            solvent = m["solvent"]
        output[
            f"{idx}"
        ] = f"data/processed/{config['exp_name']}/results/methods/{m.method}-{solvent}-{m.simtime}-{m.igamd}-NOE_fulfilled.json"
    for idx, (j, k) in enumerate(zip(conf_gens[0::2], conf_gens[1::2])):
        if not ((j == k) and j == "0"):
            output[f"confgen{idx}"] = (
                f"data/processed/{config['exp_name']}/results/methods/{j}-{k}-NOE_fulfilled.json",
            )
    return output


rule md_comp_all_methods:
    input:
        unpack(comp_all_methods_input),
    output:
        plot1=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_1.svg"
        ),
        plot1_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_1.svg"
        ),
        plot2=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_2.svg"
        ),
        plot2_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_2.svg"
        ),
        plot3=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_3.svg"
        ),
        plot3_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_3.svg"
        ),
        plot4=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_4.svg"
        ),
        plot4_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_4.svg"
        ),
        plot5=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_5.svg"
        ),
        plot5_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_5.svg"
        ),
        plot6=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_6.svg"
        ),
        plot7=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_7.svg"
        ),
        plot7_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_7.svg"
        ),
        plot8=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_8.svg"
        ),
        plot8_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_8.svg"
        ),
        plot9=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_9.svg"
        ),
        plot9_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_9.svg"
        ),
        plot10=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_10.svg"
        ),
        plot10_sig=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_sig_10.svg"
        ),
        plot_seq_length=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_seq_length.svg"
        ),
        plot_boxplot_seq_length=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_boxplot_seq_length.svg"
        ),
        plot_boxplot_seq_length_bins=report(
            "data/processed/{exp_name}/results/methods/{methods}_{conf_gens}-NOE-all_boxplot_seq_length_bins.svg"
        ),
    log:
        notebook="data/processed/{exp_name}/notebooks/methods/{methods}_{conf_gens}-NOE_method_comp.ipynb",
    conda:
        "../envs/stats_comp.yaml"
    threads: 2
    notebook:
        "../notebooks/comp_methods.py.ipynb"
