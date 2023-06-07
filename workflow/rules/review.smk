# Convergence and stepwise RMSD analysis


# rule md_GaMD_conv_analysis:
#     input:
#         weights=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_weights.dat",
#         out=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md.out",
#         traj=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_traj.netcdf",
#         traj_ncdf=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_traj.ncdf",
#         top=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_sol.prmtop",
#         parm=f"data/interim/{config['exp_name']}/{{compound_dir}}/data.json",
#         noe=f"data/interim/{config['exp_name']}/{{compound_dir}}/NOE.json",
#         ref_mol=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_gas.mol2",
#     output:
#         pca_dihe=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_pca_dihed.svg",
#             category="Compound {compound_dir}",
#             subcategory="GaMD",
#         ),
#         conv_plot=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_conv_plot.svg",
#             category="Compound {compound_dir}",
#             subcategory="GaMD",
#         ),
#         grid_cells=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_grid_cells.svg",
#             category="Compound {compound_dir}",
#             subcategory="GaMD",
#         ),
#         conv_data=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_conv_data.json",
#             category="Compound {compound_dir}",
#             subcategory="GaMD",
#         ),
#         noe_stats=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_stats.json",
#         noe_result=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_result.json",
#     params:
#         cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/",
#         rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
#         method="GaMD",
#         traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
#         weights_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_weights_short.dat",
#         dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
#         dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
#     log:
#         notebook=f"data/processed/review/notebooks/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_GaMD_processed.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 16
#     notebook:
#         "../notebooks/template-analyse-amber-aMD_convergence.py.ipynb"


# rule md_cMD_stepwise_RMSD:
#     input:
#         out=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_md.out",
#         traj=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_traj.netcdf",
#         traj_ncdf=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_traj.ncdf",
#         top=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_sol.prmtop",
#         parm=f"data/interim/{config['exp_name']}/{{compound_dir}}/data.json",
#         noe=f"data/interim/{config['exp_name']}/{{compound_dir}}/NOE.json",
#         ref_mol=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_gas.mol2",
#     output:
#         pca_dihe=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_pca_dihed.svg",
#             category="Compound {compound_dir}",
#             subcategory="cMD",
#         ),
#         conv_plot=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_conv_plot.svg",
#             category="Compound {compound_dir}",
#             subcategory="cMD",
#         ),
#         grid_cells=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_grid_cells.svg",
#             category="Compound {compound_dir}",
#             subcategory="cMD",
#         ),
#         conv_data=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_conv_data.json",
#             category="Compound {compound_dir}",
#             subcategory="cMD",
#         ),
#         noe_stats=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_noe_stats.json",
#         noe_result=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_noe_result.json",
#     params:
#         cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_clusters/",
#         rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
#         method="cMD",
#         traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
#         dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
#         dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
#     log:
#         notebook=f"data/processed/review/notebooks/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_cMD_processed.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 16
#     notebook:
#         "../notebooks/template-analyse-amber-aMD_convergence.py.ipynb"


# def conv_accumulate(wildcards):
#     output = {}
#     if wildcards.igamd == "nan" or wildcards.igamd == "3":
#         igamd_parameters = ["nan", "3"]
#     else:
#         igamd_parameters = wildcards.igamd

#     if wildcards.solvent == "native":
#         # Native solvent. Find all simulation id's matching the parameters
#         # irrespceitve of solvent
#         simulation_ids = samples.query(
#             "method == @wildcards.method and simtime == @wildcards.simtime and igamd == @igamd_parameters"
#         ).index.values.tolist()
#         sim_ids = []
#         for sim_hash in simulation_ids:
#             # Lookup parameters
#             sim_params = samples.loc[sim_hash]

#             # lookup native solvent
#             native_solvent = replace_solvent(sim_params["compound"], "native")
#             # if solvent matches native, add to simulation_ids
#             if native_solvent == sim_params["solvent"]:
#                 sim_ids.append(sim_hash)
#         simulation_ids = sim_ids
#     else:
#         # get all samples that match the method description...
#         simulation_ids = samples.query(
#             "method == @wildcards.method and solvent == @wildcards.solvent and simtime == @wildcards.simtime and igamd == @igamd_parameters"
#         ).index.values.tolist()

#     # get all outputs
#     for sim_hash in simulation_ids:
#         w = samples.loc[sim_hash]
#         output[
#             f"run_{sim_hash}"
#         ] = f"data/processed/review/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_conv_data.json"

#     return output


# # Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
# rule md_GaMD_conv_accumulate:
#     input:
#         unpack(conv_accumulate),
#     output:
#         plot=report(
#             "data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-conv_hist.svg"
#         ),
#     log:
#         notebook="data/processed/review/notebooks/methods/{method}-{solvent}-{simtime}-{igamd}-GaMD_conv_processed.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/review_convergence.py.ipynb"


# rule confgen_NOE_stepwise_RMSD:
#     input:
#         pdb=f"data/interim/{config['exp_name']}/{{compound}}/{{confgen}}/{{mode}}/mcs_aligned.pdb",
#         noe=f"data/interim/{config['exp_name']}/{{compound}}/NOE.json",
#         parm=f"data/interim/{config['exp_name']}/{{compound}}/data.json",
#         energies=f"data/interim/{config['exp_name']}/{{compound}}/{{confgen}}/{{mode}}/conf_energies.txt",
#     output:
#         best_NOE_plot="data/processed/review/results/{compound}/conf_gen/{confgen}/{mode}/best_NOE.svg",
#         NOE_violin_plot="data/processed/review/results/{compound}/conf_gen/{confgen}/{mode}/NOE_distribution.svg",
#         fulfilled="data/processed/review/results/{compound}/conf_gen/{confgen}/{mode}/NOE_fulfilled.json",
#         bundle_plot="data/processed/review/results/{compound}/conf_gen/{confgen}/{mode}/bundle_plot.svg",
#     threads: 1
#     conda:
#         "../envs/stats.yaml"
#     log:
#         notebook="data/processed/review/notebooks/{compound}/conf_gen/{confgen}_{mode}_NOE.py.ipynb",
#     notebook:
#         "../notebooks/confgen_NOE_review.py.ipynb"


# def comp_methods_input_review(wildcards):
#     output = {}
#     if wildcards.igamd == "nan" or wildcards.igamd == "3":
#         igamd_parameters = ["nan", "3"]
#     else:
#         igamd_parameters = wildcards.igamd

#     if wildcards.solvent == "native":
#         # Native solvent. Find all simulation id's matching the parameters
#         # irrespceitve of solvent
#         simulation_ids = samples.query(
#             "method == @wildcards.method and simtime == @wildcards.simtime and igamd == @igamd_parameters"
#         ).index.values.tolist()
#         sim_ids = []
#         for sim_hash in simulation_ids:
#             # Lookup parameters
#             sim_params = samples.loc[sim_hash]

#             # lookup native solvent
#             native_solvent = replace_solvent(sim_params["compound"], "native")
#             # if solvent matches native, add to simulation_ids
#             if native_solvent == sim_params["solvent"]:
#                 sim_ids.append(sim_hash)
#         simulation_ids = sim_ids
#     else:
#         # get all samples that match the method description...
#         simulation_ids = samples.query(
#             "method == @wildcards.method and solvent == @wildcards.solvent and simtime == @wildcards.simtime and igamd == @igamd_parameters"
#         ).index.values.tolist()

#     # get all outputs
#     for sim_hash in simulation_ids:
#         w = samples.loc[sim_hash]
#         output[
#             f"run_{sim_hash}"
#         ] = f"data/processed/review/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_noe_stats.json"

#     return output


# # Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
# rule md_comp_methods_review:
#     input:
#         unpack(comp_methods_input_review),
#     output:
#         plot=report(
#             "data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-NOE.png"
#         ),
#         data="data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-NOE_fulfilled.json",
#     log:
#         notebook="data/processed/review/notebooks/methods/{method}-{solvent}-{simtime}-{igamd}-NOE_method_comp.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/comp_compounds_review.py.ipynb"


# def confgen_comp_methods_input_review(wildcards):
#     output = {}
#     compounds = list(set(samples.compound.values.tolist()))
#     for c in compounds:
#         output[
#             f"{c}"
#         ] = f"data/processed/review/results/{c}/conf_gen/{wildcards.confgen}/{wildcards.mode}/NOE_fulfilled.json"
#     return output


# # Aggregate all conformer generator methods
# rule confgen_comp_methods_review:
#     input:
#         unpack(confgen_comp_methods_input_review),
#     output:
#         data="data/processed/review/results/methods/{confgen,\w+}-{mode}-NOE_fulfilled.json",
#     log:
#         notebook="data/processed/review/notebooks/methods/{confgen}-{mode}-NOE_fulfilled.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/comp_compounds_confgen_review.py.ipynb"


def comp_all_methods_input_review(wildcards):
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
        ] = f"data/processed/review/results/methods/{m.method}-{solvent}-{m.simtime}-{m.igamd}-NOE_fulfilled.json"
    for idx, (j, k) in enumerate(zip(conf_gens[0::2], conf_gens[1::2])):
        if not ((j == k) and j == "0"):
            output[f"confgen{idx}"] = (
                f"data/processed/review/results/methods/{j}-{k}-NOE_fulfilled.json",
            )
    return output


# rule md_comp_all_methods_review:
#     input:
#         unpack(comp_all_methods_input_review),
#     output:
#         plot1=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_1.svg"
#         ),
#         plot1_sig=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_sig_1.svg"
#         ),
#         plot2=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_2.svg"
#         ),
#         plot2_sig=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_sig_2.svg"
#         ),
#         plot_seq_length=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-seq_length_all.svg"
#         ),
#         plot_boxplot_seq_length=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-seq_length_boxplot_all.svg"
#         ),
#         plot_boxplot_seq_length_bins=report(
#             "data/processed/review/results/methods/{methods}_{conf_gens}-seq_length_boxplot_bins_all.svg"
#         ),
#     log:
#         notebook="data/processed/review/notebooks/methods/{methods}_{conf_gens}-NOE_method_comp.ipynb",
#     conda:
#         "../envs/stats_comp.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/comp_methods_review.py.ipynb"


# PSA, SASA and solvation properties, mc reweighting


rule md_GaMD_boltzmann_reweighting:
    input:
        weights=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_weights.dat",
        out=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md.out",
        traj=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_traj.netcdf",
        traj_ncdf=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_traj.ncdf",
        top=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_sol.prmtop",
        parm=f"data/interim/{config['exp_name']}/{{compound_dir}}/data.json",
        noe=f"data/interim/{config['exp_name']}/{{compound_dir}}/NOE.json",
        ref_mol=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_gas.mol2",
    output:
        pca_dihe=report(
            f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_pca_dihed_mc.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_stats=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_stats_bol.json",
        noe_result=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_result_bol.json",
    params:
        cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/",
        rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
        method="GaMD",
        traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
        weights_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_weights_short.dat",
        dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
        dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
    log:
        notebook=f"data/processed/review/notebooks/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_GaMD_processed_mc.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD_mc_sol_prop.py.ipynb"


# rule md_cMD_sol_prop_mc_rew:
#     input:
#         out=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_md.out",
#         traj=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_traj.netcdf",
#         traj_ncdf=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/8_cMD/{{time}}/{{repeat}}/{{index}}_traj.ncdf",
#         top=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_sol.prmtop",
#         parm=f"data/interim/{config['exp_name']}/{{compound_dir}}/data.json",
#         noe=f"data/interim/{config['exp_name']}/{{compound_dir}}/NOE.json",
#         ref_mol=f"data/interim/{config['exp_name']}/{{compound_dir}}/{{solvent}}/1_make_topology/mc_gas.mol2",
#     output:
#         pca_dihe=report(
#             f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_pca_dihed_mc.svg",
#             category="Compound {compound_dir}",
#             subcategory="cMD",
#         ),
#         noe_stats=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_noe_stats_bol.json",
#         noe_result=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_noe_result_bol.json",
#         sasa=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_sasa_mc.pkl",
#         psa=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_psa_mc.pkl",
#         solvation_properties=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_solvation_properties.json",
#     params:
#         cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_clusters/",
#         rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
#         method="cMD",
#         traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
#         dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
#         dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
#     log:
#         notebook=f"data/processed/review/notebooks/{{compound_dir}}/{{solvent}}/cMD/{{time}}/{{repeat}}/{{index}}_cMD_processed_mc.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 16
#     notebook:
#         "../notebooks/template-analyse-amber-aMD_mc_sol_prop.py.ipynb"


def solv_accumulate(wildcards):
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
        ] = f"data/processed/refactor-test/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_solvation_properties.json"

    return output


# Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
rule md_sovlent_props_accumulate:
    input:
        unpack(solv_accumulate),
    output:
        plot=report(
            "data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-solv_plot.svg"
        ),
        data="data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-solv_properties.json",
    log:
        notebook="data/processed/review/notebooks/methods/{method}-{solvent}-{simtime}-{igamd}-solv_properties.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/review_solvation_props.py.ipynb"


def confgen_comp_solvents_review(wildcards):
    output = {}
    compounds = list(set(samples.compound.values.tolist()))
    for c in compounds:
        output[
            f"{c}-sasa"
        ] = f"data/processed/refactor-test/results/{c}/conf_gen/{wildcards.confgen}/{wildcards.mode}/sasa.json"
        output[
            f"{c}-psa"
        ] = f"data/processed/refactor-test/results/{c}/conf_gen/{wildcards.confgen}/{wildcards.mode}/psa.json"
    return output


# # Aggregate all conformer generator methods
# rule confgen_comp_methods_review:
#     input:
#         unpack(confgen_comp_methods_input_review),
#     output:
#         data="data/processed/review/results/methods/{confgen,\w+}-{mode}-NOE_fulfilled.json",
#     log:
#         notebook="data/processed/review/notebooks/methods/{confgen}-{mode}-NOE_fulfilled.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/comp_compounds_confgen_review.py.ipynb"


rule comp_solvent_props:
    input:
        unpack(confgen_comp_solvents_review),
        gamd="data/processed/review/results/methods/GaMD-native-2000-3-solv_properties.json",
        cmd="data/processed/review/results/methods/cMD-native-2000-nan-solv_properties.json",
    output:
        plot_sasa="data/processed/review/results/methods/{confgen,\w+}-{mode}-solv-sasa_plot.svg",
        plot_psa="data/processed/review/results/methods/{confgen,\w+}-{mode}-solv-psa_plot.svg",
    log:
        notebook="data/processed/review/notebooks/methods/{confgen,\w+}-{mode}-solv_properties.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/review_solvation_comp.py.ipynb"


def comp_methods_input_bol(wildcards):
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
        ] = f"data/processed/review/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_noe_stats_bol.json"

    return output


# Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
rule md_comp_methods_review_bol:
    input:
        unpack(comp_methods_input_bol),
    output:
        plot=report(
            "data/processed/review/results/methods/{method}bol-{solvent}-{simtime}-{igamd}-NOE_bol.png"
        ),
        data="data/processed/review/results/methods/{method}bol-{solvent}-{simtime}-{igamd}-NOE_fulfilled_bol.json",
    log:
        notebook="data/processed/review/notebooks/methods/{method}bol-{solvent}-{simtime}-{igamd}-NOE_method_comp_bol.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/comp_compounds_review.py.ipynb"


rule md_comp_reweighting_bol:
    input:
        # unpack(comp_all_methods_input_review),
        gamd_mc="data/processed/refactor-test/results/methods/GaMD-native-2000-3-NOE_fulfilled.json",
        gamd_bol="data/processed/review/results/methods/GaMDbol-native-2000-3-NOE_fulfilled_bol.json",
    output:
        plot1=report(
            "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_1_rw.svg"
        ),
        plot1_sig=report(
            "data/processed/review/results/methods/{methods}_{conf_gens}-NOE-all_sig_1_rw.svg"
        ),
    log:
        notebook="data/processed/review/notebooks/methods/{methods}-{conf_gens}-NOE_fulfilled_reweighting.ipynb",
    conda:
        "../envs/stats_comp.yaml"
    threads: 2
    notebook:
        "../notebooks/comp_methods_review_reweight_all.py.ipynb"
