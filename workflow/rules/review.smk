rule md_GaMD_conv_analysis:
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
            f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_pca_dihed.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        conv_plot=report(
            f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_conv_plot.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        grid_cells=report(
            f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_grid_cells.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        conv_data=report(
            f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_conv_data.json",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_stats=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_stats.json",
        noe_result=f"data/processed/review/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_result.json",
    params:
        cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/",
        rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
        method="GaMD",
        traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
        weights_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_weights_short.dat",
        dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
        dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
    log:
        notebook=f"data/processed/review/notebooks/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_GaMD_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD_convergence.py.ipynb"


def conv_accumulate(wildcards):
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
        ] = f"data/processed/review/results/{w.compound}/{w.solvent}/{w.method}/{w.simtime}/{w.repeats}/{sim_hash}_conv_data.json"

    return output


# Aggregate 1 method type (e.g. all cMD,2k ns with exactly the same parameters)
rule md_GaMD_conv_accumulate:
    input:
        unpack(conv_accumulate),
    output:
        plot=report(
            "data/processed/review/results/methods/{method}-{solvent}-{simtime}-{igamd}-conv_hist.svg"
        ),
    log:
        notebook="data/processed/review/notebooks/methods/{method}-{solvent}-{simtime}-{igamd}-GaMD_conv_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/review_convergence.py.ipynb"
