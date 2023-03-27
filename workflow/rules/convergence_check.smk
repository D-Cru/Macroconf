
rule conv_check_GaMD_make_param:
    input:
        template="libs/md_parameters/GaMD.tmpl",
    output:
        param="{dir}/interim/conv_check/{compound}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_GaMD_conv.in",
    log:
        "{dir}/interim/conv_check/{compound}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_GaMD_conv_make_param.log",
    params:
        sample=lambda w: samples.loc["{}".format(w.index)].to_dict(),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/GaMD_params.py"


rule conv_check_topol:
    input:
        f"{{dir}}/{config['exp_name']}/{{dir2}}/1_make_topology/mc_sol.prmtop",
    output:
        "{dir}/conv_check/{dir2}/1_make_topology/mc_sol_2.prmtop",
    log:
        "{dir}/conv_check/{dir2}/1_make_topology/mc_sol_2_check_topol.log",
    shell:
        "cp {input} {output}"


rule conv_check_GaMD_get_coord:
    input:
        f"data/processed/{config['exp_name']}/results/{{compound}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/rst_{{cluster_id}}.rst",
    output:
        "data/interim/conv_check/{compound}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/rst_{cluster_id}.rst",
    log:
        "data/interim/conv_check/{compound}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/rst_{cluster_id}_get_coord.log",
    shell:
        "cp {input} {output}"


rule conv_check_GaMD_full:
    input:
        param="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_GaMD_conv.in",
        top="{dir}/conv_check/{dir2}/1_make_topology/mc_sol_2.prmtop",
        coord="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/rst_{cluster_id}.rst",  # 0 hardcoded.. bad...
    output:
        out="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_md.out",
        restart="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_md-1.rst",
        mdinfo="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_md-1.info",
        traj="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_traj.netcdf",
        gamd="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_gamd.log",
        gamd_restart="{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_gamd-restart.dat",
    log:
        "{dir}/conv_check/{dir2}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_md-1.log",
    shadow:
        "full"
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    resources:
        gpu=1,
        runtime=lambda w: int(
            (int("{}".format(w.time)) + 52) / float(config["ns_h"]) * 60
        ),
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -gamd {output.gamd} 2> {log} && mv gamd-restart.dat {output.gamd_restart}"


rule conv_check_GaMD_pre_ana:
    input:
        gamd_log="{data_dir}/{compound_dir}/{md_dir}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_gamd.log",
    output:
        weights="{data_dir}/{compound_dir}/{md_dir}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_weights.dat",
    log:
        "{data_dir}/{compound_dir}/{md_dir}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_weights_pre_ana.log",
    params:
        lines=lambda w: int(
            int(samples.loc["{}".format(w.index)]["simtime"]) * 1000 / 2
        ),
    shell:
        """tail -n {params.lines} {input.gamd_log} | awk 'NR%1==0' | awk '{{print ($8+$7)/(0.001987*300)" " $2  " " ($8+$7)}}' > {output.weights}"""


rule conv_check_GaMD_anal:
    input:
        weights="data/interim/conv_check/{compound_dir}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_weights.dat",
        out="data/interim/conv_check/{compound_dir}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_md.out",
        traj="data/interim/conv_check/{compound_dir}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_traj.netcdf",
        traj_ncdf="data/interim/conv_check/{compound_dir}/{solvent}/11_GaMD_full/{time}/{repeat}/{index}/{cluster_id}_traj.ncdf",
        top="data/interim/conv_check/{compound_dir}/{solvent}/1_make_topology/mc_sol_2.prmtop",
        parm=expand(
            "data/interim/{exp_name}/{{compound_dir}}/data.json",
            exp_name=config["exp_name"],
        )[0],
        noe=expand(
            "data/interim/{exp_name}/{{compound_dir}}/NOE.json",
            exp_name=config["exp_name"],
        )[0],
    output:
        pca_dihe=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_pca_dihed.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_plot=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_noe.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_pmf=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_noe_pmf.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_plot=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_cluster_plot.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_pca=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_cluster_pca.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_min_samp=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_cluster_min_samp.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_time=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_cluster_time.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_structs=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_cluster_structs.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_stat_plot=report(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_noe_stat_plot.png",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_pdb="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_clusters/clusters.pdb",
        noe_result="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_noe_result.json",
        noe_stats="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_noe_stats.json",
        dihedrals="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_dihedrals.dat",
        dPCA="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_dPCA.dat",
        dPCA_weights_MC="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_dPCA_weights_MC.dat",
        noe_dist="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_NOE_dist.dat",
        multiple="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_multiple.dat",
        short_traj="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_short-traj.netcdf",
        dihedrals_short="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_dihedrals_short.dat",
        dPCA_weights_MC_short="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_dPCA_weights_MC_short.dat",
        cluster_restart=touch(
            "data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_clusters/rst/done.done"
        ),
        # change this to directory once fixed.
    params:
        cluster_dir="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_clusters/",
        rst_dir="data/processed/conv_check/results/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_clusters/rst/",
        method="GaMD",
    log:
        notebook="data/processed/conv_check/notebooks/{compound_dir}/{solvent}/GaMD/{time}/{repeat}/{index}/{cluster_id}_GaMD_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD.py.ipynb"


def aggregate_input(wildcards):
    """
    force evaluation of checkpoint. Then aggregate input files.
    https://stackoverflow.com/questions/62876986/snakemake-how-to-use-glob-wildcards-for-newly-created-files
    The requested file that is returned here can be produced via the md_GaMD_analysis rule!
    """
    checkpoints.md_GaMD_analysis.get(
        index=wildcards.index,
        time=wildcards.time,
        repeat=wildcards.repeat,
        compound_dir=wildcards.compound,
        solvent=wildcards.solvent,
        exp_name="conv_check",
    )
    ids = glob_wildcards(
        f"data/processed/{config['exp_name']}/results/{wildcards.compound}/{wildcards.solvent}/{wildcards.method}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}_clusters/rst/rst_{{id}}.rst"
    ).id

    if wildcards.method == "GaMD":
        method_interim = "11_GaMD_full"

    output = {}
    output["dihedrals"] = expand(
        f"data/processed/conv_check/results/{wildcards.compound}/{wildcards.solvent}/{wildcards.method}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}/{{id}}_dihedrals.dat",
        id=ids,
    )
    output["weights"] = expand(
        f"data/interim/conv_check/{wildcards.compound}/{wildcards.solvent}/{method_interim}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}/{{id}}_weights.dat",
        id=ids,
    )
    output["starting_struct"] = expand(
        f"data/interim/conv_check/{wildcards.compound}/{wildcards.solvent}/{method_interim}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}/rst_{{id}}.rst",
        id=ids,
    )

    output["top"] = (
        f"data/interim/conv_check/{wildcards.compound}/{wildcards.solvent}/1_make_topology/mc_sol_2.prmtop",
    )

    output["ref_pca"] = (
        f"data/processed/{config['exp_name']}/results/{wildcards.compound}/{wildcards.solvent}/{wildcards.method}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}_dPCA.dat",
    )
    output["ref_dih"] = (
        f"data/processed/{config['exp_name']}/results/{wildcards.compound}/{wildcards.solvent}/{wildcards.method}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}_dihedrals.dat",
    )
    output["ref_pca_weights"] = (
        f"data/processed/{config['exp_name']}/results/{wildcards.compound}/{wildcards.solvent}/{wildcards.method}/{wildcards.time}/{wildcards.repeat}/{wildcards.index}_dPCA_weights_MC.dat",
    )
    return output


rule analyse_conv_check:
    input:
        unpack(aggregate_input),
    output:
        dih_pca_comparison="data/processed/conv_check/results/{compound}/{solvent}/{method}/{time}/{repeat}/{index}_comparison-plot.svg",
    log:
        notebook="data/processed/conv_check/notebooks/{compound}/{solvent}_{method}_{time}_{repeat}_{index}_convcheck.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 2
    notebook:
        "../notebooks/conv_check.py.ipynb"
