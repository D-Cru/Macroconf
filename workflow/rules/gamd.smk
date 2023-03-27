# Rules to run GaMD simulations in Amber18


rule md_GaMD_make_param:
    input:
        template="libs/md_parameters/GaMD.tmpl",
    output:
        param="{dir}/11_GaMD_full/{time}/{repeat}/{index}_GaMD.in",
    log:
        "{dir}/11_GaMD_full/{time}/{repeat}/{index}_GaMD_makeparam.log",
    params:
        sample=lambda w: samples.loc["{}".format(w.index)].to_dict(),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/GaMD_params.py"


rule md_GaMD_full:
    input:
        param=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_GaMD.in",
        top=f"{{dir}}/{config['exp_name']}/{{dir2}}/1_make_topology/mc_sol.prmtop",
        coord=f"{{dir}}/{config['exp_name']}/{{dir2}}/7_equil_3/eq_3.rst",
    output:
        out=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md.out",
        restart=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md-1.rst",
        mdinfo=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md-1.info",
        traj=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_traj.netcdf",
        gamd=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_gamd.log",
        gamd_restart=f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_gamd-restart.dat",
    log:
        f"{{dir}}/{config['exp_name']}/{{dir2}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_md-1.log",
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


def compute_lines(wildcards):
    """compute no of weight lines  to keep"""
    output = {}
    if int(samples.loc["{}".format(wildcards.index)]["simtime"]) > 1000:
        lines = (
            int(
                int(samples.loc["{}".format(wildcards.index)]["simtime"])
                * 1000
                / 2
                / (
                    int(samples.loc["{}".format(wildcards.index)]["simtime"])
                    / 1000
                )
            ),
        )
    else:
        lines = (
            int(
                int(samples.loc["{}".format(wildcards.index)]["simtime"])
                * 1000
                / 2
            ),
        )
    output["lines"] = lines
    return lines


rule md_GaMD_pre_ana:
    input:
        gamd_log=f"{{data_dir}}/{config['exp_name']}/{{compound_dir}}/{{md_dir}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_gamd.log",
    output:
        weights=f"{{data_dir}}/{config['exp_name']}/{{compound_dir}}/{{md_dir}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_weights.dat",
    log:
        f"{{data_dir}}/{config['exp_name']}/{{compound_dir}}/{{md_dir}}/11_GaMD_full/{{time}}/{{repeat}}/{{index}}_pre_ana.log",
    params:
        lines=compute_lines,
    shell:
        """tail -n {params.lines} {input.gamd_log} | awk 'NR%1==0' | awk '{{print ($8+$7)/(0.001987*300)" " $2  " " ($8+$7)}}' > {output.weights}"""


checkpoint md_GaMD_analysis:
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
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_pca_dihed.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_pmf=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_pmf.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_cluster_plot.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_pca=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_cluster_pca.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_min_samp=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_cluster_min_samp.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_time=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_cluster_time.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_structs=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_cluster_structs.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        noe_stat_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_stat_plot.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        cluster_pdb=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/clusters.pdb",
        cluster_solvated=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/clusters_solvated.pdb",
        noe_result=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_result.json",
        noe_stats=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_noe_stats.json",
        dihedrals=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dihedrals.dat",
        dPCA=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA.dat",
        dPCA_weights_MC=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC.dat",
        noe_dist=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_NOE_dist.dat",
        multiple=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_multiple.dat",
        cluster_restart=touch(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/done.done"
        ),
        # change this to directory once fixed in smk.
        NPR_shape_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_NPR_shape.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        NPR_shape_data=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_NPR_shape.dat",
        NPR_shape_weights=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_NPR_shape_weights.dat",
        overview_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_overview.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        rmsd_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_rmsd_plot.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
        omega_plot=report(
            f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_omega_plot.svg",
            category="Compound {compound_dir}",
            subcategory="GaMD",
        ),
    params:
        cluster_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/",
        rst_dir=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_clusters/rst/",
        method="GaMD",
        traj_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_traj_short.netcdf",
        weights_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_weights_short.dat",
        dihedrals_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dihedrals_short.dat",
        dPCA_weights_MC_short=f"data/processed/{config['exp_name']}/results/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_dPCA_weights_MC_short.dat",
    log:
        notebook=f"data/processed/{config['exp_name']}/notebooks/{{compound_dir}}/{{solvent}}/GaMD/{{time}}/{{repeat}}/{{index}}_GaMD_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD.py.ipynb"
