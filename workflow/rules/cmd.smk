# Rules to run conventional MD in Amber18
rule md_cMD_make_param:
    """Make parameter file for cMD simulation
    based on the template file. Inserts the parameters for the cMD simulation 
    that were specified in the samples.tsv input file. NPT ensemble is used.
    """
    input:
        template="libs/md_parameters/cMD.tmpl",
    output:
        param="{dir}/8_cMD/{time}/{repeat}/{index}_cMD.in",
    log:
        "{dir}/8_cMD/{time}/{repeat}/{index}_cMD_makeparam.log",
    params:
        sample=lambda w: samples.loc["{}".format(w.index)].to_dict(),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/cMD_params.py"


rule md_cMD_prod:
    """Run a production cMD simulation
    """
    input:
        param="{dir}/8_cMD/{time}/{repeat}/{index}_cMD.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/7_equil_3/eq_3.rst",
    output:
        out="{dir}/8_cMD/{time}/{repeat}/{index}_md.out",
        restart="{dir}/8_cMD/{time}/{repeat}/{index}_cMD.rst",
        mdinfo="{dir}/8_cMD/{time}/{repeat}/{index}_mdinfo.info",
        traj="{dir}/8_cMD/{time}/{repeat}/{index}_traj.netcdf",
    log:
        "{dir}/8_cMD/{time}/{repeat}/{index}_cMD.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    resources:
        gpu=1,
        runtime=lambda w: int(
            float("{}".format(w.time)) / float(config["ns_h"]) * 60
        ),
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} 2> {log}"


rule md_cMD_analysis:
    input:
        out="data/interim/{exp_name}/{compound_dir}/{solvent}/8_cMD/{time}/{repeat}/{index}_md.out",
        traj="data/interim/{exp_name}/{compound_dir}/{solvent}/8_cMD/{time}/{repeat}/{index}_traj.netcdf",
        traj_ncdf="data/interim/{exp_name}/{compound_dir}/{solvent}/8_cMD/{time}/{repeat}/{index}_traj.ncdf",
        top="data/interim/{exp_name}/{compound_dir}/{solvent}/1_make_topology/mc_sol.prmtop",
        parm="data/interim/{exp_name}/{compound_dir}/data.json",
        noe="data/interim/{exp_name}/{compound_dir}/NOE.json",
        ref_mol="data/interim/{exp_name}/{compound_dir}/{solvent}/1_make_topology/mc_gas.mol2",
    output:
        pca_dihe=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_pca_dihed.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        noe_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_noe.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        noe_pmf=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_noe_pmf.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cluster_plot.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_pca=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cluster_pca.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_min_samp=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cluster_min_samp.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_time=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cluster_time.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_structs=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cluster_structs.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        noe_stat_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_noe_stat_plot.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        cluster_pdb="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_clusters/clusters.pdb",
        cluster_solvated="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_clusters/clusters_solvated.pdb",
        noe_result="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_noe_result.json",
        noe_stats="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_noe_stats.json",
        dihedrals="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_dihedrals.dat",
        dPCA="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_dPCA.dat",
        dPCA_weights_MC="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_dPCA_weights_MC.dat",
        noe_dist="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_NOE_dist.dat",
        multiple="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_multiple.dat",
        cluster_restart=touch(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_clusters/rst/done.done"
        ), # change this to directory once fixed in Snakemake.
        NPR_shape_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_NPR_shape.png",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
        NPR_shape_data="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_NPR_shape.dat",
        NPR_shape_weights="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_NPR_shape_weights.dat",
        overview_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_overview.svg",
            category="Compound {compound_dir}",
            subcategory="cMD",
        ),
    params:
        cluster_dir="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_clusters/",
        rst_dir="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_clusters/rst/",
        method="cMD",
        traj_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_traj_short.netcdf",
        dihedrals_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_dihedrals_short.dat",
        dPCA_weights_MC_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_dPCA_weights_MC_short.dat",
    log:
        notebook="data/processed/{exp_name}/notebooks/{compound_dir}/{solvent}/cMD/{time}/{repeat}/{index}_cMD_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD.py.ipynb"
