# Rules to run aMD in Amber18
rule md_aMD_cMD:
    input:
        param="libs/md_parameters/aMD_cMD.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/7_equil_3/eq_3.rst",
    output:
        out="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_md.out",
        restart="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_cMD.rst",
        mdinfo="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_mdinfo.info",
        traj="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_traj.netcdf",
    log:
        "{dir}/9_aMD_cMD/{time}/{repeat}/{index}_cMD.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    resources:
        gpu=1,
        runtime=60,
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} 2> {log}"


rule md_aMD_cMD_readE:
    input:
        out="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_md.out",
    output:
        cpptraj="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_cpptraj.in",
        epot="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_Epot.dat",
        edih="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_Edih.dat",
    log:
        "{dir}/9_aMD_cMD/{time}/{repeat}/{index}_cpptraj.log",
    conda:
        "../envs/ambertools.yml"
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        'printf "readdata {input.out} name mdout \nwritedata {output.epot} mdout[EPtot] time 0.0002 \nwritedata {output.edih} mdout[DIHED] time 0.0002 \n" > {output.cpptraj} & cpptraj -i {output.cpptraj} -o {log}'


rule md_aMD_make_param:
    input:
        template="libs/md_parameters/aMD.tmpl",
        traj="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_traj.netcdf",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        epot="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_Epot.dat",
        edih="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_Edih.dat",
    output:
        param="{dir}/10_aMD_prod/{time}/{repeat}/{index}_aMD.in",
    params:
        sample=lambda w: samples.loc["{}".format(w.index)].to_dict(),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/aMD_params.py"


rule md_aMD_prod:
    input:
        param="{dir}/10_aMD_prod/{time}/{repeat}/{index}_aMD.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/9_aMD_cMD/{time}/{repeat}/{index}_cMD.rst",
    output:
        out="{dir}/10_aMD_prod/{time}/{repeat}/{index}_md.out",
        restart="{dir}/10_aMD_prod/{time}/{repeat}/{index}_aMD.rst",
        mdinfo="{dir}/10_aMD_prod/{time}/{repeat}/{index}_mdinfo.info",
        traj="{dir}/10_aMD_prod/{time}/{repeat}/{index}_traj.netcdf",
        amd_log="{dir}/10_aMD_prod/{time}/{repeat}/{index}_aMD.log",
    log:
        "{dir}/10_aMD_prod/{time}/{repeat}/{index}_cMD.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    resources:
        gpu=1,
        runtime=lambda w: int(
            float("{}".format(w.time)) / float(config["ns_h"]) * 60
        ),
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -amd {output.amd_log} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} 2> {log}"


rule md_aMD_pre_ana:
    input:
        amd_log="{data_dir}/{compound_dir}/{md_dir}/10_aMD_prod/{time}/{repeat}/{index}_aMD.log",
    output:
        weights="{data_dir}/{compound_dir}/{md_dir}/10_aMD_prod/{time}/{repeat}/{index}_weights.dat",
    log:
        "{data_dir}/{compound_dir}/{md_dir}/10_aMD_prod/{time}/{repeat}/{index}_pre_ana.log",
    shell:
        """awk 'NR%1==0' {input.amd_log} | awk '{{print ($8+$7)/(0.001987*300)" " $2  " " ($8+$7)}}' > {output.weights} && sed -i -e 1,3d {output.weights}"""


rule md_aMD_analysis:
    input:
        weights="data/interim/{exp_name}/{compound_dir}/{solvent}/10_aMD_prod/{time}/{repeat}/{index}_weights.dat",
        out="data/interim/{exp_name}/{compound_dir}/{solvent}/10_aMD_prod/{time}/{repeat}/{index}_md.out",
        traj="data/interim/{exp_name}/{compound_dir}/{solvent}/10_aMD_prod/{time}/{repeat}/{index}_traj.netcdf",
        traj_ncdf="data/interim/{exp_name}/{compound_dir}/{solvent}/10_aMD_prod/{time}/{repeat}/{index}_traj.ncdf",
        top="data/interim/{exp_name}/{compound_dir}/{solvent}/1_make_topology/mc_sol.prmtop",  #noe_csv="{data_dir}/{compound_dir}.csv",
        parm="data/interim/{exp_name}/{compound_dir}/data.json",
        noe="data/interim/{exp_name}/{compound_dir}/NOE.json",
        ref_mol="data/interim/{exp_name}/{compound_dir}/{solvent}/1_make_topology/mc_gas.mol2",
    output:
        pca_dihe=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_pca_dihed.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        noe_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_noe.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        noe_pmf=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_noe_pmf.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        cluster_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_cluster_plot.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        cluster_pca=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_cluster_pca.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        cluster_min_samp=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_cluster_min_samp.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        cluster_time=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_cluster_time.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        cluster_structs=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_cluster_structs.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        noe_stat_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_noe_stat_plot.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        conv_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_conv_plot.svg",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        grid_cells=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_grid_cells.svg",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        conv_data="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_conv_data.json",
        cluster_pdb="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_clusters/clusters.pdb",
        cluster_solvated="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_clusters/clusters_solvated.pdb",
        noe_result="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_noe_result.json",
        noe_stats="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_noe_stats.json",
        dihedrals="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_dihedrals.dat",
        dPCA="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_dPCA.dat",
        dPCA_weights_MC="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_dPCA_weights_MC.dat",
        noe_dist="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_NOE_dist.dat",
        multiple="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_multiple.dat",
        cluster_restart=touch(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_clusters/rst/done.done"
        ),
        # change this to directory once fixed in SMK.
        NPR_shape_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_NPR_shape.png",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        NPR_shape_data="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_NPR_shape.dat",
        NPR_shape_weights="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_NPR_shape_weights.dat",
        overview_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_overview.svg",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        rmsd_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_rmsd_plot.svg",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        omega_plot=report(
            "data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_omega_plot.svg",
            category="Compound {compound_dir}",
            subcategory="aMD",
        ),
        sasa=f"data/processed/{{exp_name}}/results/{{compound_dir}}/{{solvent}}/aMD/{{time}}/{{repeat}}/{{index}}_sasa.pkl",
        psa=f"data/processed/{{exp_name}}/results/{{compound_dir}}/{{solvent}}/aMD/{{time}}/{{repeat}}/{{index}}_psa.pkl",
        solvation_properties=f"data/processed/{{exp_name}}/results/{{compound_dir}}/{{solvent}}/aMD/{{time}}/{{repeat}}/{{index}}_solvation_properties.json",
    params:
        cluster_dir="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_clusters/",
        rst_dir="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_clusters/rst/",
        method="aMD",
        traj_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_traj_short.netcdf",
        weights_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_weights_short.dat",
        dihedrals_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_dihedrals_short.dat",
        dPCA_weights_MC_short="data/processed/{exp_name}/results/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_dPCA_weights_MC_short.dat",
    log:
        notebook="data/processed/{exp_name}/notebooks/{compound_dir}/{solvent}/aMD/{time}/{repeat}/{index}_aMD_processed.ipynb",
    conda:
        "../envs/stats.yaml"
    threads: 16
    notebook:
        "../notebooks/template-analyse-amber-aMD.py.ipynb"
