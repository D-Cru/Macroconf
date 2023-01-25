# A collection of currently unused rules. Either outdated or not in use rn.
rule md_GaMD_equil:
    input:
        param="libs/md_parameters/GaMD_eq.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/7_equil_3/eq_3.rst",
    output:
        out="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_md.out",
        restart="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_md-1.rst",
        mdinfo="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_md-1.info",
        traj="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_traj.netcdf",
        gamd="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_gamd.log",
    log:
        "{dir}/11_GaMD_eq/{time}/{repeat}/{index}_md-1.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    resources:
        gpu=1,
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -gamd {output.gamd} 2> {log}"


rule md_GaMD_prod:
    input:
        param="libs/md_parameters/GaMD_prod.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/11_GaMD_eq/{time}/{repeat}/{index}_md-1.rst",
    output:
        out="{dir}/12_GaMD_eq/{time}/{repeat}/{index}_md-2.out",
        restart="{dir}/12_GaMD_eq/{time}/{repeat}/{index}_md-2.rst",
        mdinfo="{dir}/12_GaMD_eq/{time}/{repeat}/{index}_md-2.info",
        traj="{dir}/12_GaMD_eq/{time}/{repeat}/{index}_traj.netcdf",
        gamd="{dir}/12_GaMD_eq/{time}/{repeat}/{index}_gamd.log",
    log:
        "{dir}/12_GaMD_eq/{time}/{repeat}/{index}_md-2.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -gamd {output.gamd} 2> {log}"


rule md_final_anal:
    input:
        NOE_input=expand(
            "data/interim/22-02-2021_MacroConf-v2/{sample}/{exp_name}/results/{md_method}/noe_result.json",
            sample=config["compounds"],
            md_method=config["md_methods"],
            exp_name=config["exp_name"],
        ),
    output:
        NOE_plot=report(
            expand(
                "data/processed/{exp_name}/noe_plot.png",
                exp_name=config["exp_name"],
            )
        ),
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/final_analysis.py.ipynb"


rule send_mail:
    input:
        file="data/interim/22-02-2021_MacroConf-v2/24/md-auto-amber-2/11_GaMD_full/gamd.log",
    output:
        touch("sent.mail"),
    script:
        "../scripts/send_mail.py"


# configfile: "snakemake-config.yaml"
# # import main workflow
# include: "Snakefile"

# workflow
#
#
# rule compounds_to_make:
#     input:
#         "data/interim/22-02-2021_MacroConf-v2/24/ob.pdb", #"data/interim/22-02-2021_MacroConf-v2/22/ob.pdb",
#          #"data/interim/22-02-2021_MacroConf-v2/56/ob.pdb"
#         "data/interim/22-02-2021_MacroConf-v2/24/md-auto-amber-2/1_make_topology/mc_gas.png",


# rule get_ref_mol_smiles:
#     input:
#         ref_mol="{data_dir}/{compound_dir}/1_make_topology/mc_gas.mol2",
#     output:
#         smiles="{data_dir}/{compound_dir}/1_make_topology/ref_smiles.smi",
#     conda:
#         "../envs/mol-maker.yaml"
#     script:
#         "../scripts/make_smiles_from_mol.py"


rule omega_confgen:
    input:
        smile="{data_dir}/{compound_dir}/smile.smi",
    output:
        omega_bin="{data_dir}/{compound_dir}/{md_dir}/omega/mcs.oeb.gz",
    params:
        omega_param="libs/omega/oeomega_config.param",
        prefix="{data_dir}/{compound_dir}/{md_dir}/omega/",
    shell:
        "oeomega macrocycle -in {input.smile} -out {output.omega_bin} -param {params.omega_param} -prefix {params.prefix}" #  -mpi_np 24


rule omega_extract_mol:
    input:
        omega_bin="{data_dir}/{compound_dir}/{md_dir}/omega/mcs.oeb.gz",
    output:
        mol2="{data_dir}/{compound_dir}/{md_dir}/omega/mcs.mol2",
    conda:
        "../envs/oechem.yaml"
    script:
        "../scripts/oechem_bin_to_mol2.py"


# rule md_comp_analysis:
#     input:
#         unpack(compare_inputs),
#     output:
#         pca_dihe=report(
#             "data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}_{confgen}_{mode}-pca_dihe.png"
#         ),
#         cluster_pca=report(
#             "data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}_{confgen}_{mode}-cluster-pca.png"
#         ), # wildcard_constraints:
#          #     index_0="\[0-9A-Fa-f]",
#          #     index_1="\[0-9A-Fa-f]",
#          #     index_2="\[0-9A-Fa-f]",
#         report_pca_comparison="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}_{confgen}_{mode}-pca_dihe_report.svg",
#         shape_comparsion="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}_{confgen}_{mode}-shape_comparison.png",
#         cluster_hbonds="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}_{confgen}_{mode}-cluster_hbonds.png",
#     params:
#         sample_0=lambda w: samples.loc["{}".format(w.index_0)].to_dict(),
#         sample_1=lambda w: samples.loc["{}".format(w.index_1)].to_dict(),
#         sample_2=lambda w: samples.loc["{}".format(w.index_2)].to_dict(),
#     log:
#         notebook="data/processed/{exp_name}/notebooks/{compound_dir}/{index_0}_{index_1}_{index_2}_{confgen}_{mode}_compar.ipynb",
#     conda:
#         "../envs/stats.yaml"
#     threads: 2
#     notebook:
#         "../notebooks/template-analyse-amber-compar.py.ipynb"
