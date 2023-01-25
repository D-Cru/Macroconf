# These rules produce AMBER18 topology files.
# TODO: add all output files to md_make_topology rule
ruleorder: md_pre_topology_exp > md_pre_topology_D > md_pre_topology


rule md_pre_topology:
    input:
        "{data_dir}/{compound_dir}/ob.pdb",
    output:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/4amb_noH.pdb",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    log:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/log.log",
    shell:
        "pdb4amber -y -i {input} -o {output} --logfile {log}"


rule md_pre_topology_D:
    input:
        "{data_dir}/{compound_dir}/ob_d.pdb",
    output:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/4amb_noH.pdb",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    log:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/log.log",
    shell:
        "pdb4amber -i {input} -o {output} --logfile {log}"


rule md_pre_topology_exp:
    input:
        "{data_dir}/{compound_dir}/ob_exp.pdb",
    output:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/4amb_noH.pdb",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    conda:
        "../envs/ambertools.yml",
    log:
        "{data_dir}/{compound_dir}/{solvent}/0_pre_topology/log.log",
    shell:
        "pdb4amber -i {input} -o {output} --logfile {log}"


rule md_build_leap:
    input:
        parm="{data_dir}/{compound_dir}/data.json",
        pdb="{data_dir}/{compound_dir}/{solvent}/0_pre_topology/4amb_noH.pdb",
        ff=config["forcefield"],
    output:
        tleap="{data_dir}/{compound_dir}/{solvent}/1_make_topology/tleap.in",
    params:
        path="{data_dir}/{compound_dir}/{solvent}/1_make_topology/",
    log:
        "{data_dir}/{compound_dir}/{solvent}/1_make_topology/leap.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/make-tleap.py"


rule md_make_topology:
    input:
        tleap="{dir}/tleap.in",
    output:
        top="{dir}/mc_sol.prmtop",
        coord="{dir}/mc_sol.inpcrd",
        pdb="{dir}/mc_sol.pdb",
        mol2="{dir}/mc_gas.mol2",
    log:
        "{dir}/leap.stdout",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "tleap -f {input.tleap} > {log}"


# Make DMSO, H2O cosolvation box. Does not fully work.
rule md_make_packmol:
    input:
        solute=f"data/interim/{config['exp_name']}/{{compound}}/H2O/1_make_topology/mc_gas.pdb", # Make sure that the resiude names in the .pdb files are different!
        solvent1="libs/md_solvents/cosolvation/h2o.pdb",
        solvent2="libs/md_solvents/cosolvation/dmso.pdb",
    output:
        packmol=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/packmol.inp",
    log:
        f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/packmol.log",
    params:
        system=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/system.pdb",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/make_packmol.py"


rule md_make_cosolvation:
    input:
        packmol=rules.md_make_packmol.output.packmol,
        solute=f"data/interim/{config['exp_name']}/{{compound}}/H2O/1_make_topology/mc_gas.pdb",
        solvent1="libs/md_solvents/cosolvation/h2o.pdb",
        solvent2="libs/md_solvents/cosolvation/dmso.pdb",
    output:
        system_solvated=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/system.pdb",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    log:
        f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/packmol.log",
    shell:
        "packmol < {input.packmol} > {log}"


rule md_make_tleap_cosolvation:
    input:
        system=rules.md_make_cosolvation.output.system_solvated,
        solute=f"data/interim/{config['exp_name']}/{{compound}}/H2O/1_make_topology/mc_gas.mol2",
        water="libs/md_solvents/cosolvation/water.pdb",
        dmso="libs/md_solvents/cosolvation/DMSO.mol2",
        dmso_mod="libs/md_solvents/cosolvation/dmso.frcmod",
        parm=f"data/interim/{config['exp_name']}/{{compound}}/data.json",
        pdb=f"data/interim/{config['exp_name']}/{{compound}}/H2O/0_pre_topology/4amb_noH.pdb",
        ff=config["forcefield"],
    output:
        tleap=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/tleap_cosol.in",
    params:
        path=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/",
        top=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.prmtop",
        coord=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.inpcrd",
        pdb=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.pdb",
        mol2=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_gas.mol2", #log:
    log:
        f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/leap.stdout",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    script:
        "../scripts/make-tleap-cosol.py"


rule md_run_tleap_cosolvation:
    input:
        tleap=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/tleap_cosol.in",
    output:
        top=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.prmtop",
        coord=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.inpcrd",
        pdb=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_sol.pdb",
        mol2=f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/mc_gas.mol2", #log:
    log:
        f"data/interim/{config['exp_name']}/{{compound}}/cosol_DMSO_H2O_/1_make_topology/leap.stdout",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "tleap -f {input.tleap} > {log}"
