# Rules to run energy minimisation & equilibration in Amber
rule md_em_1:
    """Run energy minimisation 1
    Energy minimisation 1 minimises the solvent only, restraining everything else.
    See the em_1.in file for more detailed parameters.
    """
    input:
        param="libs/md_parameters/{solvent}/em_1.in",
        top="{dir}/{solvent}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/{solvent}/1_make_topology/mc_sol.inpcrd",
    output:
        out="{dir}/{solvent}/2_minim_1/em_1.out",
        restart="{dir}/{solvent}/2_minim_1/em_1.rst",
        mdinfo="{dir}/{solvent}/2_minim_1/mdinfo.info",
    log:
        "{dir}/{solvent}/2_minim_1/em_1.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd -O -i {input.param} -o {output.out} -inf {output.mdinfo} -p {input.top} -c {input.coord} -r {output.restart} -ref {input.coord} 2> {log}"


rule md_em_2:
    """Run energy minimisation 2
    Energy minimisation 2 lets the solvent move, restraining everything else.
    See the em_2.in file for more detailed parameters.
    """
    input:
        param="libs/md_parameters/{solvent}/em_2.in",
        top="{dir}/{solvent}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/{solvent}/2_minim_1/em_1.rst",
    output:
        out="{dir}/{solvent}/3_minim_2/em_2.out",
        restart="{dir}/{solvent}/3_minim_2/em_2.rst",
        mdinfo="{dir}/{solvent}/3_minim_2/mdinfo.info",
        traj="{dir}/{solvent}/3_minim_2/traj.ncdf",
    log:
        "{dir}/{solvent}/3_minim_2/em_2.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -ref {input.coord} 2> {log}"


rule md_em_3:
    """Run energy minimisation 3
    Energy minimisation 3 minimises the whole system.
    See the em_3.in file for more detailed parameters.
    """
    input:
        param="libs/md_parameters/em_3.in",
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/3_minim_2/em_2.rst",
    output:
        out="{dir}/4_minim_3/em_3.out",
        restart="{dir}/4_minim_3/em_3.rst",
        mdinfo="{dir}/4_minim_3/mdinfo.info",
    log:
        "{dir}/4_minim_3/em_3.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd -O -i {input.param} -o {output.out} -inf {output.mdinfo} -p {input.top} -c {input.coord} -r {output.restart} -ref {input.coord} 2> {log}"


rule md_em_4:
    """Run energy minimisation 4
    Energy minimisation 4 heats the system in a NVT 0.5ps simulation to 300K.
    Heavy atoms (non-solvent) are restrained.
    See the em_4.in file for more detailed parameters.
    """
    input:
        param="libs/md_parameters/{solvent}/em_4.in",
        top="{dir}/{solvent}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/{solvent}/4_minim_3/em_3.rst",
    output:
        out="{dir}/{solvent}/5_equil_1/eq_1.out",
        restart="{dir}/{solvent}/5_equil_1/eq_1.rst",
        mdinfo="{dir}/{solvent}/5_equil_1/mdinfo.info",
        traj="{dir}/{solvent}/5_equil_1/traj.ncdf",
    log:
        "{dir}/{solvent}/5_equil_1/eq_1.log",
    resources:
        gpu=1,
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -ref {input.coord} 2> {log}"


rule md_eq_1:
    """Run equilibration 1
    Equilibration 1 equilibrates the system in a 0.5ns NPT simulation.
    Heavy atoms (non-solvent) are restrained.
    See the em_5.in file for more detailed parameters.
    """
    input:
        param="libs/md_parameters/{solvent}/em_5.in",
        top="{dir}/{solvent}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/{solvent}/5_equil_1/eq_1.rst",
    output:
        out="{dir}/{solvent}/6_equil_2/eq_2.out",
        restart="{dir}/{solvent}/6_equil_2/eq_2.rst",
        mdinfo="{dir}/{solvent}/6_equil_2/mdinfo.info",
        traj="{dir}/{solvent}/6_equil_2/traj.ncdf",
    log:
        "{dir}/{solvent}/6_equil_2/eq_2.log",
    resources:
        gpu=1,
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} -ref {input.coord} 2> {log}"


rule md_eq_2:
    """Run equilibration 2
    Equilibration 2 equilibrates the system in a 5ns NPT simulation.
    There are no restraints.
    See the eq.in file for more detailed parameters.
    """
    input:
        param=report("libs/md_parameters/eq.in", category="md-params"),
        top="{dir}/1_make_topology/mc_sol.prmtop",
        coord="{dir}/6_equil_2/eq_2.rst",
    output:
        out="{dir}/7_equil_3/eq_3.out",
        restart="{dir}/7_equil_3/eq_3.rst",
        mdinfo="{dir}/7_equil_3/mdinfo.info",
        traj="{dir}/7_equil_3/traj.ncdf",
    log:
        "{dir}/7_equil_3/eq_3.log",
    resources:
        gpu=1,
        runtime=30,
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "pmemd.cuda -O -i {input.param} -o {output.out} -inf {output.mdinfo} -x {output.traj} -p {input.top} -c {input.coord} -r {output.restart} 2> {log}"
