# Old gromacs rules.
# index = 22
# compound_dir = 'data/interim/22-02-2021_MacroConf-v2/{index}/'
# md_dir = "{data_dir}/{compound_dir}md-auto-test/"
# md_dir

# workdir: "./data/interim/22-02-2021_MacroConf-v2/{index}"

# IDS = "1 2 3 ...".split() # the list of desired ids

# a pseudo-rule that collects the target files
# rule all:
# input:  expand("otherdir/{id}.bam", id=IDS)

threads_max = 24


rule premd:
    input:
        [
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-4/md/md_fit.xtc",
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-4/eq/plots/density.png",
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-4/em/plots/em_potential_energy.png",
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-4/md/analysis.xvg",
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-4/results/pca_cart.png",
        ],


rule test2:
    input:
        ["data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-3/eq/plots/density.png"],


rule all:
    input:
        [
            "data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-2/em/plots/em_potential_energy.png"
        ],


rule test:
    input:
        ["data/interim/22-02-2021_MacroConf-v2/22/md-auto-test-2/eq/analysis.xvg"],


rule md_make_topology:
    input:
        "{data_dir}/{compound_dir}/gmx.pdb",
    output:
        gro="{data_dir}/{compound_dir}/{md_dir}/macrocycle.gro",
        top="{data_dir}/{compound_dir}/{md_dir}/topol.0.top",
        res="{data_dir}/{compound_dir}/{md_dir}/posre.itp", #log:
         #"{data_dir}/{compound_dir}/{md_dir}/logs/{rule}.log",
    conda:
        "../envs/gmx2021.yaml"
    shell:
        "gmx pdb2gmx -f {input} -o {output.gro} -p {output.top} -ff amber99sb-ildn -water tip3p -i {output.res}"


rule md_add_box:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/macrocycle.gro",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/boxed.gro",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "gmx editconf -f {input} -o {output} -c -d 1.0 -bt cubic"


rule md_solvate:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/boxed.gro",
        "{data_dir}/{compound_dir}/{md_dir}/topol.0.top",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/solvated.gro",
        "{data_dir}/{compound_dir}/{md_dir}/topol.1.top",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input[1]} {output[1]} &&"
        "gmx solvate -cp {input[0]} -cs -o {output[0]} -p {output[1]}"


rule md_add_ion1:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/solvated.gro",
        "{data_dir}/{compound_dir}/{md_dir}/topol.1.top",
        "{data_dir}/data/genion.mdp",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/genion.tpr",
        "{data_dir}/{compound_dir}/{md_dir}/topol.2.top",
        "{data_dir}/{compound_dir}/{md_dir}/mdout.mdp",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input[1]} {output[1]} &&"
        "gmx grompp -c {input[0]} -p {output[1]} -f {input[2]} -o {output[0]} -po {output[2]} -v"


rule md_add_ion2:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/genion.tpr",
        "{data_dir}/{compound_dir}/{md_dir}/topol.2.top",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/system.gro",
        "{data_dir}/{compound_dir}/{md_dir}/topol.3.top",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input[1]} {output[1]} &&"
        "printf '13' | gmx genion -s {input[0]} -conc 0.0 -neutral -pname NA -nname CL -o {output[0]} -p {output[1]}"


rule md_em_prep:
    input:
        gro="{data_dir}/{compound_dir}/{md_dir}/system.gro",
        top="{data_dir}/{compound_dir}/{md_dir}/topol.3.top",
        para="{data_dir}/data/em.mdp",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/em/em.tpr",
        "{data_dir}/{compound_dir}/{md_dir}/topol.4.top",
        "{data_dir}/{compound_dir}/{md_dir}/em/mdout.mdp",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input.top} {output[1]} &&"
        "gmx grompp -c {input.gro} -p {output[1]} -f {input.para} -o {output[0]} -po {output[2]}"


rule md_em_run:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/em/em.tpr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/em/em.edr",
        "{data_dir}/{compound_dir}/{md_dir}/em/em.trr",
        "{data_dir}/{compound_dir}/{md_dir}/em/em.gro",
    params:
        wd="{data_dir}/{compound_dir}/{md_dir}/em/",
    threads: threads_max
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "gmx mdrun -s {input} -deffnm {params.wd}em -v -update cpu -ntmpi 1"


rule md_em_anal:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/em/em.tpr",
        "{data_dir}/{compound_dir}/{md_dir}/em/em.edr",
        "{data_dir}/{compound_dir}/{md_dir}/em/em.gro",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/em/plots/em_potential_energy.xvg",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf 'Potential' | gmx energy -s {input[0]} -f {input[1]} -o {output}"


rule md_em_anal_plot:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/em/plots/em_potential_energy.xvg",
    output:
        report(
            "{data_dir}/{compound_dir}/{md_dir}/em/plots/em_potential_energy.png",
            category="Energy Minimization",
        ),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_pot_energy.py"


rule md_equi_prep:
    input:
        gro="{data_dir}/{compound_dir}/{md_dir}/em/em.gro",
        top="{data_dir}/{compound_dir}/{md_dir}/topol.4.top",
        para="{data_dir}/data/eq.mdp",
    output:
        top="{data_dir}/{compound_dir}/{md_dir}/topol.5.top",
        tpr="{data_dir}/{compound_dir}/{md_dir}/eq/equi.tpr",
        mdp="{data_dir}/{compound_dir}/{md_dir}/eq/mdout.mdp",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input.top} {output.top} &&"
        "gmx grompp -c {input.gro} -p {output.top} -r {input.gro} -f {input.para} -o {output.tpr} -po {output.mdp}"


rule md_equi_run:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/eq/equi.tpr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/eq/equi.edr",
        "{data_dir}/{compound_dir}/{md_dir}/eq/equi.xtc",
        "{data_dir}/{compound_dir}/{md_dir}/eq/equi.gro",
    params:
        wd="{data_dir}/{compound_dir}/{md_dir}/eq/",
    threads: threads_max
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "gmx mdrun -s {input} -deffnm {params.wd}equi -v -update gpu"


rule md_equi_anal:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/eq/equi.edr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/eq/analysis.xvg",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf 'Potential\nTemperature\n Density' | gmx energy -f {input} -o {output}"


rule md_equi_anal_plot:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/eq/analysis.xvg",
    output:
        pot=report(
            "{data_dir}/{compound_dir}/{md_dir}/eq/plots/potential.png",
            category="Equilibration",
        ),
        temp=report(
            "{data_dir}/{compound_dir}/{md_dir}/eq/plots/temperature.png",
            category="Equilibration",
        ),
        dens=report(
            "{data_dir}/{compound_dir}/{md_dir}/eq/plots/density.png",
            category="Equilibration",
        ),
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_equi_stats.py"


rule md_md_setup:
    input:
        gro="{data_dir}/{compound_dir}/{md_dir}/eq/equi.gro",
        top="{data_dir}/{compound_dir}/{md_dir}/topol.5.top",
        para="{data_dir}/data/md.mdp",
    output:
        top="{data_dir}/{compound_dir}/{md_dir}/topol.6.top",
        tpr="{data_dir}/{compound_dir}/{md_dir}/md/md.tpr",
        mdp="{data_dir}/{compound_dir}/{md_dir}/md/mdout.mdp",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "cp {input.top} {output.top} && "
        "gmx grompp -c {input.gro} -p {output.top} -f {input.para} -o {output.tpr} -po {output.mdp}" #"gmx grompp -c equi.gro -p ../topol.top -f ../../../data/md.mdp -o md.tpr"


rule md_md_run:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/md/md.tpr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/md/md.edr",
        protected("{data_dir}/{compound_dir}/{md_dir}/md/md.xtc"),
        "{data_dir}/{compound_dir}/{md_dir}/md/md.gro",
    params:
        wd="{data_dir}/{compound_dir}/{md_dir}/md/",
    threads: threads_max
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "gmx mdrun -s {input} -deffnm {params.wd}md -v -update gpu" #gmx mdrun -deffnm md -v -update gpu


rule md_md_conv_rm_jumps:
    input:
        traj="{data_dir}/{compound_dir}/{md_dir}/md/md.xtc",
        tpr="{data_dir}/{compound_dir}/{md_dir}/md/md.tpr",
    output:
        temp("{data_dir}/{compound_dir}/{md_dir}/md/md_nojump.xtc"),
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf '0' | gmx trjconv -f {input.traj} -s {input.tpr} -pbc nojump -o {output}"


# first remove any jumps over the box boundaries
#!gmx trjconv -f md.xtc -s md.tpr -pbc nojump -o md_nojump.xtc # type 0


rule md_md_conv_set_com:
    input:
        traj="{data_dir}/{compound_dir}/{md_dir}/md/md_nojump.xtc",
        tpr="{data_dir}/{compound_dir}/{md_dir}/md/md.tpr",
    output:
        temp("{data_dir}/{compound_dir}/{md_dir}/md/md_center.xtc"),
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf '1\n 0' | gmx trjconv -f {input.traj} -s {input.tpr} -pbc mol -ur compact -center -o {output}" # set center of mass to box center
         #! gmx trjconv -f md_nojump.xtc -s md.tpr -pbc mol -ur compact -center -o md_center.xtc # type 1 then 0


rule md_md_conv_rm_tr:
    input:
        traj="{data_dir}/{compound_dir}/{md_dir}/md/md_center.xtc",
        tpr="{data_dir}/{compound_dir}/{md_dir}/md/md.tpr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/md/md_fit.xtc",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf '1\n 0' | gmx trjconv -f {input.traj} -s {input.tpr} -fit rot+trans -o {output}" # Remove translational/rotational motions
         #!gmx trjconv -f md_center.xtc -s md.tpr -fit rot+trans -o md_fit.xtc # type 1 then 0
         # Remove intermediates
         # rm md_nojump.xtc md_center.xtc


rule md_md_anal:
    input:
        "{data_dir}/{compound_dir}/{md_dir}/md/md.edr",
    output:
        "{data_dir}/{compound_dir}/{md_dir}/md/analysis.xvg",
    envmodules:
        "gromacs/2020.3-AVX512-GPU",
    shell:
        "printf '7\n 9\n 11' | gmx energy -f {input} -o {output}" # Get temperature
         #!gmx energy -f md.edr -s md.tpr -o 1hsg_temperature.xvg  # type 15 then 0
         # Get energies
         #! gmx energy -s md.tpr -f md.edr -o 1hsg_energies.xvg # type 7 9 11 0


rule md_md_anal_mda:
    input:
        traj="{data_dir}/{compound_dir}/{md_dir}/md/md_fit.xtc",
        gro="{data_dir}/{compound_dir}/{md_dir}/eq/equi.gro",
        ener="{data_dir}/{compound_dir}/{md_dir}/md/analysis.xvg",
    output:
        noe_plot=report(
            "{data_dir}/{compound_dir}/{md_dir}/results/noe.png", category="MD"
        ),
        pca_cart=report(
            "{data_dir}/{compound_dir}/{md_dir}/results/pca_cart.png", category="MD"
        ),
        pca_dihe=report(
            "{data_dir}/{compound_dir}/{md_dir}/results/pca_dihed.png", category="MD"
        ), #"NOE"
    log:
        notebook="{data_dir}/{compound_dir}/{md_dir}/notebooks/processed.ipynb",
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/template-analyse-MD.py.ipynb"
