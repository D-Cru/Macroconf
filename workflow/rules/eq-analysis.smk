# Rules to analyse equilibration/energy minimisation in Amber18
rule md_readEpot:
    """
    Read the potential energy from a MD.out file
    """
    input:
        out="{somepath}/{file}.out",
    output:
        cpptraj=temp("{somepath}/{file}_epot_cpptraj.in"),
        epot="{somepath}/{file}_Epot.dat",
    log:
        "{somepath}/{file}_Epot.log",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        'printf "readdata {input.out} name mdout \nwritedata {output.epot} mdout[EPtot] time 0.002\n" > {output.cpptraj} & cpptraj -i {output.cpptraj} -o {log}'


# dihed: \nwritedata {output.edih} mdout[DIHED] time 0.0002
rule md_readT:
    """
    Read the temperature from a MD.out file
    """
    input:
        out="{somepath}/{file}.out",
    output:
        cpptraj=temp("{somepath}/{file}_t_cpptraj.in"),
        temp="{somepath}/{file}_temp.dat",
    log:
        "{somepath}/{file}_temp.log",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        'printf "readdata {input.out} name mdout \nwritedata {output.temp} mdout[TEMP] time 0.002\n" > {output.cpptraj} & cpptraj -i {output.cpptraj} -o {log}'


rule md_readrho:
    """
    Read the density from a MD.out file
    """
    input:
        out="{somepath}/{file}.out",
    output:
        cpptraj=temp("{somepath}/{file}_rho_cpptraj.in"),
        rho="{somepath}/{file}_rho.dat",
    log:
        "{somepath}/{file}_rho.log",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        'printf "readdata {input.out} name mdout \nwritedata {output.rho} mdout[Density] time 0.002\n" > {output.cpptraj} & cpptraj -i {output.cpptraj} -o {log}'


rule md_readPress:
    """
    Read the pressure from a MD.out file
    """
    input:
        out="{somepath}/{file}.out",
    output:
        cpptraj=temp("{somepath}/{file}_press_cpptraj.in"),
        press="{somepath}/{file}_press.dat",
    log:
        "{somepath}/{file}_press.log",
    conda:
        "../envs/ambertools.yml",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        'printf "readdata {input.out} name mdout \nwritedata {output.press} mdout[PRESS] time 0.002\n" > {output.cpptraj} & cpptraj -i {output.cpptraj} -o {log}'


rule md_em_eq_analysis:
    """
    Plot various em/eq plots to assess successful equilibration of the system
    """
    input:
        em_epot="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/4_minim_3/em_3_Epot.dat",
        eq_1_epot="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/5_equil_1/eq_1_Epot.dat",
        eq_1_T="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/5_equil_1/eq_1_temp.dat",
        eq_1_p="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/5_equil_1/eq_1_press.dat",
        eq_2_epot="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/6_equil_2/eq_2_Epot.dat",
        eq_2_T="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/6_equil_2/eq_2_temp.dat",
        eq_2_p="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/6_equil_2/eq_2_press.dat",
        eq_2_rho="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/6_equil_2/eq_2_rho.dat",
        eq_3_epot="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/7_equil_3/eq_3_Epot.dat",
        eq_3_T="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/7_equil_3/eq_3_temp.dat",
        eq_3_p="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/7_equil_3/eq_3_press.dat",
        eq_3_rho="{data_dir}/interim/{md_method}/{compound_dir}/{solvent}/7_equil_3/eq_3_rho.dat",
    output:
        em_plot=report(
            "{data_dir}/processed/{md_method}/results/{compound_dir}/{solvent}/em_eq/em.png",
            category="Compound {compound_dir}",
            subcategory="EM",
        ),
        eq_1_plot=report(
            "{data_dir}/processed/{md_method}/results/{compound_dir}/{solvent}/em_eq/eq_1.png",
            category="Compound {compound_dir}",
            subcategory="EQ",
        ),
        eq_2_plot=report(
            "{data_dir}/processed/{md_method}/results/{compound_dir}/{solvent}/em_eq/eq_2.png",
            category="Compound {compound_dir}",
            subcategory="EQ",
        ),
        eq_3_plot=report(
            "{data_dir}/processed/{md_method}/results/{compound_dir}/{solvent}/em_eq/eq_3.png",
            category="Compound {compound_dir}",
            subcategory="EQ",
        ),
    log:
        notebook="{data_dir}/processed/{md_method}/notebooks/{compound_dir}/{solvent}/em_eq/em_eq.ipynb",
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/template-analyse-em_eq.py.ipynb"
