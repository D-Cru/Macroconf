# All rules to run OMEGA in the macrocycle mode.


rule omega_flipper:
    input:
        smile="{data_dir}/{compound_dir}/ref_smiles.smi",
    output:
        flipped_smile="{data_dir}/{compound_dir}/flipped.smi",
    log:
        "{data_dir}/{compound_dir}/flipper.log",
    params:
        prefix="{data_dir}/{compound_dir}/",
    envmodules:
        "my.modules",
        "OEApplications/2022.2.1",
    shell:
        "flipper -in {input.smile} -out {output.flipped_smile} -prefix {params.prefix}"


rule omega_first_flip:
    input:
        flipped_smile="{data_dir}/{compound_dir}/flipped.smi",
    output:
        omega_smile="{data_dir}/{compound_dir}/omega.smi",
    log:
        "{data_dir}/{compound_dir}/first_flip.log",
    shell:
        "head -n 1 {input.flipped_smile} > {output.omega_smile}"


rule omega_confgen:
    input:
        omega_smile="{data_dir}/{compound_dir}/omega.smi",
    output:
        omega_bin="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.oeb.gz",
    params:
        omega_param="libs/omega/{omega_type}.param",
        prefix="{data_dir}/{compound_dir}/omega/",
    log:
        "{data_dir}/{compound_dir}/omega/{omega_type}/log.log",
    threads: 4
    envmodules:
        "my.modules",
        "OEApplications/2022.2.1",
    shell:
        "oeomega macrocycle -in {input.omega_smile} -out {output.omega_bin} -param {params.omega_param} -prefix {params.prefix} 2> {log} -mpi_np {threads}"  #


rule omega_extract_mol_pdb:
    input:
        in_file="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.oeb.gz",
    output:
        out_file="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.pdb",
        conf_energy=(
            "{data_dir}/{compound_dir}/omega/{omega_type}/conf_energies.txt"
        ),
    log:
        "{data_dir}/{compound_dir}/omega/{omega_type}/extract_mol_pdb.log",
    conda:
        "../envs/oechem.yaml"
    script:
        "../scripts/oechem_bin_to_mol2.py"


rule omega_extract_mol_mol2:
    input:
        in_file="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.oeb.gz",
    output:
        out_file="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.mol2",
    log:
        "{data_dir}/{compound_dir}/omega/{omega_type}/extract_mol_mol2.log",
    conda:
        "../envs/oechem.yaml"
    script:
        "../scripts/oechem_bin_to_mol2.py"


rule oechm_convert:
    input:
        in_file="{data_dir}/{compound_dir}/H2O/1_make_topology/mc_gas.pdb",
    output:
        out_file=(
            "{data_dir}/{compound_dir}/H2O/1_make_topology/mc_gas_omega.mol2"
        ),
    log:
        "{data_dir}/{compound_dir}/H2O/1_make_topology/oechm_convert.log",
    conda:
        "../envs/oechem.yaml"
    script:
        "../scripts/oechem_bin_to_mol2.py"


rule make_cpptraj:
    input:
        parm="{data_dir}/{compound_dir}/H2O/1_make_topology/mc_sol.prmtop",
        traj="{data_dir}/{compound_dir}/H2O/7_equil_3/traj.ncdf",
    params:
        pdb="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.pdb",
    output:
        cpptraj_file=(
            "{data_dir}/{compound_dir}/H2O/7_equil_3/get_equil_pdb.cpptraj"
        ),
    log:
        "{data_dir}/{compound_dir}/H2O/7_equil_3/make_cpptraj.log",
    conda:
        "../envs/mol-maker.yaml"
    script:
        "../scripts/make_cpptraj.py"


rule pdb_from_eq:
    input:
        cpptraj_file=(
            "{data_dir}/{compound_dir}/H2O/7_equil_3/get_equil_pdb.cpptraj"
        ),
    output:
        pdb="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.pdb",
    log:
        "{data_dir}/{compound_dir}/H2O/7_equil_3/pdb_from_eq.log",
    envmodules:
        "amber/Amber18-AT19-BF17/GCC6.2-CUDA10.1",
    shell:
        "cpptraj -i {input.cpptraj_file}"


rule oechm_convert_equil:
    input:
        in_file="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.pdb",
    output:
        out_file="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.mol2",
    log:
        "{data_dir}/{compound_dir}/H2O/7_equil_3/oechm_convert_equil.log",
    conda:
        "../envs/oechem.yaml"
    script:
        "../scripts/oechem_bin_to_mol2.py"


rule make_cheminfo_smiles:
    input:
        ref_mol="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.mol2",
    output:
        smiles="{data_dir}/{compound_dir}/ref_smiles.smi",
    conda:
        "../envs/mol-maker.yaml"
    script:
        "../scripts/make_smiles_from_mol.py"


rule omega_align_top:
    input:
        pdb="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.pdb",
        mol2="{data_dir}/{compound_dir}/omega/{omega_type}/mcs.mol2",
        ref_top="{data_dir}/{compound_dir}/H2O/1_make_topology/mc_gas.pdb",
        ref_pdb_equil="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.pdb",
        ref_equil_mol2="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.mol2",
    output:
        mol_aligned=(
            "{data_dir}/{compound_dir}/omega/{omega_type}/mcs_aligned.pdb"
        ),
    log:
        notebook="{data_dir}/{compound_dir}/omega/{omega_type}/align_top_log.py.ipynb",
    conda:
        "../envs/mol-maker.yaml"
    notebook:
        "../notebooks/omega_align.py.ipynb"


rule rdkit_confgen:
    input:
        omega_smile="{data_dir}/{compound_dir}/omega.smi",
        # omega="{data_dir}/{compound_dir}/omega/basic/mcs_aligned.pdb",
        # omega="{data_dir}/{compound_dir}/H2O/1_make_topology/mc_gas.pdb",
        omega="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.mol2",
    output:
        pdb="{data_dir}/{compound_dir}/rdkit/{rdkit_type}/mcs.pdb",
        mol2=touch("{data_dir}/{compound_dir}/rdkit/{rdkit_type}/mcs.mol2"),
        mmff_energies=(
            "{data_dir}/{compound_dir}/rdkit/{rdkit_type}/conf_energies.txt"
        ),
    log:
        notebook="{data_dir}/{compound_dir}/rdkit/{rdkit_type}/rdkit_confgen_log.py.ipynb",
    params:
        rdkit_path="libs/rdkit/",
        rdkit_param="libs/rdkit/{rdkit_type}.py",
    threads: 4
    conda:
        "../envs/mol-maker.yaml"
    notebook:
        "../notebooks/rdkit_confgen.py.ipynb"


rule rdkit_align_top:
    input:
        pdb="{data_dir}/{compound_dir}/rdkit/{omega_type}/mcs.pdb",
        mol2="{data_dir}/{compound_dir}/rdkit/{omega_type}/mcs.mol2",
        ref_top="{data_dir}/{compound_dir}/H2O/1_make_topology/mc_gas.pdb",
        ref_pdb_equil="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.pdb",
        ref_equil_mol2="{data_dir}/{compound_dir}/H2O/7_equil_3/equil.mol2",
    output:
        mol_aligned=(
            "{data_dir}/{compound_dir}/rdkit/{omega_type}/mcs_aligned.pdb"
        ),
    log:
        notebook="{data_dir}/{compound_dir}/rdkit/{omega_type}/align_top_log.py.ipynb",
    conda:
        "../envs/mol-maker.yaml"
    notebook:
        "../notebooks/omega_align.py.ipynb"


rule confgen_NOE:
    input:
        pdb="data/interim/{exp_name}/{compound}/{confgen}/{mode}/mcs_aligned.pdb",
        noe="data/interim/{exp_name}/{compound}/NOE.json",
        parm="data/interim/{exp_name}/{compound}/data.json",
        energies="data/interim/{exp_name}/{compound}/{confgen}/{mode}/conf_energies.txt",
    output:
        best_NOE_plot="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/best_NOE.svg",
        NOE_violin_plot="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/NOE_distribution.svg",
        fulfilled="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/NOE_fulfilled.json",
        bundle_plot="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/bundle_plot.svg",
        sasa="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/sasa.json",
        psa="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/psa.json",
        solvation_properties="data/processed/{exp_name}/results/{compound}/conf_gen/{confgen}/{mode}/solvation_properties.json",
    threads: 1
    conda:
        "../envs/stats.yaml"
    log:
        notebook="data/processed/{exp_name}/notebooks/{compound}/conf_gen/{confgen}_{mode}_NOE.py.ipynb",
    notebook:
        "../notebooks/confgen_NOE.py.ipynb"
