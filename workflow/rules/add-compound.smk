### Workflow for creating AMBER topology and NOEs.


rule make_Smile:
    input:
        json="{data_dir}/{compound_dir}_in.json",
    output:
        json_all="{data_dir}/{compound_dir}/data.json",
        smile="{data_dir}/{compound_dir}/smile.smi",
        mol_image="{data_dir}/{compound_dir}/structure.png",
    log:
        "{data_dir}/{compound_dir}/make_smile.log",
    conda:
        "../envs/mol-maker.yaml"
    script:
        "../scripts/make_smile.py"


rule ob_convert:
    input:
        smile="{data_dir}/{compound_dir}/smile.smi",
    output:
        pdb="{data_dir}/{compound_dir}/ob.pdb",
    log:
        "{data_dir}/{compound_dir}/ob_convert.log",
    conda:
        "../envs/mol-maker.yaml"
    shell:
        "obabel -i smi {input.smile} -o pdb -O {output.pdb} --gen3D best"


rule read_mol2:
    input:
        mol="{data_dir}/{compound_dir}/1_make_topology/mc_gas.mol2",
    output:
        png="{data_dir}/{compound_dir}/1_make_topology/mc_gas.png",
    log:
        "{data_dir}/{compound_dir}/read_mol2.log",
    conda:
        "../envs/mol-maker.yaml"
    script:
        "../scripts/read_mol2.py"
