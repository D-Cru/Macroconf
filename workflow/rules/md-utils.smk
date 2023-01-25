# Small utility rules


rule util_ncdf:
    input:
        "{some_path}.netcdf",
    output:
        "{some_path}.ncdf",
    log:
        "{some_path}.ncdf_ln.log",
    shell:
        "ln -sr {input} {output}"


rule md_nglviewer: # This rule needs to always run in an actual jupyter notebook browser session.
    # Run via: snakemake --profile lab --edit-notebook {data_dir}/{compound_dir}/{md_dir}/results/comp/visual.html
    input:
        unpack(compare_inputs), # ga_clus=(
         #     "{data_dir}/{compound_dir}/{md_dir}/results/GaMD/clusters/clusters.pdb"
         # ),
         # a_clus="{data_dir}/{compound_dir}/{md_dir}/results/aMD/clusters/clusters.pdb",
         # c_clus="{data_dir}/{compound_dir}/{md_dir}/results/cMD/clusters/clusters.pdb",
         # noe="{data_dir}/{compound_dir}/NOE.json",
    output:
        visual="data/processed/{exp_name}/results/{compound_dir}/comparison/{index_0}_{index_1}_{index_2}-visual.html", #visual="{data_dir}/{md_dir}/results/{md_dir}/comparison/visual_{index_0}_{index_1}_{index_2}.html",
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/template-ngl-visualisation.py.ipynb"


def get_paper_fig_files(wildcards):
    # remove label to extract file name
    file_name = f"{wildcards.fig.split('_', 1)[1]}.{wildcards.extension}"
    # get matching original file name from paper_fig_df
    original_file = paper_fig_df.loc[
        paper_fig_df.file_names == file_name
    ].file_path_name.tolist()[0]
    return original_file


rule paper_figs:
    input:
        get_paper_fig_files,
    output:
        "reports/paper_figures/{fig}.{extension}",
    log:
        "reports/paper_figures/log/{fig}.{extension}.log",
    shell:
        "cp {input} {output}"
