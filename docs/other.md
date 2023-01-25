
# Further details

## Running this workflow on HPC systems
This workflow can be run without adaptations on most HPC/Cloud systems, provided
that the required software is installed (especially AMBER 18) as environment
modules (or otherwise). Additional configuration files are 
[provided](../workflow/hpc) for running this workflow with a SLURM managed HPC
system.

### Interactive jupyter notebooks
To interactively run a jupyter notebook with Snakemake, invoke snakemake via
`snakemake --edit-notebook` and then specifiying the required output file. This
automatically opens a jupyter notebook with an editable copy of the notebook.