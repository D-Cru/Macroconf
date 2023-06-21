

# MacroConf - dataset & workflows for cyclic peptide solutions structure conformer generators

Welcome to MacroConf. This repository contains the MacroConf dataset and associated workflows for analysing cyclic peptide solution structures. The dataset contains NOE distance constraints of cyclic peptides and macrocycles and was assembled from existing literature. The workflows are based on [Snakemake](https://github.com/snakemake/snakemake), a workflow management system. 

------------

# Documentation References
* [Getting started](#getting-started)
* [MacroConf dataset](docs/dataset.md)
* [Computational workflows](docs/workflows.md)
* [Molecular Dynamics](docs/MD.md)
* [Cheminformatics conformer generators](docs/cheminformatics.md)
* [Other details](docs/other.md)
* [Project organization](#project-organization)
* [Copyright and licenses](#copyright-and-licenses)

------------

# Readme Contents

- [MacroConf - dataset \& workflows for cyclic peptide solutions structure conformer generators](#macroconf---dataset--workflows-for-cyclic-peptide-solutions-structure-conformer-generators)
- [Documentation References](#documentation-references)
- [Readme Contents](#readme-contents)
- [Getting started](#getting-started)
  - [Explore the dataset \& code in the browser](#explore-the-dataset--code-in-the-browser)
  - [Download workflow with 2 example compounds](#download-workflow-with-2-example-compounds)
  - [Download workflow with all results](#download-workflow-with-all-results)
  - [Using the MacroConf workflows](#using-the-macroconf-workflows)
- [Project Organization](#project-organization)
- [Copyright and licenses](#copyright-and-licenses)

# Getting started

There are several ways to get started with MacroConf. You could do any of the below steps to interact with the dataset and workflows. The steps are ordered by increasing effort and more demanding computational requirements (memory and compute power required).


## Explore the dataset & code in the browser

If you are interested in the details of the dataset that were already computed,
results are available in the [`workflow/data/processed/refactor-test/`](workflow/data/processed/refactor-test) directory. In addition, you can
access key simulation notebooks directly in your browser using [jupyterbooks](https://danie.lc/Macroconf).
This works without installing anything on your system, and some of the interactive
elements in the notebooks still work. However, you will not be able to change the code that was used to analyze the data.

## Download workflow with 2 example compounds

If you want to explore the computational workflows, you can get the code with a small subset of the simulated data
(2 compounds only). This is the easiest way to get started with the workflows.

The following prerequisites are required to run the steps below:
* Git v 2.25.0 (or higher). You can check your version with `git --version`. To [upgrade git on Ubuntu](https://unix.stackexchange.com/questions/33617/how-can-i-update-to-a-newer-version-of-git-using-apt-get), run 
   ```
   sudo add-apt-repository ppa:git-core/ppa -y
   sudo apt-get update
   sudo apt-get install git -y
   git --version
   ```
* Conda to install required software dependencies: [Installation Instructions](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)
* Amber18 to run the MD simulations: [Installation Instructions](https://ambermd.org/Installation.php)
* OpenEye toolkits and OpenEye Omega if you want to use the cheminformatics conformer generators. [OMEGA](https://www.eyesopen.com/omega) [Installation Instructions](https://docs.eyesopen.com/applications/gettingstarted.html)

Start by cloning the repository, without downloading any of the data:

```git clone --filter=blob:none --no-checkout git@github.com:d-cru/03_macroconf```

```cd 03_macroconf```

After changing into the empty directory, run the following commands:

```git sparse-checkout set --cone```

```git checkout main```

Now you find the files `LICENSE`, `README.md`, `bootstrap.sh` in the directory.

Use the `bootstrap.sh` to download the relevant part of the MacroConf dataset and workflows:

* run `sh bootstrap.sh dataset` to download only the MacroConf dataset

* run `sh bootstrap.sh workflow` to download only the workflow code

* run `sh bootstrap.sh workflow-examples` to download the workflow with 2 examples

* run `sh bootstrap.sh workflow-full` to download the full workflow with all data.
This step is equivalent to downloading the entire repository (see step 3 below).

## Download workflow with all results

If you want to reproduce the analysis results yourself, you can clone the full repository. Due to the large size of the MD trajectories (~TBs), the repository does not contain the MD trajectories. However, the github repository contains shortened versions of the trajectories that can be used to test the workflows. The exact figures in the notebooks will not be reproduced, but the general trends should be the same.

Pre-requisites are the same as in step 2

* Start by cloning the full repository


## Using the MacroConf workflows

* Change to the `workflow` sub-directory and run `conda env create -f envs/snakemake.yml` to create a conda environment with Snakemake.
* Activate the conda environment with `conda activate snakemake`.
* For all the following steps, make sure that you execute those in the activated snakemake environment
* Configure the workflow
   The samples.tsv file contains the list of compounds that were run and analyzed as part of the MacroConf publication.
   In the example version of the workflow, only the data for two compounds are provided (22,24). To not run the full workflow, you need to make the following changes to the `workflow/snakemake-config.yaml` file.


   * set `sample_file: "samples_example.tsv"`
   * set `exp_name: "example"`
   * set `hash_list` to:
      ```
         hash_list: [
         ["ed6dd3148ef9b069", "a08d58e2ed091fdb", "3471e7934b649843", "22", "0", "0"],
         ["ed6dd3148ef9b069", "568b39635b10fbba", "c41b70cc89a57167", "22", "0", "0"],
         ["ed6dd3148ef9b069", "568b39635b10fbba", "c41b70cc89a57167", "22", "omega_basic", "rdkit_ETKDGv3mmff"],
         ["ed6dd3148ef9b069", "8374da5ea9cc35c0", "10236f52ef18b2a3", "22", "0", "0"],
         ["e8dd07e57d6e2799", "056d940c3b10c68b", "056d940c3b10c68b", "24", "0", "0"],
         ]
      ```

* Run the workflow via: `snakemake --cores 1 --use-conda --configfile snakemake-config.yaml --snakefile Snakefile -n` 
(`-n` is a dry-run, remove it to actually run the workflow)

See the following pointers for more details on the workflow:

* If you are interested in reproducing the results, use the `--forceall` flag to force re-running all steps via Snakemake.
  For more details concerning the workflow, see the [workflow documentation](docs/workflows.md).
* If you are interested in running the workflow on HPC clusters, see the [HPC documentation](docs/HPC.md).

* If you are interested in e.g. adding a new conformer generator, you can find more details in the [cheminformatics documentation](docs/cheminformatics.md).

* To extend the dataset and simulate new compounds, you can follow the [dataset documentation](docs/dataset.md).



![Overview of Macroconf dataset & pipeline](workflow/docs/images/4_MacroConf%20Overview%20Fig.png "MacroConf dataset & workflows")


------------


# Project Organization

```
├── data                   <- MacroConf dataset
│   ├── inputs             <- Information about every compound
│   ├── processed NOEs     <- Processed NOE files (matching topologies)
│   ├── raw NOEs           <- Raw NOE files, extracted from papers & cleaned
│   ├── topologies         <- Reference topologies for every compound
│   └── dataset.csv        <- Overview of the MacroConf dataset
│
├── docs                   <- Documentation
│
├── workflow              <- Snakemake workflow
│   ├── data               <- data
│   │   ├── external       <- Contains datasets to run workflow on
│   │   ├── interim        <- Interim files (MD outputs, trajectories, etc.)
│   │   └── processed      <- Workflow outputs, figures, etc.
│   │
│   ├── docs               <- documentation
│   ├── envs               <- Contains conda environment specifications.
│   ├── hpc                <- Contains configuration files for HPC runs.
│   ├── libs               <- External libraries/scripts and parameter files
│   ├── notebooks          <- ipynb notebook templates and notebooks for workflow 
│   ├── reports            <- figures from the paper and other reports
│   ├── rules              <- snakemake rules
│   ├── scripts            <- scripts used by the workflow
│   ├── src                <- Source code for use in this project.
│   │   │
│   │   ├── dihedrals.py   <- Computing dihedral angles
│   │   ├── noe.py         <- Computing NOEs
│   │   ├── pca.py         <- Performing PCAs
│   │   ├── pyreweight.py  <- Reweighting
│   │   ├── Ring_Analysis.py   <- Cremer Pople analysis
│   │   ├── Ring_Reconstruction.py   <- Cremer Pople analysis II 
│   │   ├── stats.py       <- Computing statistics
│   │   └── utils.py       <- Other utilities
│   │
│   ├── tests              <- tests
│   ├── samples_old.tsv    <- archive of all MD simulations run, incl. hashes
│   ├── samples.tsv        <- specification of which MD simulations to run
│   ├── Snakefile          <- Snakefile for main Snakemake workflow (simulations
│   │                         and analysis)
│   ├── snakemake-config.yaml <- Config file to specify parameters
│   └── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
│
├── LICENSE
├── README.md          <- The top-level README for developers using this project.
```

# Copyright and licenses
Licenses for external libraries that were used in this project (& modified):

* [CycloPs: A Cyclic Peptide Library Generator / PepLibGen](https://github.com/fergaljd/cyclops) 
  \- GNU General Public License v2.0. 

  **Citation**:
  CycloPs: generating virtual libraries of cyclized and constrained peptides including nonnatural amino acids FJ Duffy, M Verniere, M Devocelle, E Bernard, DC Shields, AJ Chubb Journal of chemical information and modeling 51 (4), 829-836

  **Modifications**: 
  - Update scripts to work with python3
  - bug fixes (error in some peptide SMILES string)

  located in `workflow/libs/peplibgen/`

* [PyReweighting](http://miaolab.org/PyReweighting/)
   PyReweighting: Python scripts used to reweight accelerated molecular dynamics simulations. \
   Authors: Yinglong Miao <yinglong.miao@gmail.com> \
            Bill Sinko <wsinko@gmail.com>

   Copyright <2014-2019> Yinglong Miao and William Sinko

   **Citation**:
      Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689.

   **Modifications**:
   - Update scripts to work with most recent scipy version 
   - conversion to function based scripts that can be imported into other scripts.

   located in `workflow/src/pyreweight.py`

* [RING](https://github.com/lucianlschan/RING)
   Understanding Ring Puckering in Small Molecules and Cyclic Peptides. \
   Authors: Lucian Chan \
            Geoffrey R. Hutchison \
            Garrett M. Morris

   No code License, but paper available under CC BY NC ND 4.0 License 

   **Citation**:
      J. Chem. Inf. Model. 2021, 61, 2, 743–755

   **Modifications**:
   - Update scripts to work with mdtraj objects used here
   - small other changes

   located in `workflow/src/Ring_Analysis.py` and `workflow/src/Ring_Reconstruction.py`


* [customETKDG](https://github.com/rinikerlab/customETKDG)
   Incorporating NOE-Derived Distances in Conformer Generation of Cyclic Peptides with Distance Geometry. \
   Authors: Wang, Shuzhe \
            Krummenacher, Kajo \
            Landrum, Gregory A \
            Sellers, Benjamin D \
            Di Lello, Paola \
            Robinson, Sarah J \
            Martin, Bryan \
            Holden, Jeffrey K \
            Tom, Jeffrey YK \
            Murthy, Anastasia C \
            Popovych, Nataliya \
            Riniker, Sereina
   
   No code license, but paper available under CC BY NC ND 4.0 License
  
   **Citation**:
      J. Chem. Inf. Model. 2022, 62, 3, 472–485

   **Modifications**:
   Parts of the code were used and modified (bundling methods: NAMFIS, LICUV)
   - Update scripts to work with mdtraj objects used here
   - small other changes

   located in `workflow/notebooks/confgen_NOE.py.ipynb`


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
