## Computational Workflows
The workflows consist of the following steps:
1. Inputs & Configuration:
2. MD simulations
3. Conformer generators
4. Per compound analysis
5. Per method analysis

### Inputs
To simulate a compound, a new row needs to be added to the `samples.tsv` file. The row needs to contain the following information:

* `compound-id`: unique identifier for the compound
* `method`: MD-method used for the simulation
* `solvent`: solvent to use for the simulation. Options: aqueous, DMSO, DMSO_GAFF_BCC, DMSO_GAFF_RESP, Chloroform
* `simtime`: simulation time in ns
* `dt`: timestep in ps
* `repeats`: number of simulation repeats
* `other_param`: This is a dummy parameter, changing it from 0.2 leads to creation of a different hash, so that the simulation is run again without overwriting the previous results.
* `igamd`: boosting parameter for GaMD. 1: dihedral angles boost only, 2: total energy boost only, 3: both dihedral and total energy boost

The workflow automatically expands rows found in `samples.tsv` and produces a
unique hash for every MD simulation (see the below figure). The hash is used to name the output files.
Hashes can be found together with the parameters that were set in the 
`samples_old.tsv` file.

![Expanding rows in samples.tsv](docs/images/15_hash%20examples.png)

### Configuration
Via the configuration file `snakemake-config.yaml` a number of parameters can be set.
