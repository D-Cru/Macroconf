# MacroConf Dataset
## Overview
The MacroConf dataset is a collection of cyclic peptide and macrocycle 
structures with experimental NOE distance constraints. The dataset is assembled
from existing literature. An overview of the dataset is available 
[here](data/external/22-09-2021_MacroConf-v2.1/MacroConf-data.xlsx).

Every compound has a `{compound_id}_data.yaml` file in the `data/properties/` directory.
The `{compound_id}_data.yaml` file contains the following information:

```
{
  "index": 22,
  "seq_length": 5,
  "bonds": ["HT"],
  "sequence": "cyclo-(-Ser-Pro-Leu-Asn-Asp-)",
  "sequence_1": "SPLND",
  "natural_cyclic_peptide": 1,
  "non_natural_cyclic_peptide": 0,
  "solvent": "DMSO",
  "smile": "N2[C@@]([H])(CO)C(=O)N1[C@@]([H])(CCC1)C(=O)N[C@@]([H])(CC (C)C)C(=O)N[C@@]([H])(CC(=O)N)C(=O)N[C@@]([H])(CC(=O)O)C2(=O)"
 }

```

In addition, every compound has a `{compound_id}_noe.json` file in the `data/processed_NOEs/` directory, containing the processed NOE distance constraints. The atom indices match the atom indices of the topology file in the `data/topologies/` directory. The `{compound_id}_noe.json` files were derived from the raw NOE distance constraints in the `data/raw_NOEs/` directory. The raw NOE distance constraints were extracted from original publications. Notebooks to reproduce the topology matching can be found in the `workflow/notebooks/` directory. The NOE-topology matching process was performed for all cyclic peptides with natural (L/D) amino acids only. 

## Adding a compound to the MacroConf dataset.
Following a description of how to add a new datapoint and how the datapoints were added. To propose addition of a new compound to the dataset, please open an issue on the [GitHub repository](https://github.com/D-Cru/03_macroconf/issues/new/choose), filling in the `Add a new compound to the dataset` template.

**Step 1**: Create a file: `<compound-id>_in.json` in the dataset folder (`data/interim/refactor-test/`):

A template for this file is supplied within the folder as `template.json` or `template_cis-trans.json`. 

For compound 22 this looks like this:

```
{"index": 22, "seq_length": 5, "bonds": ["HT"],
"sequence": "cyclo-(-Ser-Pro-Leu-Asn-Asp-)",
"sequence_1": "SPLND", "natural_cyclic_peptide": 1,
"non_natural_cyclic_peptide": 0, "solvent": "DMSO"}
```

* *index*: dataset index
* *seq_length*: length of the amino acid sequence
* *bonds*: specify any kind of bonds other than the standard linear peptide bonds
  If the peptide is head to tail cyclic (most of these are...), the first element of bonds
  must be "HT".
  If the peptide is instead cyclised via a disulfide bond, specify the bond as follows:
  String with SS for disulfide bond, then followed by one character per residue in the sequence, which is 'X' for residues not involved in constraints and 'C' for the disulphide bonding residues.
  See `PepLibGen` for more details...
  Example:
  The peptide `CLLRMRSIC`, which is cyclic via a disulfide bond between r.1 and r.9
  can be specified with the string: "SSCXXXXXXXC".

* *sequence*: comprehensive sequence for human readers
* *sequence_1*: 1 letter amino acid sequence (L-aa: UPPERCASE, D-aa: lowercase) (from this SMILES will be generated) . Important here: **z** is the amino acid 2TL (D-AlloThreonine).
* *natural_cyclic_peptide*: specify whether compound only contains L-amino acids
* *non_natural_cyclic_peptide*: specify whether compound only contains L/D amino acids
* *solvent*: specify solvent the NMR experiment from which the NOEs are from
  was performed in: options: aqueous, DMSO, ...

in case of a cis/trans compound, additional fields need to be supplied:
* *multi*: specify names for the subcompounds (e.g. a,b for trans, cis), e.g.: {"52a": "trans", "52b": "cis"}
* *distinction*: specify between which amino acids the cis/trans flip occurs, e.g. [5,1] if between aa:5 - aa:1

We will later generate the SMILES string automatically from the sequence_1 string. This is done via peplibgen [CITATION].
TODO: Cite CycloPs

In case of a structure that cannot be generated with peplibgen (e.g. HT + SS bond), we can add a custom smiles via:
"custom_smiles": "SMILES"
This is then recognized and copied over to smiles.smi, data.json files in rule make_smile.

**Step 2a**: Add an entry for the compound into the `samples.tsv` file. Make sure to use actual tabs as separators.

Then run `mkdir -p data/interim/<exp-name>/<compound-id>/ && touch data/interim/<exp-name>/<compound-id>/NOE.json`.
Then run parts of the snakemake workflow via: `snakemake --until md_make_topology` to generate a SMILEs string, the Amber topology
   and other files that are required for the cheminformatics conformer generators, as well as, for the later analysis steps.

**Step 2b**: In case of a d-amino acid cyclic peptide, do the following:
Add an entry for the compound into the `samples.tsv` file.
Then run `mkdir -p data/interim/<exp-name>/<compound-id>/ && touch data/interim/<exp-name>/<compound-id>/NOE.json`.
Then run parts of the snakemake workflow via: `snakemake --until ob_convert` to generate a SMILEs string, and an initial 3d structure via OpenBabel. Often, the stereochemistry can be wrong here, unfortunately. If this is the case, you need to generate an initital structure via an alternative route (e.g. take an experimental structure if available, Omega, RDKit, ...). Manually create this file as `ob_exp.pdb`, or `ob_d.pdb`. The workflow recognizes presence of either of these files and gives priority as follows: `ob_exp.pdb` > `ob_d.pdb` > `ob.pdb`. For `ob_exp.pdb`, or `ob_d.pdb`, the following rule that uses `pdb4amber` does **not** remove the H atoms, since tleap re-adds those according to a structural template for L-amino acids, and it is out of our control whether the equilibration then leads to correct stereochemistry or not. After manual intervention, run `snakemake --until md_make_topology` to create the amber topology and other files that are required for the chemical informatics conformer generators, as well as, for the later analysis steps.
Make sure to double-check the stereochemistry in step 3 (images are R/S annotated), as well as in the analysis notebooks later (post eq-structure and ref-structure are shown with R/S annotations).

Notes on protonation states and H-atom names:
In some cases the H-atom names are different in Openbabel than what Amber expects (since the new update of Openbabel) (e.g. H1, H2 instead of H2, H3). These need to be fixed manually in the ob_d.pdb file. Sometimes, Openbabel protonates a amino
acid, e.g. adds HE2/HD2 for GLU/ASP residues. These are called GLH/ASH in AMBER. This needs to be fixed in the ob.pdb/ob_d.pdb file manually. Otherwise errors will follow. [See here for more details](http://archive.ambermd.org/201807/0152.html) 

**Step 3**: Add NOEs via jupyter notebook: make a copy of the template `2.0_template` or `2.0template_cis-trans` notebook and run through the script, changing details where appropriate.
  For this step, use a conda environment built from the `envs/mol-maker.yaml` file via:
  `conda env create --file mol-maker.yaml`

  The NOE table from a publication needs to be provided in the following format:
  `Atom name 1 | atom name 2 | available NOE distances (this can be up to 3 columns
    with minimum distance, max. distance, NOE distance in any order)`
  
  For the MacroConf dataset the `.csv` files are kept in the data/external subfolder under the respective paper entries.
  
  A simple way to extract such a `.csv` file from a paper is to use the programme
  `Tabula`, followed usually by small manual corrections.
  The file needs to be provided as .csv

  Then the function getNOE() from src can be used to semi-automatically match the
  atom names in the paper, with the atom names in the Amber topology. The function
  queries the user for atom numbers. Either, indices of the H-atoms can be input directly
  (if multiple, separted via `,`), **or** (often more convenient) the linked C-atom(s)
  numbers can be put in instead. The function automatically replaces this with
  all first-order linked H-atom numbers. (e.g. input atom number of C-alpha, script
  replaces this with the alpha-H-atom).

  To avoid typos and user error, the function queries every atom name found in the .csv file
  twice. If discrepencies are detected between the 1st and 2nd entry, a third
  (and final authoritive) query needs to be provided. At the end, the function outputs a
  NOE dataframe object, as well as a dictionary of all inputs and matched atom names.
  The dictionary should be saved in the jupyter notebook in the cell above, in case the function needs to be re-run
  or if any errors need to be fixed.


**Step 4**: Add <compound-id> to the `samples.tsv` file for simulation.
