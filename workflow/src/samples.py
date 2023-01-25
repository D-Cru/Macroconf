def compute_nstlim(dt, simtime):
    """Compute nstlim from dt and simtime

    Parameters
    ----------
    dt : float
        time step of MD simulation. dt is in ps. (usually 0.002)
    simtime : float
        length of simulation in ns. (usually 2000)

    Returns
    -------
    int
        number of simulations teps


    example: simtime=100 ns * 1000, dt = 0.002

    """

    nstlim = simtime * 1000 / dt
    return int(nstlim)


def lookup_exp_solvent(compound, dataset_file):
    """Find the experimental solvent for a compound in a dataset file

    Parameters
    ----------
    compound : string
        Compound identifier, e.g. "22", or "1_2"
    dataset_file : string
        filepath to dataset file. Must be a csv file with a column
        "solvent data" that contains the experimental solvent for each
        compound.

    Returns
    -------
    string
        Experimental solvent.

    Raises
    ------
    ValueError
        If the compound is not found in the dataset.
    """
    import pandas as pd

    dataset = pd.read_csv(dataset_file)
    try:
        compound_info = dataset[dataset["id"] == compound]
        solvent = compound_info["solvent data"].values[0]
    except IndexError:
        raise ValueError("Compound {compound} not found in dataset.")
    return solvent


def update_hist_file(sample_df, hist_path):
    """Update the history file to contain all samples that are in the
    sample.tsv file.

    Parameters
    ----------
    sample_df : pandas.DataFrame
        DataFrame containing all samples to be run. This df was expanded
        to contain all combinations of parameters as single rows and has
        already a hash as unique identifier for each sample.
    hist_path : string
        Path to history file. This file is a csv file. We append to this file,
        to not loose information about previous runs.
    """
    import pandas as pd
    import os

    if not os.path.isfile(hist_path):
        sample_df.to_csv(hist_path, index=True, sep="\t")
    else:
        old_samples = pd.read_table(hist_path, dtype=str).set_index(
            "hash", drop=True
        )
        combined_samples = pd.concat([sample_df, old_samples])
        combined_samples = combined_samples[
            ~combined_samples.index.duplicated(keep="first")
        ]
        combined_samples.to_csv(hist_path, index=True, sep="\t")


def compute_hash(df, columns):
    """Compute a hash column for each row of a dataframe. Only considers columns
    provided via columns, not all columns of the df

    Parameters
    ----------
    df : Pandas.DataFrame
        DataFrame for which the index should be computed as a hash of the
        columns provided in columns.
    columns : list(string)
        List of columns to be considered for the hash computation.
        All column names must be strings and must be present in the df.

    Returns
    -------
    Pandas.DataFrame
        Updated DataFrame with a new column "hash" that contains the hash.
        The hash column is used as index for the DataFrame.
    """
    import hashlib

    # Compute md5sum for each sample set as index
    # columns = list(df.columns.values)
    for index, row in df.iterrows():
        # do not include nan columns in hash
        columns_to_hash = [c for c in columns if row[c] != "nan"]

        # compute hash over columns_to_hash
        df.loc[index, "hash"] = hashlib.md5(
            str(row[columns_to_hash].values).encode("utf-8")
        ).hexdigest()[:16]
    df = df.set_index("hash", drop=True)
    # df = df.drop_duplicates()
    return df


def adjust_method_specific_parameters(df, config_dict):
    """From a definition of method specific parameters, replace any parameter
    entries in the sample df with valid parameter values.


    Parameters
    ----------
    df : Pandas.DataFrame
        sample DataFrame. This DataFrame contains all samples to be run.
        It can contain method specific parameters, which might be invalid.
    config_dict : dict
        allowed values and defaults for method specific parameters.

    Returns
    -------
    Pandas.DataFrame
        Updated sample DataFrame with valid method specific parameters set.

    Raises
    ------
    ValueError
        If there are NaNs in the sample DataFrame for parameters that are not
        defined in the config_dict.
    """
    # set to supplied value for the method(s), set to NaN for all others,
    for k, v in config_dict.items():
        # k is the dict key, which must correspond to a column of the sample df
        # v is a dict of properties: methods, default, others
        for m in v["methods"]:
            # samples["igamd"] = samples.apply(
            #     lambda x: "3" if (x["method"] ==
            #       "GaMD" and np.isnan(float(x["igamd"]))) else x["igamd"],
            #     axis=1,
            # )
            # Replace nan with default value, if the method corresponds to
            # method specified in dict
            df[k] = df.apply(
                lambda x: v["default"]
                if (x["method"] == m and x[k] == v["others"])
                else x[k],
                axis=1,
            )
        # Replace set value with 'other' value, if the method does not
        # correspond to the methods specified in dict
        # For example: if the samples.tsv specified a row where the igamd
        # parameter is set to 3 for cMD, replace with
        # NaN to avoid duplicate runs that are actually the same!
        df[k] = df.apply(
            lambda x: v["others"]
            if (x["method"] not in v["methods"] and x[k] != v["others"])
            else x[k],
            axis=1,
        )
    # Check if df has nans in columns that should be populated...,
    # if so return error
    nan_columns = df.isin(["nan"]).any()
    for k, v in nan_columns.items():
        if k not in config_dict.keys() and v:
            raise ValueError(
                f"""Column {k} has NaN values, please check samples.tsv file, or
                define a default value in snakemake-config,yaml, key:
                method_defaults."""
            )

    return df
