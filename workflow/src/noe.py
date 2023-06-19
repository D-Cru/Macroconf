def getNOE(mol, csv_path, atom_dict=None, style=None, scaling=None):
    """Parse a NOE-csv file and return a NOE-dict.
    Parse a csv file containing a NOE table and return a NOE-dictionary that
    matches a given molecular topology specified in mol.

    Parameters
    ----------
    mol : RDKit Mol object
        Corresponding molecule object. This needs to have all H-atoms
        explicitely.
    csv_path : string
        File path to the csv file containing the NOE table extracted from
        literature. The first two columns are the atom-names of the two atoms
        forming the NOE distance. The columns 3-5 can contain the NOE value,
        upper and lower bounds. This should be specified in the `style`
        parameter.
    atom_dict : dict, optional
        This is a dictionary containing the user inputs of previous runs.
        Importantly, the user inputs can be linked atom numbers
        (e.g. C-a atom instead of directly specifying the Ha atom).
        By default None
    style : string, optional
        A 3-character style string, that labels the columns of the supplied
        csv file in csv_path. The first two columns of the csv file are fixed
        to the atom names, but the columns 3-5 can contain NOE distances, NOE
        upper bounds, or NOE lower bounds.

        Possible values are:
            'D' : NOE distance
            'A'  : NOE upper bound (mAximum)
            'I'  : NOE lower bound (mInimum)
            '-'  : Do not use column, or column does not exist if at the end.
        The string 'DIA' would match a file with the columns:
            'atom1','atom2','NOE','NOE_lower','NOE_upper'

        By default None

    scaling : float, optional
        Scale the supplied NOE distances & bounds. Should be set to a value to
        convert to Angstroms
        By default None

    TODO: Handle interupts (e.g. if the user aborts the program).
    TODO: Handle partly filled atom_dict's.

    Returns
    -------
    NOE_dict, df
        Returns a dictionary containing the NOE-values for the given molecule
        (converted from a pandas df). Also returns the pandas df containing the
        NOE-values.
        The df has the following columns:
        Atom 1 | Atom 2 | NMR exp | lower bound | upper bound
        Atom 1, 2: tuple of atom numbers, matching the mol object (md topology)
        NMR exp: NOE distance (can be empty)
        lower bound: NOE lower bound (can be empty)
        upper bound: NOE upper bound (can be empty)
    """
    import pandas as pd
    import numpy as np

    NOE_df = pd.read_csv(csv_path, header=None)

    if style is None:
        while True:
            try:
                style = input(
                    """Enter style of the columns 3-5 of the csv file table:
                    options are: D (distance), I (min), A (max),
                    - (do not use)"""
                )
                assert len(style) == 3
                for i in style:
                    assert i in ["D", "I", "A", "-"]
                break
            except AssertionError:
                print("Invalid style. Try again.")
                continue

    df_column_names = ["atom1", "atom2"]
    df_column_names.extend([char for char in style])
    NOE_df.columns = df_column_names
    # type cast first two columns to string
    NOE_df["atom1"] = NOE_df["atom1"].astype("string")
    NOE_df["atom2"] = NOE_df["atom2"].astype("string")

    # Get all atom names provided in the csv file
    atom_names = list(NOE_df.atom1)
    atom_names.extend(list(NOE_df.atom2))
    atom_names = [str(x) for x in atom_names]
    # Remove duplicates & sort
    atom_names = set(atom_names)
    atom_names = sorted(atom_names)

    # Query the user to provide atom numbers for the atom names in the csv file
    if atom_dict is None:
        atom_dict = {}
        print(
            """Following, enter atom numbers. If multiple, separate with a `,`.
            Can be H atoms or connected atoms"""
        )
        for a in atom_names:
            while True:
                # query the user for the atom number,
                # which can be a list or a single number
                try:
                    atom_number = input(f"Enter atom number(s) for {a}: ")
                    if type(atom_number) is str:
                        if atom_number == "":
                            raise ValueError("Invalid input. Try again.")
                    elif type(atom_number) is int:
                        atom_number = atom_number
                    else:
                        raise ValueError("Invalid input. Try again.")

                    # convert the atom numbers to a list
                    atom_number = list(map(int, atom_number.split(",")))

                    # Check if the input numbers are part of the molecule
                    for i in atom_number:
                        if i > mol.GetNumAtoms():
                            raise ValueError(
                                f"""Atom {i} is not part of the molecule.
                                Try again"""
                            )

                    # Save to atom_dict
                    atom_dict[a] = atom_number
                    break
                except ValueError:
                    print("Invalid atom number. Try again.")
                    continue

        # Now ask again to ensure we did not make any errors:
        for a in atom_names:
            atom_number = input(f"Insert atom number of main atom {a}:")
            atom_number = list(map(int, atom_number.split(",")))

            # Check if this matches the previous input
            if atom_dict[a] == atom_number:
                print("matches previous input!")
            else:
                print("Records do not match!")
                atom_number = input(
                    f"""Insert atom number of main atom {a}:
                    check again carefuly! this will overwrite:"""
                )
                atom_number = list(map(int, atom_number.split(",")))
                atom_dict[a] = atom_number
    # atom_dict contains the user inputs.
    # These can be H atoms or connected atoms. Now retrieve the actual
    # H atom numbers
    h_dict = {}

    # Conversion of atoms numbers via RDKit
    for atom_name, atom_number in atom_dict.items():
        h_atoms = []
        for a in atom_number:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetAtomicNum() != 1:
                for x in atom.GetNeighbors():
                    typ = x.GetAtomicNum()
                    idx = x.GetIdx()
                    if typ == 1:
                        h_atoms.append(idx)
            else:  # supplied atom is H!
                h_atoms.append(atom.GetIdx())
        # Now save this to h_dict:
        h_dict[atom_name] = str(tuple(h_atoms))

    # Now replace atom names by H atom numbers
    NOE_df = NOE_df.replace(h_dict)

    # Set columns to 0 if not determined experimentally
    if "-" in NOE_df.columns:
        NOE_df = NOE_df.drop(columns=["-"])
    if "D" not in NOE_df.columns:
        NOE_df["D"] = 0
    if "I" not in NOE_df.columns:
        NOE_df["I"] = 0
    if "A" not in NOE_df.columns:
        NOE_df["A"] = 0
    cols = ["atom1", "atom2", "D", "I", "A"]
    NOE_df = NOE_df[cols]
    NOE_df.columns = [
        "Atom 1",
        "Atom 2",
        "NMR exp",
        "lower bound",
        "upper bound",
    ]

    # ask for re-scaling:
    if scaling is None:
        scaling = float(input("Please supply a unit-scaling factor:"))

    # Deal with multiple NOE distances & scale to Angstroms
    if NOE_df["NMR exp"].dtype == "object":
        exp_list = np.array(
            NOE_df["NMR exp"].str.split(",")
        )  # split strings into list
        exp_list = NOE_df["NMR exp"].to_numpy()

        for idx, v in enumerate(exp_list):
            l_item = v.split(",")
            l_item = [float(x) for x in l_item]
            exp_list[idx] = np.array(l_item) * scaling
            exp_list[idx] = list(exp_list[idx])

        NOE_df["NMR exp"] = exp_list  # *scaling
    else:
        NOE_df["NMR exp"] = NOE_df["NMR exp"] * scaling
    NOE_df["lower bound"] = NOE_df["lower bound"] * scaling
    NOE_df["upper bound"] = NOE_df["upper bound"] * scaling

    print("NOE dataframe:")
    print(NOE_df)
    print("Dict of user inputs. Save as a variable for potential re-runs.")
    print(atom_dict)
    return NOE_df.to_dict(orient="index"), NOE_df


def NOE_average(x, weight=None):
    """Perform 1/r^6 (NOE) averaging.

    Parameters
    ----------
    x : numpy.array
        1d numpy array of NOE distances.
    weight : numpy.array, optional
        a weight vector (must not be normalized), by default None
        must be of the same size as x.
    Returns
    -------
    float
        NOE average <x^-6>^(-1/6)

    Raises
    ------
    ValueError
        If no numpy array is supplied

    AssertionError
        If x and weight are not of the same size
    """
    import numpy as np

    try:
        type(x) is np.array
        type(weight) is np.array  # or type(weight) is None
    except ValueError:
        raise ValueError("Value Error. Expected numpy array!")

    if weight is not None:
        weight = np.array(weight)
        assert x.shape == weight.shape
        tmp = np.power(x, (-6), dtype=float) * weight
        return (np.sum(tmp) / np.sum(weight)) ** (-1 / 6)
    else:
        tmp = np.power(x, (-6), dtype=float)
        return (np.sum(tmp) / len(x)) ** (-1 / 6)


def compute_NOE_mdtraj(
    NOE_dict,
    traj,
    reweigh_type=0,
    slicer=None,
    weight_data=None,
):
    """Compute NOEs for a mdtraj trajectory

    Computes NOEs for a mdtraj trajectory and corresponding NOE dictionary.

    Parameters
    ----------
    NOE_dict : dict
        a NOE dictionary, created by the src.noe.getNOE() function.
    traj : _type_
        the corresponding mdtraj trajectory (aligned and without solvent).
    weightfile : _type_, optional
        file path of a pyreweight compatible weight file for
        re-weighting. if no weightfile is supplied: run without reweighting,
        by default None
    reweigh_type : int, optional
        0: MacLauren Series expansion. Over all frames. This is not
        a correct way of reweighting NOEs!

        1: 1d PMF of NOE distance via Cummulant expansion to
        the second order., by default 0
    slicer : np.array, optional
        used for cis/trans. Sliced the trajectory according to entries, so only
        consider parts of the trajectory. Bool array.
        , by default None
    weight_data : np.array, optional
        array of weights, used if NOEs should be reweighted, by default None
    cluster_NOE : bool, optional
        _description_, by default False

    Returns
    -------
    dict
        returns the NOEs for MDtraj trajectory, and computed lower and upper
        errors.
    """
    import mdtraj as md
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn import metrics

    # storage for computed NOEs, upper/lower errors, distribution of NOEs
    # and plots
    NOE_md = []  # NOE distances
    NOE_md_up = []  # upper errors
    NOE_md_low = []  # lower errors
    NOE_dist = []  # NOE distributions
    NOE_plot_data = []

    # iterate over individual NOEs
    for i, dic in NOE_dict.items():
        # get atom indices for NOE i
        atom1 = eval(dic["Atom 1"])  # Atom 1, 2 can be multiple atom indices.
        atom2 = eval(dic["Atom 2"])  # in case of amb. NOEs.

        # Check that atom indices are valid
        if len(atom1) == 0 or len(atom2) == 0:
            # at least 1 atom is undefined.
            NOE_md.append(np.nan)
            NOE_dist.append([np.nan])
            print(
                f"NOE {i} error: at least 1 of the \
            atoms is undefined! Set value to NaN"
            )
            NOE_md_up.append(np.nan)
            NOE_md_low.append(np.nan)
        else:
            # ambiguous NOE and unique NOE combined
            NOE_amg = []  # stores the ambiguous NOE values
            NOE_amg_dist = (
                []
            )  # stores the distribution of ambiguous NOE values
            # storage for md upper / lower bounds:
            NOE_up_store = []
            NOE_low_store = []

            # iterate over all combinations of atom indices, in case of
            # ambiguous NOEs, atom1,2 are lists of atom indices.
            for a1 in atom1:
                for a2 in atom2:
                    # Compute distances in Angstrom between atoms a1, a2
                    # over all frames of the traj.
                    dist_1 = (
                        md.compute_distances(traj, np.array([[a1, a2]])) * 10
                    )

                    # Catch case where all distance values are 0(same atom)
                    if (dist_1.max() == dist_1.min()) and dist_1.max() == 0:
                        continue
                    else:
                        # atoms are not the same. -> distance > 0
                        # Save distance distribution
                        NOE_amg_dist.append(dist_1.flatten())

                        # Compute NOE values

                        if reweigh_type == 0:
                            # No reweighting, default. Still uses weights if
                            # supplied from external source.
                            NOE_new = NOE_average(
                                dist_1.flatten(), weight_data
                            )

                        elif reweigh_type == 1:
                            # Reweight with 1d PMF of NOE distance via
                            # Cummulant expansion to the second order.
                            from .pyreweight import reweight_1d_pmf

                            # Compute 1d PMF.
                            pmf, distances = reweight_1d_pmf(
                                dist_1.flatten(),
                                None,
                                "amdweight_CE",
                                slicer,
                                weight_data,
                            )
                            # Convert PMF to weight via Boltzmann factor.
                            weights = np.exp(-1 * pmf / 0.5961)
                            # Normalize weights (should not be necessary..)
                            weights_norm = weights / np.sum(weights)
                            # Perform NOE average.
                            NOE_new = NOE_average(distances[:-1], weights_norm)
                            # save pmf for plotting
                            NOE_plot_data.append((pmf, distances[:-1], i))

                        elif reweigh_type == 2:
                            # Reweigh by using 1d MC expansion. This should not
                            # be used.
                            from .pyreweight import reweight_1d

                            weight_distances = reweight_1d(
                                dist_1.flatten(), None, "amdweight_MC", slicer
                            )
                            # Perform NOE average.
                            NOE_new = NOE_average(
                                dist_1.flatten(), weight_distances
                            )

                        elif reweigh_type == 3:
                            # Reweight by using 1d PMF of NOE distance via
                            # Cummulant expansion to the second order,
                            # if weight_data is supplied. Else don't reweigh
                            # However, don't use NOE averaging. Instead, use
                            # a standard average. This is experimental for use
                            # of cluster bundles.
                            from .pyreweight import reweight_1d_pmf

                            # No reweighting. Use std averaging, instead of
                            # r^-6 average. Only for cheminformatics clusters!
                            if weight_data is None:
                                weights = np.ones(dist_1.flatten().shape)
                                weights = weights / np.sum(weights)

                            else:
                                # Compute 1d PMF.
                                pmf, distances = reweight_1d_pmf(
                                    dist_1.flatten(),
                                    None,
                                    "amdweight_CE",
                                    slicer,
                                    weight_data,
                                )
                                #    Convert PMF to weight via Boltzmann factor
                                weights = np.exp(-1 * pmf / 0.5961)
                                # Normalize weights (should not be necessary..)

                            NOE_new = np.average(
                                dist_1.flatten(), weights=weights
                            )

                        elif reweigh_type == 4:
                            # Reweight with 1d PMF of NOE distance via
                            # Maclaurin Series expansion
                            from .pyreweight import reweight_1d_pmf

                            # Compute 1d PMF.
                            pmf, distances = reweight_1d_pmf(
                                dist_1.flatten(),
                                None,
                                "amdweight_MC",
                                slicer,
                                weight_data,
                            )
                            # Convert PMF to weight via Boltzmann factor.
                            weights = np.exp(-1 * pmf / 0.5961)
                            # Normalize weights (should not be necessary..)
                            weights_norm = weights / np.sum(weights)
                            # Perform NOE average.
                            NOE_new = NOE_average(distances[:-1], weights_norm)
                            # save pmf for plotting
                            NOE_plot_data.append((pmf, distances[:-1], i))

                    # Save computed NOE average to ambiguous NOE list.
                    NOE_amg.append(NOE_new)

                    # Compute upper / lower bounds
                    if len(dist_1.flatten()) > 1:
                        # upper error bound
                        NOE_i = dist_1.flatten()[dist_1.flatten() > NOE_new]
                        # if there are no values that fulfill above criterium
                        # -> error = 0 (which means NOE average is equal to
                        # smallest/biggest value)
                        if NOE_i.size == 0:
                            NOE_i = np.array([NOE_new])
                        NOE_up = (
                            metrics.mean_squared_error(
                                np.full(len(NOE_i), NOE_new),
                                NOE_i,
                                squared=False,
                            )
                            + NOE_new
                        )
                        # do I do this with the reweighted values? or with all
                        # of the values as observed in (G)aMD?
                        # My interpretation of Kamenik et al.: not reweighted..

                        # lower error bound
                        NOE_i = dist_1.flatten()[dist_1.flatten() < NOE_new]
                        # if there are no values that fulfill above criterium
                        # -> error=0
                        if NOE_i.size == 0:
                            NOE_i = np.array([NOE_new])
                        NOE_low = NOE_new - metrics.mean_squared_error(
                            np.full(len(NOE_i), NOE_new), NOE_i, squared=False
                        )

                        # Save upper / lower bound of computes NOE
                        NOE_up_store.append(NOE_up)
                        NOE_low_store.append(NOE_low)
                    else:
                        # catch cases where NOE is computed for 1 frame only.
                        NOE_up_store.append(None)
                        NOE_low_store.append(None)

            # We do not average the ambiguous NOEs.
            # (This could be done via again using NOE_average())
            # NOE_md.append(NOE_average(NOE_amg))
            NOE_md.append(NOE_amg)
            NOE_dist.append(NOE_amg_dist)

            NOE_md_up.append(NOE_up_store)
            NOE_md_low.append(NOE_low_store)

    # Produce a PMF plot for every NOE
    if reweigh_type == 1 or reweigh_type == 4:
        fig, axs = plt.subplots(
            int((len(NOE_plot_data) / 5) + 1), 5, squeeze=True
        )
        fig.set_size_inches(6, int((len(NOE_plot_data) / 5) + 1))
        for ax, data in zip(axs.flat, NOE_plot_data):
            pmf, distances, name = data
            ax.plot(distances, pmf)
            ax.set_title(f"NOE {name}")
        for i in range(len(NOE_plot_data), len(axs.flat)):
            axs.flatten()[i].axis("off")
        fig.supxlabel("Distance (Angstrom)")
        fig.supylabel("PMF (kcal/mol)")
        fig.tight_layout()
        NOE_plot = fig
    else:
        NOE_plot = None

    return (NOE_md, NOE_md_low, NOE_md_up, NOE_dist, NOE_plot)


def plot_NOE(NOE_df, fig=None, ax=None, default=None):
    """Make a NOE plot from a NOE dataframe.

    Parameters
    ----------
    NOE_dict : pandas.DataFrame
        a NOE dataframe, created by the src.noe.getNOE() function.
    fig : matplotlib.figure, optional
        Can supply a matplotlib figure object to use for plotting, by default
        None
    ax : matplotlib.axes, optional
        Can supply a matplotlib ax object to use for plotting, by default None
    default : any, optional
        Set this to use an alternative plotting option, by default None

    Returns
    -------
    (matplotlib.figure, matplotlib.axes)
        Returns a matplotlib figure and axes object, containing the NOE plot.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Check if a figure and axes object are supplied. Use supplied objects only
    # if both are supplied. If not, create new ones.
    if (fig is None) or (ax is None):
        if (fig is not None) or (ax is not None):
            raise Warning(
                "Supply both fig and corresponding axs. Created new ones."
            )
        fig, ax = plt.subplots(figsize=(6.7323, 3.2677))  # figsize=(7, 3)

    if default is None:
        # Default plotting option.
        # Plot difference between simulation average and experiment

        # Generate x labels. Repeat ambiguous NOE labels, e.g. 1,1,1,2,3...
        xlabels = [str(x) for x in NOE_df.index]
        xvalues = np.arange(len(NOE_df.index))

        # Plot experimental reference line
        ax.plot(xvalues, np.zeros(len(NOE_df.index)), color="black")

        # If no experimental values exist, use middle between upper and lower
        # bound as reference line.
        if all(v == 0 for v in NOE_df["NMR exp"]):
            NOE_df["NMR exp"] = (
                NOE_df["upper bound"] + NOE_df["lower bound"]
            ) * 0.5

        # Plot experimental values
        ax.scatter(
            xvalues,
            NOE_df["md"] - np.abs(NOE_df["NMR exp"]),
            color="darkred",
            marker="x",
            label="deviation",
        )

        # Plot experimental NOE upper bound, if exists
        if not all(v == 0 for v in NOE_df["upper bound"]):
            # should not need an if here, b/c if exp not there, values are 0!
            # if not all(v == 0 for v in NOE_dict['NMR exp']):
            ax.plot(
                xvalues,
                NOE_df["upper bound"] - np.abs(NOE_df["NMR exp"]),
                label="NMR upper bound",
                color="lightgrey",
            )

        # Plot experimental NOE lower bound, if exists
        if not all(v == 0 for v in NOE_df["lower bound"]):
            ax.plot(
                xvalues,
                NOE_df["lower bound"] - np.abs(NOE_df["NMR exp"]),
                label="NMR lower bound",
                color="lightgrey",
            )

        # Plot MD errors:
        if "upper" in NOE_df.columns:
            if not all(v == 0 for v in NOE_df["upper"]):
                ax.scatter(
                    xvalues,
                    NOE_df["upper"] - np.abs(NOE_df["NMR exp"]),
                    label="MD upper bound",
                    color="blue",
                    marker="_",
                )
        if "lower" in NOE_df.columns:
            if not all(v == 0 for v in NOE_df["lower"]):
                ax.scatter(
                    xvalues,
                    NOE_df["lower"] - np.abs(NOE_df["NMR exp"]),
                    label="MD lower bound",
                    color="blue",
                    marker="_",
                )
        # Overwrite xticks to account for ambiguous NOEs
        ax.set_xticks(xvalues, xlabels)
        # Make x-ticks staggered
        for tick in ax.xaxis.get_major_ticks()[1::2]:
            tick.set_pad(15)
        ax.set_xlabel("NOE no.")
        ax.set_ylabel(r"MD - NMR(exp): distance / $\AA$")

    else:
        # Use alternative plotting option
        # Plot NOE values
        ax.scatter(
            [str(x) for x in NOE_df.index],
            NOE_df["md"],
            color="darkred",
            marker="x",
            label="MD distance",
        )

        # Experimental NOE value
        if not all(v == 0 for v in NOE_df["NMR exp"]):
            ax.plot(
                [str(x) for x in NOE_df.index],
                NOE_df["NMR exp"],
                linewidth=0,
                marker="x",
                label="NMR distance",
            )

        # Experimental NOE upper bound
        if not all(v == 0 for v in NOE_df["upper bound"]):
            d = np.zeros(len(NOE_df["upper bound"])) + 10
            ax.fill_between(
                [str(x) for x in NOE_df.index],
                d,
                NOE_df["upper bound"],
                color="lightgrey",
                alpha=0.5,
            )
            ax.plot(
                [str(x) for x in NOE_df.index],
                NOE_df["upper bound"],
                label="NMR upper bound",
            )

        # Experimental NOE lower bound
        if not all(v == 0 for v in NOE_df["lower bound"]):
            e = np.zeros(len(NOE_df["lower bound"]))
            ax.fill_between(
                [str(x) for x in NOE_df.index],
                e,
                NOE_df["lower bound"],
                color="lightgrey",
                alpha=0.5,
            )
            ax.plot(
                [str(x) for x in NOE_df.index],
                NOE_df["lower bound"],
                label="NMR lower bound",
            )

        ax.set_xlabel("NOE no.")
        ax.set_ylabel(r"distance / $\AA$")

    ax.legend(
        mode="expand",
        ncol=2,
        bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
        loc="lower left",
        borderaxespad=0.0,
    )
    return fig, ax


def read_NOE(filepath):
    """Read NOE json and return a (tuple of) pandas.DataFrame object(s)

    Parameters
    ----------
    filepath : string
        NOE json filepath. The outermost dict. must contain a key
        and the NOE-dict as a value.
        e.g. {'24a': {NOE-dict for 24a}, '24b': {NOE-dict for 24b}, ...}

    Returns
    -------
    (pandas.DataFrame,)
        NOE dataframe(s)
    """
    import json
    import pandas as pd

    with open(filepath) as f:
        data = json.load(f)

    keys = data.keys()

    if len(keys) > 1:  # if contains multiple:
        li = []
        for k in keys:
            NOE_df = pd.DataFrame.from_dict(data[k], orient="index")
            li.append(NOE_df)
        return tuple(li)
    else:
        for k in keys:
            NOE_df = pd.DataFrame.from_dict(data[k], orient="index")
        return NOE_df
