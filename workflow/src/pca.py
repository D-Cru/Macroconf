def make_PCA(traj, type, dimensions=None):
    """Perform Principal Component Analysis (PCA) with different input types

    Automatically perform Principal Component Analysis (PCA) with various
    different input types.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        mdtraj trajectory object
    type : string
        input type. one of
          "cartesian": x,y,z coordinates of all atoms
          "dihedral": reduced dihedral angles
            (src.dihedrals.getReducedDihedrals)
          "pairwise_d": pairwise distances of all atoms
          "pairwise_NO": pairwise distances between all N, O atoms
            (N-N, N-O, O-O)
          "pairwise_N_O": pairwise distances of all N-O distances
          "dihe_pNO": reduced dihedral angles and pairwise distances of all
            N-O distances
    dimensions : int, optional
        Principal components to project on, by default 2

    Returns
    -------
    sklearn.decomposition.PCA, np.array
        Returns sklearn PCA-object and the transformed output
    """
    from sklearn.decomposition import PCA
    import mdtraj as md
    from src.dihedrals import getReducedDihedrals

    if dimensions is None:
        dimensions = 2

    pca = PCA(n_components=dimensions)

    # Cartesian
    if type == "cartesian":
        pca_input = traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)

    elif type == "dihedral":
        pca_input = getReducedDihedrals(traj)

    elif type == "pairwise_d":
        # This python function gives all unique pairs of elements from a list
        from itertools import combinations

        atom_pairs = list(combinations(range(traj.n_atoms), 2))
        pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
        pca_input = pairwise_distances

    elif type == "pairwise_NO":
        # Between all N & O atoms -> (N-N, N-O, O-O)
        from itertools import combinations

        N_O = traj.top.select("name N O")
        atom_pairs = list(combinations(N_O, 2))
        pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
        pca_input = pairwise_distances

    elif type == "pairwise_N_O":
        # Between all N-O distances -> (N-O only)
        from itertools import combinations

        N = traj.top.select("name N")
        Ox = traj.top.select("name O")
        atom_pairs = []
        for n in N:
            for o in Ox:
                atom_pairs.append([n, o])
        pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
        pca_input = pairwise_distances

    elif type == "dihe_pNO":
        # Dihedral angles and pairwise distances of all N-O distances
        import numpy as np

        N = traj.top.select("name N")
        Ox = traj.top.select("name O")
        atom_pairs = []
        for n in N:
            for o in Ox:
                atom_pairs.append([n, o])
        pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
        dihedrals = getReducedDihedrals(traj)
        pca_input = np.concatenate((pairwise_distances, dihedrals), axis=1)

    reduced_output = pca.fit_transform(pca_input)
    return pca, reduced_output


def plot_PCA(
    reduced_variables,
    type,
    compound,
    weight=None,
    cbar_label=None,
    fig=None,
    ax=None,
    cbar_plot=None,
    explained_variance=None,
):
    """Produce a PCA plot (2d scatter plot)

    Parameters
    ----------
    reduced_variables : np.array
        variables in PCA space (transformed)
    type : string
        type of PCA performed, for title of plot. Options are:
        ['cartesian', 'dihedral', 'pairwise']
    compound : string
        Compound name
    weight : np.array, optional
        weight vector to color scatterplot, by default None
    cbar_label : string, optional
        Colorbar legend, by default None
    fig : matplotlib.figure, optional
        can supply an already existing figure & axis to plot into.
            If only 1 of the two, or None is supplied, produce fig, ax via
            plt.subplots() from scratch., by default None
    ax : matplotlib.axes, optional
        can supply an already existing figure & axis to plot into.
            If only 1 of the two, or None is supplied, produce fig, ax via
            plt.subplots() from scratch., by default None
    cbar_plot : None, optional
        set to anything but None to disable colorbar, by default None
    explained_variance : array, optional
        array([0.0, 0.0]). If supplied, adds proportion of explained variance
        to the x,y axes, by default None

    Returns
    -------
    matplotlib.axes
        matplotlib axis object of a PCA plot.
    """
    import matplotlib.pyplot as plt

    plot_titles = {
        "cartesian": "Cartesian coordinate",
        "dihedral": "Dihedral",
        "pairwise": "Pairwise distance",
        "CP": "Cremer-Pople",
    }

    if (fig is None) or (ax is None):
        fig, ax = plt.subplots(figsize=(3.2677, 3.2677))
    im = ax.scatter(
        reduced_variables[:, 0],
        reduced_variables[:, 1],
        marker=".",
        s=0.5,
        c=weight,
        cmap="Spectral_r",
        vmin=0,
        vmax=8,
        rasterized=True,
    )
    if explained_variance is None:
        explained_variance = [0.0, 0.0]
    ax.set_xlabel(f"Principal Component 1 ({explained_variance[0]:.2f})")
    ax.set_ylabel(f"Principal Component 2 ({explained_variance[1]:.2f})")
    ax.set_title(f"{plot_titles[type]} PCA: compound {compound}")
    if (weight is not None) and (cbar_plot is None):
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(cbar_label)
    return ax


def plot_PCA_citra(c1, c2, type, compound, label=None, fig=None, ax=None):
    """Produce a PCA plot for a cis/trans compound.
    TODO: add reweighting / colorbar
    Parameters
    ----------
    c1 : np.array
        variables in PCA space (transformed) for 1/2
    c2 : np.array
        variables in PCA space (transformed) for 2/2
    type : string
        type of PCA performed, for title of plot. Options are:
        ['cartesian', 'dihedral', 'pairwise']
    compound : string
        Compound name
    label : array([string, string]), optional
        Labels to label plots, e.g. a,b or cis/trans, by default None
    fig : matplotlib.figure, optional
        can supply an already existing figure & axis to plot into.
            If only 1 of the two, or None is supplied, produce fig, ax via
            plt.subplots() from scratch., by default None
    ax : matplotlib.axes, optional
        can supply an already existing figure & axis to plot into.
            If only 1 of the two, or None is supplied, produce fig, ax via
            plt.subplots() from scratch., by default None

    Returns
    -------
    matplotlib.axes
        Returns a pca plot for a cis trans compound
    """
    import matplotlib.pyplot as plt

    plot_titles = {
        "cartesian": "Cartesian coordinate",
        "dihedral": "Dihedral",
        "pairwise": "Pairwise distance",
        "CP": "Cremer-Pople",
    }

    if label is None:
        label = ["Compound 1", "Compound 2"]

    if (fig is None) or (ax is None):
        fig, ax = plt.subplots()
    ax.scatter(
        c1[:, 0], c1[:, 1], marker=".", s=0.5, label=label[0], rasterized=True
    )
    ax.scatter(
        c2[:, 0], c2[:, 1], marker=".", s=0.5, label=label[1], rasterized=True
    )
    ax.legend()
    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_title(f"{plot_titles[type]} PCA: compound {compound}")
    return ax
