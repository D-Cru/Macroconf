import numpy as np
from sklearn.metrics import pairwise_distances


def getDihedrals(t):
    """Compute dihedral angles Phi, Psi, Omega for a cyclic peptide
    This is required because the MDTraj library does not support cyclic
    peptides, when computing dihedral angles.

    Parameters
    ----------
    t : mdtraj.Trajectory
        MDTraj trajectory

    Returns
    -------
    (numpy.ndarray, numpy.ndarray, numpy.ndarray)
        tuple of dihedral angles Phi, Psi, Omega
    """
    import mdtraj as md

    phi_idx = []
    psi_idx = []
    omega_idx = []

    # Determine atom indices
    for i in range(0, t.n_residues):
        n_ca_c = t.top.select(f"resid {i} and name N CA C")

        if i - 1 < 0:
            i_tmp = t.n_residues - 1
            c_pre = t.top.select(f"resid {i_tmp} and name C")
        else:
            c_pre = t.top.select(f"resid {i-1} and name C")

        if i + 1 >= t.n_residues:
            i_tmp = 0
            n_next = t.top.select(f"resid {i_tmp} and name N")
        else:
            n_next = t.top.select(f"resid {i+1} and name N")

        phi_i = np.append(c_pre, n_ca_c)
        psi_i = np.append(n_ca_c, n_next)

        phi_idx.append(phi_i)
        psi_idx.append(psi_i)

        ca_c = t.top.select(f"resid {i} and name CA C")

        if i + 1 >= t.n_residues:
            i_tmp = 0
            n_ca_next = t.top.select(f"resid {i_tmp} and name N CA")
        else:
            n_ca_next = t.top.select(f"resid {i+1} and name N CA")

        omega_i = np.append(ca_c, n_ca_next)
        omega_idx.append(omega_i)

    # compute angles
    phi = md.compute_dihedrals(t, phi_idx)
    psi = md.compute_dihedrals(t, psi_idx)
    omega = md.compute_dihedrals(t, omega_idx)
    return phi, psi, omega


def getReducedDihedrals(t):
    """Compute reduced dihedral angles Phi, Psi, Omega for a cyclic peptide
    Reduced dihedral angles are the sin/cos of the dihedral angles.

    Parameters
    ----------
    t : mdtraj.Trajectory
        MDTraj trajectory object

    Returns
    -------
    np.hstack
        sin/cos of the dihedral angles as hstack:
        cos(phi)|sin(phi)|cos(psi)|sin(psi)|cos(omega)|sin(omega)
    """
    # compute dihedral angles
    phi, psi, omega = getDihedrals(t)

    # transform angles
    d1 = np.cos(phi)
    d2 = np.sin(phi)
    d3 = np.cos(psi)
    d4 = np.sin(psi)
    d5 = np.cos(omega)
    d6 = np.sin(omega)

    reduced_dihedrals = np.hstack((d1, d2, d3, d4, d5, d6))
    return reduced_dihedrals


def angle_mean(theta):
    """Computes the degree average of an angle

    Computes the degree average of an angle (taking into account periodicity).
    For details, see here: https://en.wikipedia.org/wiki/Circular_mean

    Parameters
    ----------
    theta : numpy.array
        Angle in radians
    axis : int, optional
        _description_, by default 0

    Returns
    -------
    float
        mean angle in radians
    """
    tmp1 = np.mean(np.sin(theta), axis=0)
    tmp2 = np.mean(np.cos(theta), axis=0)

    return np.arctan2(tmp1, tmp2)


def miao_ramachandran(phi, psi):
    """Compute the structural motive according to Miao et al 2021.
    Experimental. For details on structural motives, see:
    https://doi.org/10.1039/D1SC05562C

    Parameters
    ----------
    phi : numpy.array
        Array of phi angles
    psi : numpy.array
        Array of phi angles

    Returns
    -------
    _type_
        _description_
    """
    # Derived from Miao et al 2021, via Inkscape (360x360px image, then hover
    # over cluster centroids.) from supplement figure
    centroids = np.array(
        [
            [-115.1, 147.6],
            [-61.2, 118.8],
            [-75.6, 54.0],
            [-67.7, -28.8],
            [-101.5, -132.4],
            [104.4, 129.6],
            [72.0, 25.2],
            [79.2, -57.6],
            [64.8, -121.6],
            [118.8, -150.4],
        ]
    )

    # Compute euclidean distance, account for periodicity of angles
    # src: https://stackoverflow.com/a/10405273
    dimension = [360, 360]

    def periodic_distance(p1, p2):
        total = 0
        for i, (a, b) in enumerate(zip(p1, p2)):
            delta = abs(b - a)
            if delta > dimension[i] - delta:
                delta = dimension[i] - delta
            total += delta**2
        return total**0.5

    if len(phi) == 1:
        X = np.array([phi, psi]).reshape(1, -1)
    elif len(phi) > 1:
        X = np.column_stack((phi, psi))

    # compute distances to centroids
    # src http://ethen8181.github.io/machine-learning/
    # clustering/kmeans.html#K-means
    distances_to_centroids = pairwise_distances(
        X, centroids, metric=periodic_distance
    )
    cluster_assignment = np.argmin(distances_to_centroids, axis=1)

    cluster_names = {
        "0": "B",
        "1": "Π",
        "2": "Γ",
        "3": "Λ",
        "4": "Z",
        "9": "β",
        "8": "π",
        "7": "γ",
        "6": "λ",
        "5": "ζ",
    }
    motive = []
    for i in cluster_assignment:
        motive.append(cluster_names[f"{i}"])

    return motive
