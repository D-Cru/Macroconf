def mol_with_atom_index(mol):
    """
    Draw a (rdkit.mol) molecule with atom indices. (Atom 0 has no label)
    """
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


class dotdict(dict):
    # source: https://stackoverflow.com/questions/2352181/how-to-use-a-dot-
    # to-access-members-of-dictionary
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def json_load(path):
    """Loads a json file located at 'path'"""
    import json

    with open(path, "r") as f:
        compound = json.load(f)
    return dotdict(compound)


def json_dump(path, object):
    """Write 'object' to json file, given by 'path'"""
    import json

    with open(path, "w") as f:
        json.dump(object, f)


def pickle_dump(path, object):
    """Write 'object' to pickle file, given by 'path'"""
    import pickle

    with open(path, "wb") as f:
        pickle.dump(object, f)


def pickle_load(path):
    """Load 'object' from pickle file, given by 'path'"""
    import pickle

    with open(path, "rb") as f:
        object = pickle.load(f)
    return object


def link_ngl_wdgt_to_ax_pos(ax, pos, ngl_widget):
    # source: https://github.com/gph82/nglview_notebooks
    from matplotlib.widgets import AxesWidget
    from scipy.spatial import cKDTree

    """
    Link an ngl widget to a matplotlib axis object. Clicking on the matplotlib
    axis changes the frame in the ngl widget.
    Initial idea for this function comes from @arose,
    the rest is @gph82 and @clonker
    """

    kdtree = cKDTree(pos)
    assert ngl_widget.trajectory_0.n_frames == pos.shape[0]
    x, y = pos.T

    lineh = ax.axhline(ax.get_ybound()[0], c="black", ls="--")
    linev = ax.axvline(ax.get_xbound()[0], c="black", ls="--")
    (dot,) = ax.plot(pos[0, 0], pos[0, 1], "o", c="red", ms=7)

    ngl_widget.isClick = False

    def onclick(event):
        linev.set_xdata((event.xdata, event.xdata))
        lineh.set_ydata((event.ydata, event.ydata))
        data = [event.xdata, event.ydata]
        _, index = kdtree.query(x=data, k=1)
        dot.set_xdata((x[index]))
        dot.set_ydata((y[index]))
        ngl_widget.isClick = True
        ngl_widget.frame = index

    def my_observer(change):
        r"""Here comes the code that you want to execute"""
        ngl_widget.isClick = False
        _idx = change["new"]
        try:
            dot.set_xdata((x[_idx]))
            dot.set_ydata((y[_idx]))
        except IndexError:
            dot.set_xdata((x[0]))
            dot.set_ydata((y[0]))
            print(
                "caught index error with index %s (new=%s, old=%s)"
                % (_idx, change["new"], change["old"])
            )

    # Connect axes to widget
    axes_widget = AxesWidget(ax)
    axes_widget.connect_event("button_release_event", onclick)

    # Connect widget to axes
    ngl_widget.observe(my_observer, "frame", "change")


def color_cycle():
    """Cycle through colors in Bokeh."""
    # source: https://stackoverflow.com/questions/39839409/when-plotting-with-
    # bokeh-how-do-you-automatically-cycle-through-a-color-pallett
    # select a palette
    from bokeh.palettes import Dark2_5 as palette
    import itertools

    # itertools handles the cycling
    # create a color iterator
    colors = itertools.cycle(palette)
    return colors


def pymol_image(
    mol, ref, mol_state=1, ref_state=1, plotHbonds=True, label=False
):
    """Plot a pymol image of a CP molecule and its reference.

    Parameters
    ----------
    mol : filepath
        File path to the CP molecule.
    ref : filepath
        File path to the CP reference molecule.
    mol_state : int, optional
        pymol state (which structure to consider if multiple), by default 1
    ref_state : int, optional
        pymol state (which structure to consider if multiple), by default 1
    plotHbonds : bool, optional
        Show H bond contacts, by default True

    Returns
    -------
    mpimg.imread image
        Returns an image of the pymol plot, via mpimg.imread(). This can then
        be used in matplotlib plots, e.g. as part of a subplot.
    """
    from pymol import cmd
    import mdtraj as md
    import matplotlib.image as mpimg
    import tempfile

    # check if mol and ref are mdtraj trajectories
    if isinstance(mol, md.Trajectory):
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            mol.save_pdb(tmp.name)
            mol = tmp.name

    if isinstance(ref, md.Trajectory):
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            ref.save_pdb(tmp.name)
            ref = tmp.name

    if ref is None:
        ref = mol
    cmd.reinitialize()
    cmd.load(mol, object="mol")
    cmd.load(ref, object="ref")
    cmd.hide("everything", "mol")
    cmd.show("sticks", "mol")
    cmd.set("state", str(mol_state))
    # Align mol to ref
    cmd.align(
        "polymer and name CA and (mol)",
        "polymer and name CA and (ref)",
        object="aln_mol_to_ref",
        reset=1,
    )
    cmd.disable("aln_mol_to_ref")
    if plotHbonds:
        # Plot H bond contacts
        cmd.dist(
            "mol_polar_contacts", "mol", "mol", mode=2, label=0, reset=1
        )  # ;cmd.enable("mol_polar_contacts")
    cmd.disable("ref")
    cmd.enable("mol")
    cmd.center("ref", animate=-1)
    # Orient view relative to reference, this ensures standardised view
    cmd.orient("polymer and backbone and (ref)", state=ref_state)
    cmd.zoom("mol", 1, state=str(mol_state))
    if label:
        cmd.set("label_size", 40)
        cmd.set("label_color", "black")
        cmd.set("label_bg_color", "white")
        cmd.set("label_bg_transparency", 0.7)
        cmd.set("label_bg_outline", 0)
        cmd.set("label_shadow_mode", 0)
        cmd.set("label_position", [0, 0, 5])
        # Label CA atoms with residue number
        # cmd.label("polymer and name CA and (mol)", "'%s' % resn")
        cmd.label("mol and name CA", expression="resn")  # resn
        # cmd.label('name CA','resn')  #  and n. CA  and (byres((all)))
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_color", "black")
    cmd.set("ray_trace_mode", 3)
    cmd.set("antialias", 2)
    cmd.ray(1000, 1000)
    with tempfile.NamedTemporaryFile(suffix=".png") as tmp:
        cmd.save(tmp.name)
        image = mpimg.imread(tmp.name)

    return image


def determine_no_plots(n, x=4, y=5):
    """Determine no of rows, columns and figsize for a subplot plot.

    Parameters
    ----------
    n : int
        required no. of plots
    x : int, optional
        X aspect ratio (column width), by default 4
    y : int, optional
        Y aspect ratio (row height), by default 5

    Returns
    -------
    (int, int, (float, float))
        Tuple of (rows, cols, figsize).
    """
    if n <= 3:
        rows = 1
        cols = n
    else:
        cols = 3
        if n % 3 == 0:
            rows = n // 3
        else:
            rows = (n // 3) + 1

    figsize = (x * cols, y * rows)
    return (rows, cols, figsize)


def generate_symlink(source, dst):
    import os

    try:
        os.symlink(source, dst)
    except FileExistsError:
        os.remove(dst)
        os.symlink(source, dst)


def remove_symlinks(directory):
    import os

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)


def show_svg(filename):
    """Display an SVG in a Jupyter notebook."""
    from IPython.display import SVG, display

    return display(SVG(filename))
