# TODO: Finish documentation for this!
def pyreweight_2d(input, args, weight_data=None):
    """adapted from Miao et al.
    input: 2D np.array of coordinates (e.g. reduced dihedrals from PCA)
    weight: file pointer to weights file
    job: jobtypes are amdweight_MC, amdweight_CE, noweight and others

    possible arguments are:
    args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff, T, Emax, fit, order}

    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy
    import scipy.special

    def loadfiletoarray(file):
        loaded = np.loadtxt(file, usecols=[0, 1])
        print("DATA LOADED:    " + file)
        return loaded

    def weightparse(rows, args):
        if args.job == "weighthist":
            data = np.loadtxt(args.weight)
            weights = data[:, 0]
            dV = np.zeros(rows)
        elif (
            args.job == "amd_time"
            or args.job == "amd_dV"
            or args.job == "amdweight"
            or args.job == "amdweight_MC"
            or args.job == "amdweight_CE"
        ):
            if weight_data is None:
                data = np.loadtxt(args.weight)
            else:
                data = weight_data
            weights = np.exp(data[:, 0])
            dV = data[:, 2]
        elif args.job == "noweight" or args.job == "histo":
            weights = np.zeros(rows)
            weights = weights + 1
            dV = np.zeros(rows)
        else:
            print("ERROR JOBTYPE" + args.job + " NOT RECOGNIZED")
            del data
            del weights
            del dV
        return weights, dV

    def histo(data, hist_min, binsX, discX, binsY, discY):
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=None
        )
        return hist2, newedgesX, newedgesY

    def assignbins(dim, disc):
        minimum = float(dim[0])
        maximum = float(dim[1])
        bins = np.arange(minimum, (maximum + disc), disc)
        return bins

    def normalize2D(pmf, cb_max):
        pmf = pmf - np.min(pmf)  # zero value to lowest energy state
        temphist = pmf
        # set infinity free energy values to is cb_max
        for jy in range(len(temphist[0, :])):
            for jx in range(len(temphist[:, 0])):
                if np.isinf(temphist[jx, jy]):
                    temphist[jx, jy] = cb_max
        return temphist

    def prephist(hist2, T, cb_max):
        hist2 = np.add(hist2, 0.000000000000000001)  # so that distrib
        hist2 = (0.001987 * T) * np.log(
            hist2
        )  # Convert to free energy in Kcal/mol
        hist2 = np.max(hist2) - hist2  # zero value to lowest energy state
        temphist2 = hist2
        # set infinity free energy values to is cb_max
        for jy in range(len(temphist2[0, :])):
            for jx in range(len(temphist2[:, 0])):
                if np.isinf(temphist2[jx, jy]):
                    temphist2[jx, jy] = cb_max
        return temphist2

    # memory usage is much reduced with multidimensional list for dV_mat;
    # pretty fast ~ O(N)
    def reweight_CE(data, hist_min, binsX, discX, binsY, discY, dV, T, fit):
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=None
        )

        beta = 1.0 / (0.001987 * T)
        nf = len(data[:, 0])
        nbinsX = len(hist2[:, 0])
        nbinsY = len(hist2[0, :])

        c1 = np.zeros((nbinsX, nbinsY))
        c2 = np.zeros((nbinsX, nbinsY))
        c3 = np.zeros((nbinsX, nbinsY))

        binfX = np.zeros(nf)  # array for storing assigned bin of each frame
        binfY = np.zeros(nf)  # array for storing assigned bin of each frame
        nA = np.zeros(
            (nbinsX, nbinsY), dtype=np.int
        )  # nA is equivalent to hist here
        dV_avg = np.zeros((nbinsX, nbinsY))
        dV_avg2 = np.zeros((nbinsX, nbinsY))
        dV_avg3 = np.zeros((nbinsX, nbinsY))
        dV_std = np.zeros((nbinsX, nbinsY))

        dV_avg_all = np.average(dV)
        dV_std_all = np.std(dV)
        print("dV all: avg = ", dV_avg_all, "std = ", dV_std_all)

        dV_mat = [
            [[[] for i in range(1)] for i in range(nbinsY)]
            for i in range(nbinsX)
        ]
        for i in range(len(data[:, 0])):
            jx = int((data[i, 0] - binsX[0]) / discX)
            jy = int((data[i, 1] - binsY[0]) / discY)
            if jx < nbinsX and jy < nbinsY:
                binfX[i] = jx
                binfY[i] = jy
                dV_mat[jx][jy].append(dV[i])
                nA[jx, jy] = nA[jx, jy] + 1

        for jx in range(nbinsX):
            for jy in range(nbinsY):
                if nA[jx, jy] >= hist_min:
                    num = int(nA[jx, jy])
                    atemp = np.asarray(dV_mat[jx][jy][1 : num + 1])
                    atemp2 = np.power(atemp, 2)
                    atemp3 = np.power(atemp, 3)
                    dV_avg[jx, jy] = np.average(atemp)
                    dV_std[jx, jy] = np.std(atemp)
                    dV_avg2[jx, jy] = np.average(atemp2)
                    dV_avg3[jx, jy] = np.average(atemp3)
                    del atemp
                    del atemp2
                    del atemp3
                    c1[jx, jy] = beta * dV_avg[jx, jy]
                    c2[jx, jy] = 0.5 * beta**2 * dV_std[jx, jy] ** 2
                    c3[jx, jy] = (
                        (1.0 / 6.0)
                        * beta**3
                        * (
                            dV_avg3[jx, jy]
                            - 3.0 * dV_avg2[jx, jy] * dV_avg[jx, jy]
                            + 2.0 * dV_avg[jx, jy] ** 3
                        )
                    )

        del dV_mat
        del dV_avg
        del dV_avg2
        del dV_avg3
        del dV_std
        return hist2, newedgesX, newedgesY, c1, c2, c3

    def reweight_dV(data, hist_min, binsX, binsY, discX, discY, dV, T):
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=None
        )

        nf = len(data[:, 0])
        nbinsX = len(hist2[:, 0])
        nbinsY = len(hist2[0, :])

        binfX = np.zeros(nf)  # array for storing assigned bin of each frame
        binfY = np.zeros(nf)  # array for storing assigned bin of each frame
        nA = np.zeros(
            (nbinsX, nbinsY), dtype=np.int
        )  # nA is equivalent to hist here
        dV_avg = np.zeros((nbinsX, nbinsY))
        dV_std = np.zeros((nbinsX, nbinsY))
        dV_anharm = np.zeros((nbinsX, nbinsY))

        dV_mat = [
            [[[] for i in range(1)] for i in range(nbinsY)]
            for i in range(nbinsX)
        ]
        for i in range(len(data[:, 0])):
            jx = int((data[i, 0] - binsX[0]) / discX)
            jy = int((data[i, 1] - binsY[0]) / discY)
            if jx < nbinsX and jy < nbinsY:
                binfX[i] = jx
                binfY[i] = jy
                dV_mat[jx][jy].append(dV[i])
                nA[jx, jy] = nA[jx, jy] + 1

        for jx in range(nbinsX):
            for jy in range(nbinsY):
                dV_anharm[jx, jy] = 100
                if nA[jx, jy] >= hist_min:
                    num = int(nA[jx, jy])
                    atemp = np.asarray(dV_mat[jx][jy][1 : num + 1])
                    dV_avg[jx, jy] = np.average(atemp)
                    dV_std[jx, jy] = np.std(atemp)
                    dV_anharm[jx, jy] = anharm(atemp)
                    del atemp
        return (
            hist2,
            newedgesX,
            newedgesY,
            binfX,
            binfY,
            dV_avg,
            dV_std,
            dV_anharm,
            dV_mat,
        )

    ##  Convert histogram to free energy in Kcal/mol
    def hist2pmf2D(hist, hist_min, T):
        nbinsX = len(hist[:, 0])
        nbinsY = len(hist[0, :])
        pmf = np.zeros((nbinsX, nbinsY))
        pmf_min = 100
        for jx in range(len(hist[:, 0])):
            for jy in range(len(hist[0, :])):
                if hist[jx, jy] >= hist_min:
                    pmf[jx, jy] = -(0.001987 * T) * np.log(hist[jx, jy])
                if pmf_min > pmf[jx, jy]:
                    pmf_min = pmf[jx, jy]
        ##        pmf=pmf-pmf_min  ## zero value to lowest energy state
        return pmf

    def output_pmf2D(pmffile, hist, binsX, binsY):
        fpmf = open(pmffile, "w")
        strpmf = '#RC1\tRC2\tPMF(kcal/mol)\n\n@    xaxis  label "RC1"\n@    yaxis  label "RC2"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(hist[:, 0])):
            for jy in range(len(hist[0, :])):
                strpmf = (
                    str(binsX[jx])
                    + " \t"
                    + str(binsY[jy])
                    + " \t"
                    + str(hist[jx, jy])
                    + "\n"
                )
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV(pmffile, dV):
        fpmf = open(pmffile, "w")
        strpmf = '#dV \tp(dV) \n\n@    xaxis  label "dV"\n@    yaxis  label "p(dV)"\n@TYPE xy\n'
        hist_dV, bin_dV = np.histogram(dV, bins=50)
        for k in range(len(hist_dV)):
            strpmf = strpmf + str(bin_dV[k]) + " \t" + str(hist_dV[k]) + " \n"
        fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_anharm2D(pmffile, binsX, binsY, dV_anharm):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tdV_anharm \tError\n\n@    xaxis  label "RC"\n@    yaxis  label "dV_anmarm"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(dV_anharm[:, 0])):
            for jy in range(len(dV_anharm[0, :])):
                strpmf = (
                    str(binsX[jx])
                    + " \t"
                    + str(binsY[jy])
                    + " \t"
                    + str(dV_anharm[jx, jy])
                    + "\n"
                )
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_stat2D(pmffile, binsX, binsY, dV_avg, dV_std, dV_anharm):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label "RC"\n@    yaxis  label "dV(kcal/mol)"\n@TYPE xydy\n'
        fpmf.write(strpmf)
        for jx in range(len(dV_anharm[:, 0])):
            for jy in range(len(dV_anharm[0, :])):
                strpmf = (
                    str(binsX[jx])
                    + " \t"
                    + str(binsY[jy])
                    + " \t"
                    + str(dV_avg[jx, jy])
                    + " \t"
                    + str(dV_std[jx, jy])
                    + " \t"
                    + str(dV_anharm[jx, jy])
                    + "\n"
                )
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_mat2D(
        pmffile, binsX, binsY, hist, dV_avg, dV_std, dV_anharm, dV_mat
    ):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tNf \tdV_avg \tdV_std \tdV_ij \n\n@    xaxis  label "RC"\n@    yaxis  label "dV(kcal/mol)"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for jx in range(len(hist[:, 0])):
            for jy in range(len(hist[0, :])):
                nf_j = int(hist[jx, jy])
                strpmf = (
                    str(binsX[jx])
                    + " \t"
                    + str(binsY[jy])
                    + " \t"
                    + str(hist[jx, jy])
                    + " \t"
                    + str(dV_avg[jx, jy])
                    + " \t"
                    + str(dV_std[jx, jy])
                    + " \t"
                    + str(dV_anharm[jx, jy])
                )
                strpmf = strpmf + " \t" + str(dV_mat[jx][jy][1 : nf_j + 1])
                strpmf = strpmf + "\n"
                fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def anharm(data):
        var = np.var(data)
        hist, edges = np.histogram(data, 50, normed=True)
        hist = np.add(hist, 0.000000000000000001)  ###so that distrib
        dx = edges[1] - edges[0]
        S1 = -1 * np.trapz(np.multiply(hist, np.log(hist)), dx=dx)
        S2 = 0.5 * np.log(
            2.00 * np.pi * np.exp(1.0) * var + 0.000000000000000001
        )
        alpha = S2 - S1
        if np.isinf(alpha):
            alpha = 100
        return alpha

    plt_figs = 0

    data = input

    rows = len(data[:, 0])
    weights, dV = weightparse(rows, args)

    if args.discX:
        discX = float(args.discX)
    else:
        discX = 6

    if args.discY:
        discY = float(args.discY)
    else:
        discY = 6

    if args.Xdim:
        binsX = assignbins(args.Xdim, discX)
    else:
        max_data = discX * (int(np.amax(data[:, 0]) / discX) + 1)
        min_data = discX * (int(np.amin(data[:, 0]) / discX) - 1)
        binsX = assignbins([min_data, max_data], discX)  ## Default bin size

    if args.Ydim:
        binsY = assignbins(args.Ydim, discY)
    else:
        max_data = discY * (int(np.amax(data[:, 1]) / discY) + 1)
        min_data = discY * (int(np.amin(data[:, 1]) / discY) - 1)
        binsY = assignbins([min_data, max_data], discY)  ## Default bin size

    ##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        cb_max = float(args.Emax)
    else:
        cb_max = 8

    ##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min = int(args.cutoff)
    else:
        hist_min = 10  # minimum number of configurations in one bin

    ##  SET ORDER of McLaurin series expansion
    if args.order:
        order = int(args.order)
    else:
        order = 10  # default

    ##  SET TEMPERATURE
    if args.T:
        T = float(args.T)
    else:
        T = 300  # simulation temperature
    beta = 1.0 / (0.001987 * T)

    ##REWEIGHTING
    ##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit = args.fit
    else:
        fit = False  # simulation temperature

    ##REWEIGHTING
    if args.job == "amdweight_CE":
        hist2, newedgesX, newedgesY, c1, c2, c3 = reweight_CE(
            data, hist_min, binsX, discX, binsY, discY, dV, T, fit
        )
        pmf = hist2pmf2D(hist2, hist_min, T)
        c1 = -np.multiply(1.0 / beta, c1)
        c2 = -np.multiply(1.0 / beta, c2)
        c3 = -np.multiply(1.0 / beta, c3)

        c12 = np.add(c1, c2)
        c123 = np.add(c12, c3)
        pmf_c1 = np.add(pmf, c1)
        print("pmf_min-c1 = ", np.min(pmf_c1))
        pmf_c1 = normalize2D(pmf_c1, cb_max)
        pmf_c2 = np.add(pmf, c12)
        print("pmf_min-c2 = ", np.min(pmf_c2))
        pmf_c2 = normalize2D(pmf_c2, cb_max)
        pmf_c3 = np.add(pmf, c123)
        print("pmf_min-c3 = ", np.min(pmf_c3))
        pmf_c3 = normalize2D(pmf_c3, cb_max)
    elif args.job == "amdweight_MC":
        n = order
        MCweight = np.zeros(len(dV))
        beta_dV = np.multiply(dV, beta)
        for x in range(0, n + 1):
            MCweight = np.add(
                MCweight,
                (
                    np.divide(
                        np.power(beta_dV, x), float(scipy.special.factorial(x))
                    )
                ),
            )  # changed from misc to special!
        weights = MCweight
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=weights
        )
        hist2 = prephist(hist2, T, cb_max)
    elif args.job == "amdweight":
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=weights
        )
        hist2 = prephist(hist2, T, cb_max)
    else:
        hist2, newedgesX, newedgesY = np.histogram2d(
            data[:, 0], data[:, 1], bins=(binsX, binsY), weights=None
        )
        hist2 = prephist(hist2, T, cb_max)

    ##SAVE FREE ENERGY DATA INTO A FILE
    if (
        args.job == "amdweight_MC"
        or args.job == "amdweight"
        or args.job == "noweight"
    ):
        pmffile = "pmf-" + str(args.input) + ".xvg"
        return hist2, binsX, binsY
        # output_pmf2D(pmffile,hist2,binsX,binsY)
    if args.job == "amdweight_CE":
        hist2 = pmf_c1
        pmffile = "pmf-c1-" + str(args.input) + ".xvg"
        # output_pmf2D(pmffile,hist2,binsX,binsY)

        hist2 = pmf_c3
        pmffile = "pmf-c3-" + str(args.input) + ".xvg"
        # output_pmf2D(pmffile,hist2,binsX,binsY)

        hist2 = pmf_c2
        pmffile = "pmf-c2-" + str(args.input) + ".xvg"
        # output_pmf2D(pmffile,hist2,binsX,binsY)
        return hist2, binsX, binsY
    if args.job == "histo":
        hist2, newedgesX, newedgesY = histo(
            data, hist_min, binsX, discX, binsY
        )
        pmffile = "histo-" + str(args.input) + ".xvg"
        # output_dV_anharm2D(pmffile,binsX,binsY,hist2)

    if args.job == "amd_dV":
        plt_figs = 0
        (
            hist2,
            newedgesX,
            newedgesY,
            binfX,
            binfY,
            dV_avg,
            dV_std,
            dV_anharm,
            dV_mat,
        ) = reweight_dV(data, hist_min, binsX, binsY, discX, discY, dV, T)

        pmffile = "dV-hist-2D-" + str(args.input) + ".xvg"
        # output_dV(pmffile,dV)

        alpha = anharm(dV)
        print("Anharmonicity of all dV = " + str(alpha))

        pmffile = "dV-anharm-2D-" + str(args.input) + ".xvg"
        # output_dV_anharm2D(pmffile,binsX,binsY,dV_anharm)

        pmffile = "dV-stat-2D-" + str(args.input) + ".xvg"
        # output_dV_stat2D(pmffile,binsX,binsY,dV_avg,dV_std,dV_anharm)

        pmffile = "dV-mat-2D-" + str(args.input) + ".xvg"
        # output_dV_mat2D(pmffile,binsX,binsY,hist2,dV_avg,dV_std,dV_anharm,dV_mat)

    ###PLOTTING FUNCTION FOR FREE ENERGY FIGURE
    if plt_figs:
        cbar_ticks = [0, cb_max * 0.25, cb_max * 0.5, cb_max * 0.75, 8.0]
        plt.figure(2, figsize=(11, 8.5))
        extent = [newedgesX[0], newedgesX[-1], newedgesY[-1], newedgesY[0]]
        print(extent)
        plt.imshow(hist2.transpose(), extent=extent, interpolation="gaussian")
        cb = plt.colorbar(
            ticks=cbar_ticks, format=("% .1f"), aspect=10
        )  # grab the Colorbar instance
        imaxes = plt.gca()
        plt.sca(cb.ax)
        plt.clim(vmin=0, vmax=8.0)
        plt.yticks(fontsize=18)
        plt.sca(imaxes)
        axis = (min(binsX), max(binsX), min(binsY), max(binsY))
        plt.axis(axis)
        plt.xticks(size="18")
        plt.yticks(size="18")
        plt.xlabel("RC1", fontsize=18)
        plt.ylabel("RC2", fontsize=18)
        ##    	plt.xlabel(r'$\phi$',fontsize=18)
        ##    	plt.ylabel(r'$\psi$',fontsize=18)
        ##    	plt.xlabel(r'$\chi$1',fontsize=18)
        ##    	plt.ylabel(r'$\chi$2',fontsize=18)
        plt.savefig("2D_Free_energy_surface.png", bbox_inches=0)
        print("FIGURE SAVED 2D_Free_energy_surface.png")

        ###PLOTTING FUNCTION FOR WEIGHTS histogram
        [hist, edges] = np.histogram(weights, bins=100)
        width = np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11, 8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale(
            "log"
        )  ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize="18")
        plt.yticks(fontsize="18")
        plt.savefig("weights.png", bbox_inches=0)
        print("FIGURE SAVED weights.png")

    print(" ")
    print("END")


def reweight(reduced_dihedrals, weightf, jobtype, weight_data=None):
    """2D Reweighting wrapper.

    Calls pyreweight_2d, which was adapted from Miao et al
    Py-reweighting scripts.

    Parameters
    ----------
    reduced_dihedrals : np.array
        reduced dihedrals array, 2d quantity to reweigh.
    weightf : filepath
        weight-file
    jobtype : string
        Type of reweighting to perform. Possible options can be found in
        pyreweight_2d
    weight_data : np.array, optional
        array of weights, by default None

    Returns
    -------
    np.array
        weights for input array reduced dihedrals
    """

    from .utils import dotdict
    import numpy as np

    # Determine bounds and make 100 steps per dimension
    max_X = reduced_dihedrals[:, 0].max()
    min_X = reduced_dihedrals[:, 0].min()

    max_Y = reduced_dihedrals[:, 1].max()
    min_Y = reduced_dihedrals[:, 1].min()

    discX = (max_X - min_X) / 100
    discY = (max_Y - min_Y) / 100

    # now do reweighting...
    # possible arguments are:
    # args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff, T, Emax, fit, order}
    # -job amdweight_MC -weight weights.dat -Xdim -3 3 -discX 0.1 -Ydim -3 3 -discY 0.1 -Emax 20
    arguments = {
        "job": jobtype,
        "weight": weightf,
        "Xdim": [min_X, max_X],
        "discX": str(discX),
        "Ydim": [min_Y, max_Y],
        "discY": str(discY),
        "Emax": 10,
        "cutoff": 10,
    }  # , "cutoff":200} discX,Y: 0.02
    arguments = dotdict(arguments)
    hist2, binsX, binsY = pyreweight_2d(
        reduced_dihedrals, arguments, weight_data
    )

    # assign every input value the corresponding pmf.
    c1 = np.digitize(reduced_dihedrals[:, 0], binsX) - 1  # shape (50000,)
    c1[c1 == c1.max()] = c1.max() - 1
    c2 = np.digitize(reduced_dihedrals[:, 1], binsY) - 1  # shape (50000,)
    c2[c2 == c2.max()] = c2.max() - 1
    # c2 = [c2.max()-1 if i==c2.max() else i for i in c2]

    weights = hist2[c1, c2]
    # set weights over threshold to const. value
    # weights[weights > 15] = 15
    return weights


def pyreweight_1d(input, args, slicer=None, weight_data=None):
    import math
    import scipy
    import scipy.stats as stats
    import numpy as np
    import sys
    import matplotlib.pyplot as plt
    import csv
    from argparse import ArgumentParser
    from scipy.optimize import curve_fit

    ## from scipy.optimize import *

    def loadfiletoarray(file):
        loaded = np.loadtxt(file, usecols=[0])
        print("DATA LOADED:    " + file)
        return loaded

    def weightparse(rows, args, slice_mask=None):
        if args.job == "weighthist":
            data = np.loadtxt(args.weight)
            weights = data[:, 0]
            if slice_mask is not None:
                weights = weights[slice_mask]
            dV = np.zeros(rows)
        elif (
            args.job == "amd_time"
            or args.job == "amd_dV"
            or args.job == "amdweight"
            or args.job == "amdweight_MC"
            or args.job == "amdweight_CE"
        ):
            if weight_data is None:
                data = np.loadtxt(args.weight)
            else:
                data = weight_data
            weights = np.exp(data[:, 0])
            if slice_mask is not None:
                weights = weights[slice_mask]
            dV = data[:, 2]
            if slice_mask is not None:
                dV = dV[slice_mask]
        elif args.job == "noweight" or args.job == "histo":
            weights = np.zeros(rows)
            weights = weights + 1
            dV = np.zeros(rows)
        else:
            print("ERROR JOBTYPE" + args.job + " NOT RECOGNIZED")
            del data
            del weights
            del dV
        return weights, dV

    def reweight_CE(data, hist_min, binsX, discX, dV, T, fit):
        hist, newedgesX = np.histogram(data, bins=binsX, weights=None)

        beta = 1.0 / (0.001987 * T)
        nf = len(data)
        nbins = len(hist)

        c1 = np.zeros(nbins)
        c2 = np.zeros(nbins)
        c3 = np.zeros(nbins)

        binf = np.zeros(nf)  # array for storing assigned bin of each frame
        nA = np.zeros(nbins, dtype=np.int)  # nA is equivalent to hist here
        dV_avg = np.zeros(nbins)
        dV_avg2 = np.zeros(nbins)
        dV_avg3 = np.zeros(nbins)
        dV_std = np.zeros(nbins)

        dV_avg_all = np.average(dV)
        dV_std_all = np.std(dV)
        # print ('dV all: avg = ', dV_avg_all, 'std = ', dV_std_all)

        diff_tol_avg = dV_std_all * 3
        diff_tol_std = dV_std_all
        dV_binsize = 50

        dV_mat = [[[] for i in range(1)] for i in range(nbins)]
        for i in range(len(data)):
            j = int((data[i] - binsX[0]) / discX)
            if j == nbins:
                j = nbins - 1
            binf[i] = j
            dV_mat[j].append(dV[i])
            nA[j] = nA[j] + 1

        for j in range(nbins):
            if nA[j] >= hist_min:
                num = int(nA[j])
                atemp = np.asarray(dV_mat[j][1 : num + 1])
                atemp2 = np.power(atemp, 2)
                atemp3 = np.power(atemp, 3)
                if fit:
                    ## calculate average/std through gaussian fitting
                    hist_temp, bin_edges_temp = np.histogram(
                        atemp, bins=dV_binsize
                    )
                    bin_centres_temp = (
                        bin_edges_temp[:-1] + bin_edges_temp[1:]
                    ) / 2
                    ## output original histograms
                    pmffile = (
                        "dV-hist-forFit-RC"
                        + str("%#08.2f" % binsX[j])
                        + ".xvg"
                    )
                    output_pmf(pmffile, hist_temp, bin_centres_temp)
                    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
                    mean = np.average(atemp)
                    std = np.std(atemp)
                    p0 = [0.0, 1.0, mean, std]
                    ## coeff, var_matrix = curve_fit(gauss, bin_centres_temp, hist_temp, p0=p0)
                    coeff, var_matrix = curve_fit(
                        gauss, bin_centres_temp, hist_temp, p0=p0
                    )
                    # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
                    # print (binsX[j], ': mean = ', coeff[2], 'sigma = ', coeff[3])
                    dV_avg[j] = coeff[2]
                    dV_std[j] = coeff[3]
                    # Get the fitted curve
                    hist_fit = gauss(bin_centres_temp, *coeff)
                    ## output fitted histograms
                    pmffile = (
                        "dV-hist-gaussFit-RC"
                        + str("%#08.2f" % binsX[j])
                        + ".xvg"
                    )
                    output_pmf(pmffile, hist_fit, bin_centres_temp)
                else:
                    ## calculate average/std directly
                    dV_avg[j] = np.average(atemp)
                    dV_std[j] = np.std(atemp)

                ##	  if np.absolute(dV_avg[j]-dV_avg_all)>diff_tol_avg or np.absolute(dV_std[j]-dV_std_all)>diff_tol_std :
                ##	     dV_avg[j]=0
                ##	     dV_std[j]=0

                dV_avg2[j] = np.average(atemp2)
                dV_avg3[j] = np.average(atemp3)
                del atemp
                del atemp2
                del atemp3
                c1[j] = beta * dV_avg[j]
                c2[j] = 0.5 * beta**2 * dV_std[j] ** 2
                c3[j] = (
                    (1.0 / 6.0)
                    * beta**3
                    * (
                        dV_avg3[j]
                        - 3.0 * dV_avg2[j] * dV_avg[j]
                        + 2.0 * dV_avg[j] ** 3
                    )
                )

        del dV_mat
        del dV_avg
        del dV_avg2
        del dV_avg3
        del dV_std
        return hist, newedgesX, c1, c2, c3

    def reweight_dV(data, hist_min, binsX, discX, dV, T):
        hist, newedgesX = np.histogram(data, bins=binsX, weights=None)

        nf = len(data)
        nbins = len(hist)

        binf = np.zeros(nf)  # array for storing assigned bin of each frame
        nA = np.zeros(nbins, dtype=np.int)  # nA is equivalent to hist here
        dV_avg = np.zeros(nbins)
        dV_std = np.zeros(nbins)
        dV_anharm = np.zeros(nbins)

        dV_mat = [[[] for i in range(1)] for i in range(nbins)]
        for i in range(len(data)):
            j = int((data[i] - binsX[0]) / discX)
            if j >= nbins:
                j = nbins - 1
            binf[i] = j
            dV_mat[j].append(dV[i])
            nA[j] = nA[j] + 1

        for j in range(nbins):
            dV_anharm[j] = 100
            if nA[j] > 0:
                num = int(nA[j])
                atemp = np.asarray(dV_mat[j][1 : num + 1])
                dV_avg[j] = np.average(atemp)
                dV_std[j] = np.std(atemp)
                dV_anharm[j] = anharm(atemp)
                del atemp
        return hist, newedgesX, binf, dV_avg, dV_std, dV_anharm, dV_mat

    def histo(data, hist_min, binsX, discX):
        hist, newedgesX = np.histogram(data, bins=binsX, weights=None)
        return hist, newedgesX

    def assignbins(dim, disc):
        minimum = float(dim[0])
        maximum = float(dim[1])
        bins = np.arange(minimum, (maximum + disc), disc)
        return bins

    def normalize(pmf, cb_max):
        pmf = pmf - np.min(pmf)  ## zero value to lowest energy state
        temphist = pmf
        # set infinity free energy values to is cb_max
        for x in range(len(temphist[:])):
            if np.isinf(temphist[x]):
                temphist[x] = cb_max
        return temphist

    def prephist(hist, T, cb_max):
        hist = np.add(hist, 0.000000000000000001)  ###so that distrib
        hist = (0.001987 * T) * np.log(
            hist
        )  ####Convert to free energy in Kcal/mol
        # print ("PMF_min = ", -np.max(hist))
        hist = np.max(hist) - hist  ## zero value to lowest energy state
        temphist = hist
        # set infinity free energy values to is cb_max
        for x in range(len(temphist[:])):
            if np.isinf(temphist[x]):
                temphist[x] = cb_max
        return temphist

    def prepdV(hist, cb_max):
        for x in np.where(hist == 0)[0]:
            hist[x] = cb_max
        return hist

    ##  Convert histogram to free energy in Kcal/mol
    def hist2pmf(hist, hist_min, T):
        nbins = len(hist)
        pmf = np.zeros(nbins)
        pmf_min = 100
        for j in range(len(hist)):
            if hist[j] >= hist_min:
                pmf[j] = -(0.001987 * T) * np.log(hist[j])
                if pmf_min > pmf[j]:
                    pmf_min = pmf[j]
        return pmf

    def output_pmf(pmffile, hist, binsX):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tPMF(kcal/mol)\n\n@    xaxis  label "RC"\n@    yaxis  label "PMF(kcal/mol)"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for j in range(len(hist[:])):
            strpmf = str(binsX[j]) + " \t" + str(hist[j]) + "\n"
            fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV(pmffile, dV):
        fpmf = open(pmffile, "w")
        strpmf = '#dV \tp(dV) \n\n@    xaxis  label "dV"\n@    yaxis  label "p(dV)"\n@TYPE xy\n'
        hist_dV, bin_dV = np.histogram(dV, bins=50)
        for k in range(len(hist_dV)):
            strpmf = strpmf + str(bin_dV[k]) + " \t" + str(hist_dV[k]) + " \n"
        fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_anharm(pmffile, binsX, dV_anharm):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tdV_anharm \tError\n\n@    xaxis  label "RC"\n@    yaxis  label "dV_anmarm"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for j in range(len(dV_anharm[:])):
            strpmf = str(binsX[j]) + " \t" + str(dV_anharm[j]) + "\n"
            fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_stat(pmffile, binsX, dV_avg, dV_std, dV_anharm):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label "RC"\n@    yaxis  label "dV(kcal/mol)"\n@TYPE xydy\n'
        fpmf.write(strpmf)
        for j in range(len(dV_avg[:])):
            strpmf = (
                str(binsX[j])
                + " \t"
                + str(dV_avg[j])
                + " \t"
                + str(dV_std[j])
                + " \t"
                + str(dV_anharm[j])
                + "\n"
            )
            fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def output_dV_mat(pmffile, binsX, hist, dV_avg, dV_std, dV_anharm, dV_mat):
        fpmf = open(pmffile, "w")
        strpmf = '#RC \tNf \tdV_avg \tdV_std \tdV_anharm \tdV_ij \n\n@    xaxis  label "RC"\n@    yaxis  label "dV(kcal/mol)"\n@TYPE xy\n'
        fpmf.write(strpmf)
        for j in range(len(dV_avg[:])):
            nf_j = int(hist[j])
            strpmf = (
                str(binsX[j])
                + " \t"
                + str(hist[j])
                + " \t"
                + str(dV_avg[j])
                + " \t"
                + str(dV_std[j])
                + " \t"
                + str(dV_anharm[j])
            )
            strpmf = strpmf + " \t" + str(dV_mat[j][1 : nf_j + 1])
            strpmf = strpmf + "\n"
            fpmf.write(strpmf)
        fpmf.closed
        return fpmf

    def accel_amd_window(c12, itb, ite):
        avg = np.average(dV)
        std = np.std(dV)
        accel = np.exp(c12)
        return accel

    def accel_amd(dV, T):
        avg = np.average(dV)
        std = np.std(dV)
        accel = np.exp(
            avg / (0.001987 * T) + (std / (0.001987 * T)) ** 2 / 2.0
        )
        return accel

    def gauss(x, *p):
        y0, A, mu, sigma = p
        return y0 + A * np.exp(-((x - mu) ** 2) / (2.0 * sigma**2))

    def anharm(data):
        var = np.var(data)
        # hist, edges=np.histogram(data, 50, normed=True)
        hist, edges = np.histogram(data, 50, density=True)
        hist = np.add(hist, 0.000000000000000001)  ###so that distrib
        dx = edges[1] - edges[0]
        S1 = -1 * np.trapz(np.multiply(hist, np.log(hist)), dx=dx)
        S2 = 0.5 * np.log(
            np.add(2.00 * np.pi * np.exp(1) * var, 0.000000000000000001)
        )
        alpha = S2 - S1
        if np.isinf(alpha):
            alpha = 100
        return alpha

    def anharmND(datafull):
        print("Performing error estimation")
        width = datafull[0, :]
        lendata = len(datafull[:, 0])
        for x in range(len(width)):
            var = np.var(datafull[:, x])
            std = np.std(datafull[:, x])
            print(
                "variance of "
                + str(x)
                + " is : "
                + str(var)
                + " Standard Deviation:  "
                + str(std)
            )
            hist, edges = np.histogram(datafull[:, x], 100, normed=True)
            hist = np.add(hist, 0.000000000000000001)  ###so that distrib
            dx = edges[1] - edges[0]
            S1 = -1 * np.trapz(np.multiply(hist, np.log(hist)), dx=dx)
            S2 = 0.5 * np.log(
                np.add(2.00 * np.pi * np.exp(1) * var, 0.000000000000000001)
            )
            alpha = S2 - S1
            print(
                str(x)
                + "dimension   S1 = "
                + str(S1)
                + "  S2 = "
                + str(S2)
                + " Alpha = "
                + str(alpha)
            )
        return var, std, alpha

    ## Set control parameters
    plt_figs = 0

    data = input

    rows = len(data[:])
    weights, dV = weightparse(rows, args, slicer)
    if args.disc:
        discX = float(args.disc)
    else:
        discX = 6

    if args.Xdim:
        binsX = assignbins(args.Xdim, discX)
    else:
        max_data = discX * (int(np.amax(data) / discX) + 1)
        min_data = discX * (int(np.amin(data) / discX) - 1)
        binsX = assignbins([min_data, max_data], discX)  ## Default bin size

    ##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        cb_max = float(args.Emax)
    else:
        cb_max = 8

    ##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min = int(args.cutoff)
    else:
        hist_min = 10  # minimum number of configurations in one bin

    ##  SET ORDER of McLaurin series expansion
    if args.order:
        order = int(args.order)
    else:
        order = 10  # default

    ##  SET TEMPERATURE
    if args.T:
        T = float(args.T)
    else:
        T = 300  # simulation temperature
    beta = 1.0 / (0.001987 * T)

    ##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit = args.fit
    else:
        fit = False  # simulation temperature

    ##REWEIGHTING
    if args.job == "amdweight_CE":
        ##CALCULATE effective acceleration factor
        accel = accel_amd(dV, T)
        # print ("Effective acceleration factor with Gaussian approximation: ", accel)

        hist, newedgesX, c1, c2, c3 = reweight_CE(
            data, hist_min, binsX, discX, dV, T, fit
        )
        pmf = prephist(hist, T, cb_max)
        c1 = -np.multiply(1.0 / beta, c1)
        c2 = -np.multiply(1.0 / beta, c2)
        c3 = -np.multiply(1.0 / beta, c3)
        c12 = np.add(c1, c2)
        c123 = np.add(c12, c3)
        pmf_c1 = np.add(pmf, c1)
        # print ("pmf_min-c1 = ", np.min(pmf_c1))
        pmf_c1 = normalize(pmf_c1, cb_max)
        pmf_c2 = np.add(pmf, c12)
        # print ("pmf_min-c2 = ", np.min(pmf_c2))
        pmf_c2 = normalize(pmf_c2, cb_max)
        pmf_c3 = np.add(pmf, c123)
        # print ("pmf_min-c3 = ", np.min(pmf_c3))
        pmf_c3 = normalize(pmf_c3, cb_max)
    elif args.job == "amd_time":
        ##CALCULATE effective acceleration factor
        accel = accel_amd(dV, T)
        print(
            "Effective acceleration factor with Gaussian approximation: ",
            accel,
        )

        hist, newedgesX, c1, c2, c3 = reweight_CE(
            data, hist_min, binsX, discX, dV, T, fit
        )
        c1 = np.multiply(1.0 / beta, c1)
        c2 = np.multiply(1.0 / beta, c2)
        c3 = np.multiply(1.0 / beta, c3)
        c12 = np.add(c1, c2)
        c123 = np.add(c12, c3)

        c1 = prepdV(c1, cb_max)
        print("accel_min-c1 = ", np.min(c1))
        c1 = c1 - np.min(c1)
        c12 = prepdV(c12, cb_max)
        print("accel_min-c2 = ", np.min(c12))
        c12 = c12 - np.min(c12)
        c123 = prepdV(c123, cb_max)
        print("accel_min-c3 = ", np.min(c123))
        c123 = c123 - np.min(c123)

        itb = 10
        ite = len(c12) - 5
        dV = c12[itb:ite]
        accel = accel_amd(dV, T)
        print(
            "Time-averaged acceleration factor with Gaussian approximation: ",
            accel,
        )

    elif args.job == "amdweight_MC":
        n = order
        MCweight = np.zeros(len(dV))
        beta_dV = np.multiply(dV, beta)
        for x in range(0, n + 1):
            MCweight = np.add(
                MCweight,
                (
                    np.divide(
                        np.power(beta_dV, x), float(scipy.special.factorial(x))
                    )
                ),
            )
        weights = MCweight
        hist, newedgesX = np.histogram(data, bins=binsX, weights=weights)
        hist = prephist(hist, T, cb_max)
    elif args.job == "amdweight":
        hist, newedgesX = np.histogram(data, bins=binsX, weights=weights)
        hist = prephist(hist, T, cb_max)
    else:
        hist, newedgesX = np.histogram(data, bins=binsX, weights=None)
        hist = prephist(hist, T, cb_max)

    ##SAVE FREE ENERGY DATA INTO A FILE
    if (
        args.job == "amdweight_MC"
        or args.job == "amdweight"
        or args.job == "noweight"
    ):
        pmffile = "pmf-" + str(args.input) + ".xvg"
        return hist, binsX, weights, data
        # output_pmf(pmffile,hist,binsX)

    if args.job == "amdweight_CE":
        hist = pmf_c1
        pmffile = "pmf-c1-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        hist = pmf_c3
        pmffile = "pmf-c3-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        hist = pmf_c2
        pmffile = "pmf-c2-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)
        return hist, binsX, weights, data

    if args.job == "amd_time":
        hist = c1
        pmffile = "accel-c1-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        hist = c12
        pmffile = "accel-c2-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        hist = c123
        pmffile = "accel-c3-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

    ##SAVE WEIGHTS
    if args.job == "amdweight_MC" or args.job == "amdweight":
        pmffile = "weights-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,weights,data)
    if args.job == "amdweight_CE":
        # hist = np.exp(c1)
        pmffile = "weights-c1-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        # hist = np.exp(c12)
        pmffile = "weights-c2-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

        # hist = np.exp(c123)
        pmffile = "weights-c3-" + str(args.input) + ".xvg"
        # output_pmf(pmffile,hist,binsX)

    if args.job == "histo":
        hist, newedgesX = histo(data, hist_min, binsX, discX)
        pmffile = "histo-" + str(args.input) + ".xvg"
        # output_dV_anharm(pmffile,binsX,hist)

    if args.job == "amd_dV":
        hist, newedgesX, binf, dV_avg, dV_std, dV_anharm, dV_mat = reweight_dV(
            data, hist_min, binsX, discX, dV, T
        )

        pmffile = "dV-hist-" + str(args.input) + ".xvg"
        # output_dV(pmffile,dV)

        dV_avg_all = np.average(dV)
        dV_std_all = np.std(dV)
        print("dV all: avg = ", dV_avg_all, "std = ", dV_std_all)

        alpha = anharm(dV)
        print("Anharmonicity of all dV = " + str(alpha))

        pmffile = "dV-anharm-" + str(args.input) + ".xvg"
        # output_dV_anharm(pmffile,binsX,dV_anharm)

        pmffile = "dV-stat-" + str(args.input) + ".xvg"
        # output_dV_stat(pmffile,binsX,dV_avg,dV_std,dV_anharm)

        pmffile = "dV-mat-" + str(args.input) + ".xvg"
        # output_dV_mat(pmffile,binsX,hist,dV_avg,dV_std,dV_anharm,dV_mat)

    ###PLOTTING FUNCTION FOR WEIGHTS histogram
    if plt_figs:
        [hist, edges] = np.histogram(weights, bins=100)
        width = np.absolute(np.subtract(edges[0], edges[1]))
        plt.figure(1, figsize=(11, 8.5))
        plt.bar(edges[:100], hist, width=width, log=True)
        plt.yscale(
            "log"
        )  ###if typerror is thrown delete .matplotlib/fontList.cache  file
        plt.xticks(fontsize="18")
        plt.yticks(fontsize="18")
        # plt.savefig('weights.png',bbox_inches=0)
        print("FIGURE SAVED weights.png")
        plt.show()

    return hist, binsX, weights, data


def reweight_1d(input_1d, weightf, jobtype, slice_mask=None):
    """Reweighting 1D Uses an adapted version of Miao et al py-reweighting
    under the hood.
    returns actual weights
    input: 1d quantity to reweight (np array), weightfile, jobtype
    """
    from .utils import dotdict
    import numpy as np

    # Determine bounds and make 100 steps per dimension
    max_X = input_1d.max()
    min_X = input_1d.min()
    disc = (max_X - min_X) / 100

    # now do reweighting...
    # possible arguments are:
    # args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff, T, Emax, fit, order}
    # -job amdweight_MC -weight weights.dat -Xdim -3 3 -discX 0.1 -Ydim -3 3 -discY 0.1 -Emax 20
    arguments = {
        "job": jobtype,
        "weight": weightf,
        "Xdim": [min_X, max_X],
        "disc": str(disc),
        "Emax": 20,
    }
    arguments = dotdict(arguments)
    hist, binsX, weights, data = pyreweight_1d(input_1d, arguments, slice_mask)

    # c1 = np.digitize(input_1d, binsX) - 1  # shape (50000,)
    # #c2 = np.digitize(reduced_dihedrals[:,1], binsY) - 1 # shape (50000,)
    # print(c1)
    # weights = hist[c1]
    # set weights over threshold to const. value
    # weights[weights > 15] = 15
    return weights


def reweight_1d_pmf(
    input_1d, weightf, jobtype, slice_mask=None, weight_data=None
):
    """Reweighting 1D Uses an adapted version of Miao et al py-reweighting
    under the hood.
    returns 1d pmf

    input: 1d quantity to reweight (np array), weightfile, jobtype,
    bool array to slice traj if necessary for cis/trans
    """
    from .utils import dotdict
    import numpy as np

    # Determine bounds and make 100 steps per dimension
    max_X = input_1d.max()
    min_X = input_1d.min()
    disc = (max_X - min_X) / 100

    # now do reweighting...
    # possible arguments are:
    # args={input: "input file", job, weight, Xdim, Ydim, discX, discY, cutoff, T, Emax, fit, order}
    # -job amdweight_MC -weight weights.dat -Xdim -3 3 -discX 0.1 -Ydim -3 3 -discY 0.1 -Emax 20
    arguments = {
        "job": jobtype,
        "weight": weightf,
        "Xdim": [min_X, max_X],
        "disc": disc,
        "Emax": 20,
    }
    arguments = dotdict(arguments)
    hist, binsX, weights, data = pyreweight_1d(
        input_1d, arguments, slice_mask, weight_data
    )

    # c1 = np.digitize(input_1d, binsX) - 1  # shape (50000,)
    # #c2 = np.digitize(reduced_dihedrals[:,1], binsY) - 1 # shape (50000,)
    # print(c1)
    # weights = hist[c1]
    # set weights over threshold to const. value
    # weights[weights > 15] = 15
    return hist, binsX
