import sklearn
import scipy
import numpy as np


def compute_MAE(exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the mean absolute error

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns MAE, upper CI, lower CI
    """
    MAE = sklearn.metrics.mean_absolute_error(exp, md)

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(exp.to_list(), md.to_list())

        # perform statistical eval
        result = sklearn.metrics.mean_absolute_error(resample[0], resample[1])
        boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = np.percentile(boot_stats, p)
    return MAE, upper, lower


def compute_MSE(dev, CI=0.95, n_bootstrap=10000):
    """Compute the mean signed error

    Parameters
    ----------
    dev : np.array
        deviation between experiment and simulation
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns MSE, upper CI, lower CI
    """
    # MSE - Mean signed error
    MSE = dev.mean()

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(dev)

        # perform statistical eval
        result = resample.mean()
        boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = np.percentile(boot_stats, p)
    # plt.hist(boot_stats)
    return MSE, upper, lower


def compute_RMSD(exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the root mean squared deviation

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns RMSD, upper CI, lower CI
    """
    # RMSD - root mean squared error
    RMSD = sklearn.metrics.mean_squared_error(exp, md, squared=False)

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(exp.to_list(), md.to_list())

        # perform statistical eval
        result = sklearn.metrics.mean_squared_error(
            resample[0], resample[1], squared=False
        )
        boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = np.percentile(boot_stats, p)
    return RMSD, upper, lower


def compute_pearsonr(exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the pearson r correlation coefficient

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns pearsonr, upper CI, lower CI
    """
    # pearsonr
    pearsonr = scipy.stats.pearsonr(exp, md)[0]
    # slope, intercept, pearsonr, p_value, std_err = scipy.stats.
    # linregress(exp.to_list(), md.to_list())

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(exp.to_list(), md.to_list())

        # perform statistical eval
        result = scipy.stats.pearsonr(resample[0], resample[1])
        # _, _2, result, _3, _4 =
        # scipy.stats.linregress(resample[0], resample[1])
        boot_stats.append(result[0])
        # boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = min(1.0, np.percentile(boot_stats, p))
    return pearsonr, upper, lower


def compute_kendalltau(exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the kendall tau

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns kendall tau, upper CI, lower CI
    """
    kendalltau = scipy.stats.kendalltau(exp, md)[0]

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(exp.to_list(), md.to_list())

        # perform statistical eval
        result = scipy.stats.kendalltau(resample[0], resample[1])
        boot_stats.append(result[0])

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = min(1.0, np.percentile(boot_stats, p))
    return kendalltau, upper, lower


def compute_chisquared(exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the chi squared statistic

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns chi squared, upper CI, lower CI
    """
    chi_elements = (md - exp) ** 2 / exp
    chisq = chi_elements.sum()

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(chi_elements)

        # perform statistical eval
        result = resample.sum()
        boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = np.percentile(boot_stats, p)
    return chisq, upper, lower


def compute_fulfilled_percentage(NOE_df):
    """Compute the percentage of fulfilled NOE metric.

    Parameters
    ----------
    NOE_df : pandas.DataFrame
        a NOE dataframe, produced by the noe.get_NOE() function

    Returns
    -------
    _type_
        % NOE fulfilled
    """
    # I: exp_low/exp_high exist
    # exp_low <= md <= exp_high        -> fulfilled
    # md < exp_low OR exp_high < md  -> not fulfilled

    # II: no exp_low:
    # md <= exp_high     -> fulfilled?

    # md > exp_high      -> not fulfilled

    # III: no exp_high:
    # 20% of exp value ~ exp_high
    # md < exp_low       -> not fulfilled

    # IV: no exp_low, no exp_high:
    # 20% of exp value ~ exp high

    # treat non-existence of lower bound as a lower bound of 0
    # (ignore / always fulfilled) in previous scripts
    # so this should work for no lower bound available.
    # if there is only bounds, e.g. no NMR exp, this should still work,
    # since the if np.all evaluates to false,
    # and the 'NMR exp' column is not required..

    # upper bound does not exist:
    if np.all(NOE_df["upper bound"] == 0):
        # set higher bound to 20% of experimental value
        high_bound_value = 0.2
        NOE_df["upper bound"] = NOE_df["NMR exp"] + (
            NOE_df["NMR exp"] * high_bound_value
        )

    fulfilled = (NOE_df["md"] <= NOE_df["upper bound"]) & (
        NOE_df["md"] >= NOE_df["lower bound"]
    )

    return sum(fulfilled) / len(fulfilled)


def compute_RMSD_stepwise(NOE_df, exp, md, CI=0.95, n_bootstrap=10000):
    """Compute the root mean squared deviation

    Parameters
    ----------
    exp : numpy.array
        Experimental values
    md : numpy.array
        simulation values
    CI : float, optional
        Confidence Interval, by default 0.95
    n_bootstrap : int, optional
        number of bootstrap intervals, by default 10000

    Returns
    -------
    (float, float, float)
        Returns RMSD, upper CI, lower CI
    """
    # create deepcopy of md and exp
    md = md.copy()
    exp = exp.copy()

    # upper bound does not exist:
    if np.all(NOE_df["upper bound"] == 0):
        # set higher bound to 20% of experimental value
        high_bound_value = 0.2
        NOE_df["upper bound"] = NOE_df["NMR exp"] + (
            NOE_df["NMR exp"] * high_bound_value
        )

    fulfilled = (NOE_df["md"] <= NOE_df["upper bound"]) & (
        NOE_df["md"] >= NOE_df["lower bound"]
    )

    # set the fulfilled NOEs to the experimental value -> RMSD_{fulfilled} = 0
    md[fulfilled] = exp[fulfilled]

    # set the not fulfilled NOEs to the upper bound -> RMSD only considers violation
    md[~fulfilled] = NOE_df["upper bound"][~fulfilled]

    # compute RMSD 
    RMSD = sklearn.metrics.mean_squared_error(exp, md, squared=False)

    # bootstrap CI's
    # Compute errors/CI:
    boot_stats = []

    for i in range(n_bootstrap):
        resample = sklearn.utils.resample(exp.to_list(), md.to_list())

        # perform statistical eval
        result = sklearn.metrics.mean_squared_error(
            resample[0], resample[1], squared=False
        )
        boot_stats.append(result)

    # Confidence intervals

    p = ((1.0 - CI) / 2.0) * 100
    lower = max(0.0, np.percentile(boot_stats, p))
    p = (CI + ((1.0 - CI) / 2.0)) * 100
    upper = np.percentile(boot_stats, p)
    return RMSD, upper, lower